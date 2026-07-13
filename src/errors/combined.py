"""Combined Rydberg CZ error-budget simulation.

This module is the integration point for the five isolated error channels.  It
uses a 9-state basis that keeps the dark qubit state, the laser-coupled ground
state, both single-Rydberg states, and the doubly excited ``|rr>`` state.  That
one basis is large enough for finite blockade, Rydberg decay, intermediate-state
scattering, Doppler detuning, and shot-to-shot Rabi-amplitude noise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from qutip import Qobj, mesolve

from ..analytical import epsilon_total_additive
from ..fidelity import CZ_TARGET, LOCAL_Z_PRODUCT, pedersen_fidelity, pedersen_fidelity_kraus
from ..params import (
    DEFAULT_OMEGA_RAD_PER_US,
    DEFAULT_TEMPERATURE_K,
    K_EFF_RAD_PER_UM,
    RB87_MASS_KG,
    get_rydberg_params,
    thermal_velocity_rms_um_per_us,
)
from .amplitude import DEFAULT_SIGMA_OMEGA, run_amplitude_noise_gate
from .blockade import run_blockade_gate
from .decay import _choi_from_process_outputs, _kraus_from_choi, run_decay_gate
from .doppler import run_doppler_gate_mc
from .scattering import run_scattering_gate

COMBINED_BASIS: Final[tuple[str, ...]] = ("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg", "rr")
"""Two-atom basis including dark states and the doubly excited Rydberg state."""

COMPUTATIONAL_BASIS: Final[tuple[str, ...]] = ("00", "01", "10", "11")
COMPUTATIONAL_TO_COMBINED: Final[dict[str, str]] = {
    "00": "00",
    "01": "0g",
    "10": "g0",
    "11": "gg",
}

BASIS_INDEX: Final[dict[str, int]] = {label: index for index, label in enumerate(COMBINED_BASIS)}
COMPUTATIONAL_INDICES: Final[tuple[int, ...]] = tuple(
    BASIS_INDEX[COMPUTATIONAL_TO_COMBINED[label]] for label in COMPUTATIONAL_BASIS
)

DEFAULT_INTERMEDIATE_DETUNING_MHZ: Final[float] = 1000.0
DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US: Final[float] = 2.0 * np.pi * DEFAULT_INTERMEDIATE_DETUNING_MHZ
DEFAULT_COMBINED_SAMPLES: Final[int] = 24

DEFAULT_MESOLVE_OPTIONS: Final[dict[str, float | int | bool | str]] = {
    "store_states": True,
    "nsteps": 100_000,
    "atol": 1e-11,
    "rtol": 1e-11,
    "method": "bdf",
}


@dataclass(frozen=True)
class CombinedGateResult:
    """Monte Carlo result for all modeled errors switched on together."""

    omega: float
    gamma: float
    blockade_shift: float
    gamma_e: float
    delta_p: float
    sigma_omega: float
    k_eff_rad_per_um: float
    temperature_K: float
    mass_kg: float
    omega1_over_omega2: float
    n_samples: int
    seed: int | None
    average_fidelity: float
    std_fidelity: float
    per_shot_fidelities: np.ndarray
    rabi_scales: np.ndarray
    detunings_rad_per_us: np.ndarray
    analytical_errors: dict[str, float]
    additive_sum: float

    @property
    def average_gate_fidelity(self) -> float:
        return self.average_fidelity

    @property
    def infidelity(self) -> float:
        return 1.0 - self.average_fidelity

    @property
    def standard_error(self) -> float:
        if self.n_samples <= 1:
            return 0.0
        return self.std_fidelity / np.sqrt(self.n_samples)


@dataclass(frozen=True)
class ErrorBudgetResult:
    """Numerical and analytical error budget at one operating point."""

    combined: CombinedGateResult
    individual_errors: dict[str, float]
    individual_analytical_errors: dict[str, float]
    additive_numerical_sum: float
    additive_analytical_sum: float

    @property
    def combined_numerical_error(self) -> float:
        return self.combined.infidelity

    @property
    def additive_relative_difference(self) -> float:
        numerical = self.combined_numerical_error
        if numerical == 0.0:
            return 0.0 if self.additive_analytical_sum == 0.0 else np.inf
        return abs(self.additive_analytical_sum - numerical) / numerical


def _check_positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
    return value


def _check_nonnegative_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be non-negative and finite, got {value!r}")
    return value


def _check_n_samples(n_samples: int) -> int:
    if not isinstance(n_samples, int) or n_samples < 1:
        raise ValueError("n_samples must be an integer >= 1")
    return n_samples


def _check_blockade_shift(blockade_shift: float) -> float:
    blockade_shift = float(blockade_shift)
    if np.isposinf(blockade_shift):
        return blockade_shift
    if not np.isfinite(blockade_shift) or blockade_shift <= 0.0:
        raise ValueError(f"blockade_shift must be positive finite or +inf, got {blockade_shift!r}")
    return blockade_shift


def _qobj(data: np.ndarray) -> Qobj:
    return Qobj(np.asarray(data, dtype=complex))


def _scattering_rate_for_actual_omega(
    omega_actual: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float,
) -> float:
    """Return the driven-pulse scattering dephasing rate for this shot."""

    if gamma_e == 0.0 or omega_actual == 0.0:
        return 0.0
    beam_imbalance = 0.5 * (omega1_over_omega2 + 1.0 / omega1_over_omega2)
    return float(gamma_e * abs(omega_actual) * beam_imbalance / delta_p)


def _hamiltonian(
    omega_actual: float,
    blockade_shift: float,
    atom_index: int,
    detuning: float,
) -> Qobj:
    """Return the 9-state Hamiltonian for one driven pulse segment."""

    if not np.isfinite(omega_actual):
        raise ValueError(f"omega_actual must be finite, got {omega_actual!r}")
    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")
    if not np.isfinite(detuning):
        raise ValueError(f"detuning must be finite, got {detuning!r}")

    data = np.zeros((len(COMBINED_BASIS), len(COMBINED_BASIS)), dtype=complex)
    finite_blockade = np.isfinite(blockade_shift)

    if atom_index == 1:
        couplings = [("g0", "r0"), ("gg", "rg")]
        if finite_blockade:
            couplings.append(("gr", "rr"))
        detuned_labels = ("r0", "rg", "rr")
    else:
        couplings = [("0g", "0r"), ("gg", "gr")]
        if finite_blockade:
            couplings.append(("rg", "rr"))
        detuned_labels = ("0r", "gr", "rr")

    for left, right in couplings:
        i = BASIS_INDEX[left]
        j = BASIS_INDEX[right]
        data[i, j] = data[j, i] = omega_actual / 2.0
    for label in detuned_labels:
        data[BASIS_INDEX[label], BASIS_INDEX[label]] += -float(detuning)
    if finite_blockade:
        data[BASIS_INDEX["rr"], BASIS_INDEX["rr"]] += blockade_shift
    return _qobj(data)


def _decay_collapse_operators(gamma: float) -> list[Qobj]:
    """Return atom-resolved ``sqrt(gamma)|g><r|`` operators in COMBINED_BASIS."""

    if gamma == 0.0:
        return []

    amplitude = np.sqrt(gamma)
    size = len(COMBINED_BASIS)
    atom1 = np.zeros((size, size), dtype=complex)
    atom2 = np.zeros((size, size), dtype=complex)

    for destination, source in (("g0", "r0"), ("gg", "rg"), ("gr", "rr")):
        atom1[BASIS_INDEX[destination], BASIS_INDEX[source]] = amplitude
    for destination, source in (("0g", "0r"), ("gg", "gr"), ("rg", "rr")):
        atom2[BASIS_INDEX[destination], BASIS_INDEX[source]] = amplitude

    return [_qobj(atom1), _qobj(atom2)]


def _scattering_collapse_operators(rate: float, atom_index: int) -> list[Qobj]:
    """Return driven-atom ``sqrt(rate)|r><r|`` dephasing operators."""

    if rate == 0.0:
        return []
    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")

    size = len(COMBINED_BASIS)
    data = np.zeros((size, size), dtype=complex)
    labels = ("r0", "rg", "rr") if atom_index == 1 else ("0r", "gr", "rr")
    amplitude = np.sqrt(rate)
    for label in labels:
        data[BASIS_INDEX[label], BASIS_INDEX[label]] = amplitude
    return [_qobj(data)]


def _embed_computational_operator(operator: np.ndarray) -> Qobj:
    matrix = np.asarray(operator, dtype=complex)
    if matrix.shape != (4, 4):
        raise ValueError(f"computational operator must be 4x4, got {matrix.shape}")

    embedded = np.zeros((len(COMBINED_BASIS), len(COMBINED_BASIS)), dtype=complex)
    for row, combined_row in enumerate(COMPUTATIONAL_INDICES):
        for col, combined_col in enumerate(COMPUTATIONAL_INDICES):
            embedded[combined_row, combined_col] = matrix[row, col]
    return _qobj(embedded)


def _project_to_computational(operator: Qobj) -> np.ndarray:
    full = np.asarray(operator.full(), dtype=complex)
    return np.asarray(full[np.ix_(COMPUTATIONAL_INDICES, COMPUTATIONAL_INDICES)], dtype=complex)


def _mesolve_final(hamiltonian: Qobj, rho0: Qobj, duration: float, c_ops: list[Qobj]) -> Qobj:
    result = mesolve(
        hamiltonian,
        rho0,
        np.array([0.0, duration]),
        c_ops=c_ops,
        e_ops=[],
        options=DEFAULT_MESOLVE_OPTIONS,
    )
    return result.states[-1]


def _run_three_pulse_sequence(
    rho0: Qobj,
    *,
    omega_nominal: float,
    rabi_scale: float,
    blockade_shift: float,
    gamma: float,
    scattering_rate: float,
    delta1: float,
    delta2: float,
) -> Qobj:
    t_pi = np.pi / omega_nominal
    omega_actual = omega_nominal * rabi_scale
    decay_ops = _decay_collapse_operators(gamma)

    h1 = _hamiltonian(omega_actual, blockade_shift, atom_index=1, detuning=delta1)
    h2 = _hamiltonian(omega_actual, blockade_shift, atom_index=2, detuning=delta2)
    c_ops_atom1 = decay_ops + _scattering_collapse_operators(scattering_rate, atom_index=1)
    c_ops_atom2 = decay_ops + _scattering_collapse_operators(scattering_rate, atom_index=2)

    rho = _mesolve_final(h1, rho0, t_pi, c_ops_atom1)
    rho = _mesolve_final(h2, rho, 2.0 * t_pi, c_ops_atom2)
    rho = _mesolve_final(h1, rho, t_pi, c_ops_atom1)
    return rho


def _propagate_process_outputs_for_sample(
    *,
    omega: float,
    rabi_scale: float,
    blockade_shift: float,
    gamma: float,
    scattering_rate: float,
    delta1: float,
    delta2: float,
) -> dict[tuple[int, int], np.ndarray]:
    outputs: dict[tuple[int, int], np.ndarray] = {}
    for row in range(4):
        for col in range(4):
            basis_operator = np.zeros((4, 4), dtype=complex)
            basis_operator[row, col] = 1.0
            rho0 = _embed_computational_operator(basis_operator)
            outputs[(row, col)] = _project_to_computational(
                _run_three_pulse_sequence(
                    rho0,
                    omega_nominal=omega,
                    rabi_scale=rabi_scale,
                    blockade_shift=blockade_shift,
                    gamma=gamma,
                    scattering_rate=scattering_rate,
                    delta1=delta1,
                    delta2=delta2,
                )
            )
    return outputs


def _fidelity_for_sample(
    *,
    omega: float,
    rabi_scale: float,
    blockade_shift: float,
    gamma: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float,
    delta1: float,
    delta2: float,
) -> float:
    omega_actual = omega * rabi_scale
    scattering_rate = _scattering_rate_for_actual_omega(omega_actual, gamma_e, delta_p, omega1_over_omega2)
    outputs = _propagate_process_outputs_for_sample(
        omega=omega,
        rabi_scale=rabi_scale,
        blockade_shift=blockade_shift,
        gamma=gamma,
        scattering_rate=scattering_rate,
        delta1=delta1,
        delta2=delta2,
    )
    choi = _choi_from_process_outputs(outputs)
    kraus_ops = _kraus_from_choi(choi)
    corrected_kraus_ops = tuple(LOCAL_Z_PRODUCT @ op for op in kraus_ops)
    return pedersen_fidelity_kraus(corrected_kraus_ops, CZ_TARGET)


def _sample_classical_errors(
    *,
    sigma_omega: float,
    k_eff_rad_per_um: float,
    temperature_K: float,
    mass_kg: float,
    n_samples: int,
    seed: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    if sigma_omega == 0.0:
        rabi_scales = np.ones(n_samples, dtype=float)
    else:
        rabi_scales = 1.0 + rng.normal(loc=0.0, scale=sigma_omega, size=n_samples)

    velocity_rms = thermal_velocity_rms_um_per_us(temperature_K, mass_kg)
    if k_eff_rad_per_um == 0.0 or velocity_rms == 0.0:
        detunings = np.zeros((n_samples, 2), dtype=float)
    else:
        velocities = rng.normal(loc=0.0, scale=velocity_rms, size=(n_samples, 2))
        detunings = k_eff_rad_per_um * velocities
    return rabi_scales, detunings


def _default_parameters() -> tuple[float, float, float, float]:
    params = get_rydberg_params()
    return (
        params.omega_rad_per_us,
        params.rydberg_decay_rate_per_us,
        params.blockade_shift_rad_per_us,
        params.intermediate_decay_rate_per_us,
    )


def run_combined_gate(
    omega: float | None = None,
    gamma: float | None = None,
    blockade_shift: float | None = None,
    gamma_e: float | None = None,
    delta_p: float = DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
    sigma_omega: float = DEFAULT_SIGMA_OMEGA,
    k_eff_rad_per_um: float = K_EFF_RAD_PER_UM,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    mass_kg: float = RB87_MASS_KG,
    *,
    omega1_over_omega2: float = 1.0,
    n_samples: int = DEFAULT_COMBINED_SAMPLES,
    seed: int | None = 2024,
) -> CombinedGateResult:
    """Run the full combined-error simulation at one operating point."""

    default_omega, default_gamma, default_blockade, default_gamma_e = _default_parameters()
    omega = default_omega if omega is None else omega
    gamma = default_gamma if gamma is None else gamma
    blockade_shift = default_blockade if blockade_shift is None else blockade_shift
    gamma_e = default_gamma_e if gamma_e is None else gamma_e

    omega = _check_positive_finite("omega", omega)
    gamma = _check_nonnegative_finite("gamma", gamma)
    blockade_shift = _check_blockade_shift(blockade_shift)
    gamma_e = _check_nonnegative_finite("gamma_e", gamma_e)
    delta_p = _check_positive_finite("delta_p", delta_p)
    sigma_omega = _check_nonnegative_finite("sigma_omega", sigma_omega)
    k_eff_rad_per_um = _check_nonnegative_finite("k_eff_rad_per_um", k_eff_rad_per_um)
    temperature_K = _check_nonnegative_finite("temperature_K", temperature_K)
    mass_kg = _check_positive_finite("mass_kg", mass_kg)
    omega1_over_omega2 = _check_positive_finite("omega1_over_omega2", omega1_over_omega2)
    n_samples = _check_n_samples(n_samples)

    additive_sum, analytical_errors = epsilon_total_additive(
        omega=omega,
        gamma=gamma,
        blockade_shift=blockade_shift,
        gamma_e=gamma_e,
        delta_p=delta_p,
        sigma_omega=sigma_omega,
        k_eff_rad_per_um=k_eff_rad_per_um,
        temperature_K=temperature_K,
        mass_kg=mass_kg,
        omega1_over_omega2=omega1_over_omega2,
    )

    rabi_scales, detunings = _sample_classical_errors(
        sigma_omega=sigma_omega,
        k_eff_rad_per_um=k_eff_rad_per_um,
        temperature_K=temperature_K,
        mass_kg=mass_kg,
        n_samples=n_samples,
        seed=seed,
    )

    fidelities = np.empty(n_samples, dtype=float)
    for index, rabi_scale in enumerate(rabi_scales):
        fidelities[index] = _fidelity_for_sample(
            omega=omega,
            rabi_scale=float(rabi_scale),
            blockade_shift=blockade_shift,
            gamma=gamma,
            gamma_e=gamma_e,
            delta_p=delta_p,
            omega1_over_omega2=omega1_over_omega2,
            delta1=float(detunings[index, 0]),
            delta2=float(detunings[index, 1]),
        )

    return CombinedGateResult(
        omega=omega,
        gamma=gamma,
        blockade_shift=blockade_shift,
        gamma_e=gamma_e,
        delta_p=delta_p,
        sigma_omega=sigma_omega,
        k_eff_rad_per_um=k_eff_rad_per_um,
        temperature_K=temperature_K,
        mass_kg=mass_kg,
        omega1_over_omega2=omega1_over_omega2,
        n_samples=n_samples,
        seed=seed,
        average_fidelity=float(np.mean(fidelities)),
        std_fidelity=float(np.std(fidelities, ddof=1)) if n_samples > 1 else 0.0,
        per_shot_fidelities=fidelities,
        rabi_scales=rabi_scales,
        detunings_rad_per_us=detunings,
        analytical_errors=analytical_errors,
        additive_sum=additive_sum,
    )


def evaluate_error_budget(
    *,
    n_samples: int = DEFAULT_COMBINED_SAMPLES,
    seed: int | None = 2024,
    individual_n_samples: int = 200,
) -> ErrorBudgetResult:
    """Evaluate individual channels and the full combined channel at baseline."""

    params = get_rydberg_params()
    combined = run_combined_gate(n_samples=n_samples, seed=seed)

    decay = run_decay_gate(params.omega_rad_per_us, params.rydberg_decay_rate_per_us).infidelity
    blockade_result = run_blockade_gate(params.omega_rad_per_us, params.blockade_shift_rad_per_us, n_steps_per_pi=160)
    blockade = 1.0 - pedersen_fidelity(blockade_result["unitary"], CZ_TARGET, correct_local_z=True)
    doppler = run_doppler_gate_mc(
        omega=params.omega_rad_per_us,
        k_eff_rad_per_um=K_EFF_RAD_PER_UM,
        temperature_K=DEFAULT_TEMPERATURE_K,
        mass_kg=RB87_MASS_KG,
        n_samples=individual_n_samples,
        seed=None if seed is None else seed + 101,
    ).infidelity
    scattering = run_scattering_gate(
        omega=params.omega_rad_per_us,
        gamma_e=params.intermediate_decay_rate_per_us,
        delta_p=DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
    ).infidelity
    amplitude = run_amplitude_noise_gate(
        omega=params.omega_rad_per_us,
        sigma_omega=DEFAULT_SIGMA_OMEGA,
        n_samples=individual_n_samples,
        seed=None if seed is None else seed + 202,
    ).infidelity

    individual_errors = {
        "decay": float(decay),
        "blockade": float(blockade),
        "doppler": float(doppler),
        "scattering": float(scattering),
        "amplitude": float(amplitude),
    }
    return ErrorBudgetResult(
        combined=combined,
        individual_errors=individual_errors,
        individual_analytical_errors=combined.analytical_errors,
        additive_numerical_sum=float(sum(individual_errors.values())),
        additive_analytical_sum=combined.additive_sum,
    )
