"""Shot-to-shot Rabi-amplitude-noise simulation for the pi-2pi-pi CZ gate.

Amplitude noise is modeled as a quasi-static fractional error for one gate shot:
``Omega_actual = Omega_nominal * (1 + epsilon)`` with
``epsilon ~ N(0, sigma_omega^2)``.  The pulse durations stay calibrated to the
nominal ``Omega``.  That distinction matters: changing both the Hamiltonian and
pulse durations would hide the pulse-area error entirely.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from qutip import Qobj, basis, sesolve

from ..fidelity import CZ_TARGET, LOCAL_Z_PRODUCT, pedersen_fidelity
from ..params import DEFAULT_OMEGA_RAD_PER_US

DEFAULT_SIGMA_OMEGA: Final[float] = 0.02
"""Evered-like fractional Rabi-frequency noise used as the reference point."""

DEFAULT_N_SAMPLES: Final[int] = 500

DEFAULT_SESOLVE_OPTIONS: Final[dict[str, float | int | bool]] = {
    "store_states": True,
    "nsteps": 10_000,
    "atol": 1e-11,
    "rtol": 1e-11,
}


@dataclass(frozen=True)
class AmplitudeNoiseResult:
    """Monte Carlo result for isolated Rabi-amplitude noise."""

    omega: float
    sigma_omega: float
    n_samples: int
    seed: int | None
    t_pi: float
    average_fidelity: float
    std_fidelity: float
    per_shot_fidelities: np.ndarray
    fractional_rabi_errors: np.ndarray

    @property
    def average_gate_fidelity(self) -> float:
        """Alias matching the other error-channel result classes."""

        return self.average_fidelity

    @property
    def infidelity(self) -> float:
        return 1.0 - self.average_fidelity

    @property
    def standard_error(self) -> float:
        if self.n_samples <= 1:
            return 0.0
        return self.std_fidelity / np.sqrt(self.n_samples)

    @property
    def fidelity_std(self) -> float:
        """Compatibility alias for the Monte Carlo standard deviation."""

        return self.std_fidelity

    @property
    def fidelity_stderr(self) -> float:
        """Compatibility alias for the Monte Carlo standard error."""

        return self.standard_error

    @property
    def sample_fidelities(self) -> np.ndarray:
        """Compatibility alias for per-shot fidelities."""

        return self.per_shot_fidelities


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


def _check_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value):
        raise ValueError(f"{name} must be finite, got {value!r}")
    return value


def _check_n_samples(n_samples: int) -> int:
    if not isinstance(n_samples, int) or n_samples < 1:
        raise ValueError("n_samples must be an integer >= 1")
    return n_samples


def _qobj(data: np.ndarray) -> Qobj:
    return Qobj(np.asarray(data, dtype=complex))


def _single_atom_hamiltonian_signed(omega_actual: float) -> Qobj:
    """Return a two-level Hamiltonian allowing signed Rabi amplitudes."""

    omega_actual = _check_finite("omega_actual", omega_actual)
    data = np.array(
        [[0.0, omega_actual / 2.0], [omega_actual / 2.0, 0.0]],
        dtype=complex,
    )
    return _qobj(data)


def _ideal_hamiltonian_signed(omega_actual: float, atom_index: int) -> Qobj:
    """Return the infinite-blockade Hamiltonian allowing signed Rabi amplitudes."""

    omega_actual = _check_finite("omega_actual", omega_actual)
    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")

    data = np.zeros((3, 3), dtype=complex)
    if atom_index == 1:
        data[0, 2] = data[2, 0] = omega_actual / 2.0
    else:
        data[0, 1] = data[1, 0] = omega_actual / 2.0
    return _qobj(data)


def _solve_final(hamiltonian: Qobj, psi0: Qobj, duration: float) -> Qobj:
    result = sesolve(
        hamiltonian,
        psi0,
        np.array([0.0, duration]),
        e_ops=[],
        options=DEFAULT_SESOLVE_OPTIONS,
    )
    return result.states[-1]


def _amplitude(state: Qobj, index: int) -> complex:
    return complex(state.full()[index, 0])


def projected_gate_for_rabi_scale(omega: float, rabi_scale: float) -> np.ndarray:
    """Return the projected computational gate for one amplitude-noise shot.

    ``omega`` is the nominal angular Rabi frequency that sets the calibrated
    pulse durations.  ``rabi_scale`` multiplies the Hamiltonian only.  Values
    below zero are allowed: a negative Rabi amplitude is a pi phase shift of the
    drive and still gives a well-defined coherent evolution.
    """

    omega = _check_positive_finite("omega", omega)
    rabi_scale = _check_finite("rabi_scale", rabi_scale)
    omega_actual = omega * rabi_scale
    t_pi = np.pi / omega

    psi_g = basis(2, 0)
    h_single = _single_atom_hamiltonian_signed(omega_actual)

    # |01>: atom 2 receives only the nominal 2pi pulse.
    state_01 = _solve_final(h_single, psi_g, 2.0 * t_pi)
    amp_01 = _amplitude(state_01, 0)

    # |10>: atom 1 receives two nominal pi pulses separated by a laser-off hold.
    state_10 = _solve_final(h_single, psi_g, t_pi)
    state_10 = _solve_final(h_single, state_10, t_pi)
    amp_10 = _amplitude(state_10, 0)

    # |11>: infinite blockade basis |gg>, |gr>, |rg>.  The middle pulse is
    # blocked for |rg>, but any residual |gg> amplitude from an imperfect first
    # pulse still rotates under the atom-2 drive.
    psi_gg = basis(3, 0)
    h_atom1 = _ideal_hamiltonian_signed(omega_actual, atom_index=1)
    h_atom2 = _ideal_hamiltonian_signed(omega_actual, atom_index=2)
    state_11 = _solve_final(h_atom1, psi_gg, t_pi)
    state_11 = _solve_final(h_atom2, state_11, 2.0 * t_pi)
    state_11 = _solve_final(h_atom1, state_11, t_pi)
    amp_11 = _amplitude(state_11, 0)

    return np.diag(np.array([1.0 + 0.0j, amp_01, amp_10, amp_11], dtype=complex))


def _zero_amplitude_result(
    *,
    omega: float,
    sigma_omega: float,
    n_samples: int,
    seed: int | None,
) -> AmplitudeNoiseResult:
    fidelities = np.ones(n_samples, dtype=float)
    fractional_errors = np.zeros(n_samples, dtype=float)
    return AmplitudeNoiseResult(
        omega=omega,
        sigma_omega=sigma_omega,
        n_samples=n_samples,
        seed=seed,
        t_pi=np.pi / omega,
        average_fidelity=1.0,
        std_fidelity=0.0,
        per_shot_fidelities=fidelities,
        fractional_rabi_errors=fractional_errors,
    )


def run_amplitude_noise_gate(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    sigma_omega: float = DEFAULT_SIGMA_OMEGA,
    *,
    n_samples: int = DEFAULT_N_SAMPLES,
    seed: int | None = None,
) -> AmplitudeNoiseResult:
    """Run a Monte Carlo average over quasi-static Rabi-amplitude noise.

    One fractional error is sampled per shot and held fixed across all three
    pulses, matching global shot-to-shot laser-intensity noise.  The known raw
    local-Z phases of the pi-2pi-pi sequence are corrected with the fixed
    ``Z_1 \\otimes Z_2`` frame before evaluating Pedersen fidelity.
    """

    omega = _check_positive_finite("omega", omega)
    sigma_omega = _check_nonnegative_finite("sigma_omega", sigma_omega)
    n_samples = _check_n_samples(n_samples)

    if sigma_omega == 0.0:
        return _zero_amplitude_result(
            omega=omega,
            sigma_omega=sigma_omega,
            n_samples=n_samples,
            seed=seed,
        )

    rng = np.random.default_rng(seed)
    fractional_errors = rng.normal(loc=0.0, scale=sigma_omega, size=n_samples)
    fidelities = np.empty(n_samples, dtype=float)

    for index, fractional_error in enumerate(fractional_errors):
        gate = projected_gate_for_rabi_scale(
            omega,
            rabi_scale=1.0 + float(fractional_error),
        )
        fidelities[index] = pedersen_fidelity(LOCAL_Z_PRODUCT @ gate, CZ_TARGET)

    return AmplitudeNoiseResult(
        omega=omega,
        sigma_omega=sigma_omega,
        n_samples=n_samples,
        seed=seed,
        t_pi=np.pi / omega,
        average_fidelity=float(np.mean(fidelities)),
        std_fidelity=float(np.std(fidelities, ddof=1)) if n_samples > 1 else 0.0,
        per_shot_fidelities=fidelities,
        fractional_rabi_errors=fractional_errors,
    )
