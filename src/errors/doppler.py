"""Monte Carlo Doppler dephasing simulation for the pi-2pi-pi CZ gate.

Doppler shifts are quasi-static classical detunings during one gate shot.  Each
trial samples independent axial velocities for the two atoms, runs the coherent
three-pulse sequence with those detunings, and averages the projected CZ
fidelity over the velocity ensemble.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from qutip import Qobj, basis, sesolve

from ..fidelity import CZ_TARGET, pedersen_fidelity
from ..hamiltonian import build_ideal_hamiltonian, single_atom_hamiltonian
from ..params import (
    DEFAULT_OMEGA_RAD_PER_US,
    DEFAULT_TEMPERATURE_K,
    K_EFF_RAD_PER_UM,
    RB87_MASS_KG,
    thermal_velocity_rms_um_per_us,
)

DEFAULT_SESOLVE_OPTIONS: Final[dict[str, float | int | bool]] = {
    "store_states": True,
    "nsteps": 10_000,
    "atol": 1e-11,
    "rtol": 1e-11,
}


@dataclass(frozen=True)
class DopplerGateResult:
    """Numerical result for isolated Doppler dephasing."""

    omega: float
    k_eff_rad_per_um: float
    temperature_K: float
    mass_kg: float
    n_samples: int
    seed: int | None
    t_pi: float
    velocity_rms_um_per_us: float
    detuning_rms_rad_per_us: float
    average_fidelity: float
    std_fidelity: float
    per_shot_fidelities: np.ndarray

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


def projected_gate_for_detunings(omega: float, delta1: float, delta2: float) -> np.ndarray:
    """Return the computational-space gate for one Doppler realization.

    ``delta1`` and ``delta2`` are atom-resolved Doppler detunings in rad/us.
    The returned 4x4 matrix is the computational projection of the coherent
    evolution, so diagonal magnitudes below one represent leakage left in
    Rydberg states after imperfect detuned pulses.
    """

    omega = _check_positive_finite("omega", omega)
    if not np.isfinite(delta1) or not np.isfinite(delta2):
        raise ValueError("detunings must be finite")

    t_pi = np.pi / omega

    psi_g = basis(2, 0)
    h_atom1_single = single_atom_hamiltonian(omega, detuning=delta1)
    h_atom2_single = single_atom_hamiltonian(omega, detuning=delta2)

    # |01>: atom 2 performs only the middle 2pi pulse.
    state_01 = _solve_final(h_atom2_single, psi_g, 2.0 * t_pi)
    amp_01 = _amplitude(state_01, 0)

    # |10>: atom 1 receives the first and third pi pulses; the middle segment is
    # a laser-off hold and contributes no Doppler detuning in this rotating frame.
    state_10 = _solve_final(h_atom1_single, psi_g, t_pi)
    state_10 = _solve_final(h_atom1_single, state_10, t_pi)
    amp_10 = _amplitude(state_10, 0)

    # |11>: evolve in the infinite-blockade basis |gg>, |gr>, |rg>.  The atom-2
    # pulse is blocked for the |rg> branch but still acts on any residual |gg>
    # amplitude created by a detuned first pi pulse.
    psi_gg = basis(3, 0)
    h_atom1_pair = build_ideal_hamiltonian(omega, atom_index=1, detuning=delta1)
    h_atom2_pair = build_ideal_hamiltonian(omega, atom_index=2, detuning=delta2)
    state_11 = _solve_final(h_atom1_pair, psi_gg, t_pi)
    state_11 = _solve_final(h_atom2_pair, state_11, 2.0 * t_pi)
    state_11 = _solve_final(h_atom1_pair, state_11, t_pi)
    amp_11 = _amplitude(state_11, 0)

    return np.diag(np.array([1.0 + 0.0j, amp_01, amp_10, amp_11], dtype=complex))


def _zero_doppler_result(
    *,
    omega: float,
    k_eff_rad_per_um: float,
    temperature_K: float,
    mass_kg: float,
    n_samples: int,
    seed: int | None,
) -> DopplerGateResult:
    fidelities = np.ones(n_samples, dtype=float)
    velocity_rms = thermal_velocity_rms_um_per_us(temperature_K, mass_kg)
    return DopplerGateResult(
        omega=omega,
        k_eff_rad_per_um=k_eff_rad_per_um,
        temperature_K=temperature_K,
        mass_kg=mass_kg,
        n_samples=n_samples,
        seed=seed,
        t_pi=np.pi / omega,
        velocity_rms_um_per_us=velocity_rms,
        detuning_rms_rad_per_us=k_eff_rad_per_um * velocity_rms,
        average_fidelity=1.0,
        std_fidelity=0.0,
        per_shot_fidelities=fidelities,
    )


def run_doppler_gate_mc(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    k_eff_rad_per_um: float = K_EFF_RAD_PER_UM,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    mass_kg: float = RB87_MASS_KG,
    *,
    n_samples: int = 500,
    seed: int | None = None,
) -> DopplerGateResult:
    """Run a Monte Carlo average over thermal Doppler detunings.

    Velocities are sampled from the 1D Maxwell-Boltzmann distribution
    ``N(0, sqrt(k_B T / m))``.  With velocities in um/us and ``k_eff`` in
    rad/um, each shot uses ``delta = k_eff * v`` in rad/us.
    """

    omega = _check_positive_finite("omega", omega)
    k_eff_rad_per_um = _check_nonnegative_finite("k_eff_rad_per_um", k_eff_rad_per_um)
    temperature_K = _check_nonnegative_finite("temperature_K", temperature_K)
    mass_kg = _check_positive_finite("mass_kg", mass_kg)
    n_samples = _check_n_samples(n_samples)

    if temperature_K == 0.0 or k_eff_rad_per_um == 0.0:
        return _zero_doppler_result(
            omega=omega,
            k_eff_rad_per_um=k_eff_rad_per_um,
            temperature_K=temperature_K,
            mass_kg=mass_kg,
            n_samples=n_samples,
            seed=seed,
        )

    velocity_rms = thermal_velocity_rms_um_per_us(temperature_K, mass_kg)
    rng = np.random.default_rng(seed)
    velocities = rng.normal(loc=0.0, scale=velocity_rms, size=(n_samples, 2))

    fidelities = np.empty(n_samples, dtype=float)
    for index, (v1, v2) in enumerate(velocities):
        gate = projected_gate_for_detunings(
            omega,
            delta1=k_eff_rad_per_um * float(v1),
            delta2=k_eff_rad_per_um * float(v2),
        )
        fidelities[index] = pedersen_fidelity(gate, CZ_TARGET, correct_local_z=True)

    return DopplerGateResult(
        omega=omega,
        k_eff_rad_per_um=k_eff_rad_per_um,
        temperature_K=temperature_K,
        mass_kg=mass_kg,
        n_samples=n_samples,
        seed=seed,
        t_pi=np.pi / omega,
        velocity_rms_um_per_us=velocity_rms,
        detuning_rms_rad_per_us=k_eff_rad_per_um * velocity_rms,
        average_fidelity=float(np.mean(fidelities)),
        std_fidelity=float(np.std(fidelities, ddof=1)) if n_samples > 1 else 0.0,
        per_shot_fidelities=fidelities,
    )
