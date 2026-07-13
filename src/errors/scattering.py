"""Intermediate-state scattering simulation for the pi-2pi-pi CZ gate.

Rubidium two-photon excitation only visits the 5P intermediate state virtually,
but that virtual population can still scatter photons.  In the large-detuning
limit this is an effective pure-dephasing channel on the atom being driven.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from qutip import Qobj

from ..fidelity import CZ_TARGET, LOCAL_Z_PRODUCT, pedersen_fidelity_kraus
from ..params import DEFAULT_OMEGA_RAD_PER_US
from .decay import (
    COMPUTATIONAL_BASIS,
    DECAY_BASIS,
    _BASIS_INDEX,
    _check_gamma,
    _check_omega,
    _choi_from_process_outputs,
    _embed_computational_operator,
    _hamiltonian,
    _kraus_from_choi,
    _mesolve_final,
    _project_to_computational,
    _qobj,
    run_decay_gate,
)


@dataclass(frozen=True)
class ScatteringGateResult:
    """Numerical result for the isolated intermediate-state-scattering channel."""

    omega: float
    gamma_e: float
    delta_p: float
    omega1_over_omega2: float
    omega1: float
    omega2: float
    scattering_rate_per_us: float
    t_pi: float
    average_gate_fidelity: float
    raw_average_gate_fidelity: float
    kraus_ops: tuple[np.ndarray, ...]
    corrected_kraus_ops: tuple[np.ndarray, ...]
    choi_matrix: np.ndarray
    final_density_matrices: dict[str, np.ndarray]

    @property
    def infidelity(self) -> float:
        return 1.0 - self.average_gate_fidelity


_COMPUTATIONAL_DIM: Final[int] = 4


def _check_delta_p(delta_p: float) -> float:
    delta_p = float(delta_p)
    if not np.isfinite(delta_p) or delta_p <= 0.0:
        raise ValueError(f"delta_p must be positive and finite, got {delta_p!r}")
    return delta_p


def _check_beam_ratio(omega1_over_omega2: float) -> float:
    omega1_over_omega2 = float(omega1_over_omega2)
    if not np.isfinite(omega1_over_omega2) or omega1_over_omega2 <= 0.0:
        raise ValueError(
            "omega1_over_omega2 must be positive and finite, "
            f"got {omega1_over_omega2!r}"
        )
    return omega1_over_omega2


def single_photon_rabi_frequencies(
    omega: float,
    delta_p: float,
    omega1_over_omega2: float = 1.0,
) -> tuple[float, float]:
    """Return ``(Ω1, Ω2)`` for ``Ω_eff = Ω1 Ω2 / (2 Δ_p)``.

    All frequencies use angular units.  ``omega1_over_omega2`` is ``Ω1/Ω2``.
    """

    omega = _check_omega(omega)
    delta_p = _check_delta_p(delta_p)
    omega1_over_omega2 = _check_beam_ratio(omega1_over_omega2)
    omega1 = np.sqrt(2.0 * omega * delta_p * omega1_over_omega2)
    omega2 = np.sqrt(2.0 * omega * delta_p / omega1_over_omega2)
    return float(omega1), float(omega2)


def scattering_rate_per_us(
    omega: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float = 1.0,
) -> float:
    """Return the effective dephasing collapse rate during a driven pulse.

    The eliminated 5P population receives contributions from both two-photon
    legs: ``Γ_e (Ω1² + Ω2²)/(4 Δ_p²)``.  With ``L=sqrt(rate)|r><r|``, the
    ground-Rydberg coherence decays at half this collapse rate, as expected for
    a projective scattering event.
    """

    omega = _check_omega(omega)
    gamma_e = _check_gamma(gamma_e)
    delta_p = _check_delta_p(delta_p)
    omega1, omega2 = single_photon_rabi_frequencies(omega, delta_p, omega1_over_omega2)
    return float(gamma_e * (omega1**2 + omega2**2) / (4.0 * delta_p**2))


def _dephasing_collapse_operators(rate: float, atom_index: int) -> list[Qobj]:
    """Return ``sqrt(rate) |r_i><r_i|`` for the driven atom only."""

    rate = _check_gamma(rate)
    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")
    if rate == 0.0:
        return []

    size = len(DECAY_BASIS)
    data = np.zeros((size, size), dtype=complex)
    labels = ("r0", "rg") if atom_index == 1 else ("0r", "gr")
    amplitude = np.sqrt(rate)
    for label in labels:
        data[_BASIS_INDEX[label], _BASIS_INDEX[label]] = amplitude
    return [_qobj(data)]


def _run_three_pulse_sequence(
    rho0: Qobj,
    omega: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float,
) -> Qobj:
    t_pi = np.pi / omega
    h1 = _hamiltonian(omega, atom_index=1)
    h2 = _hamiltonian(omega, atom_index=2)
    rate = scattering_rate_per_us(omega, gamma_e, delta_p, omega1_over_omega2)
    c_ops_atom1 = _dephasing_collapse_operators(rate, atom_index=1)
    c_ops_atom2 = _dephasing_collapse_operators(rate, atom_index=2)

    rho = _mesolve_final(h1, rho0, t_pi, c_ops_atom1)
    rho = _mesolve_final(h2, rho, 2.0 * t_pi, c_ops_atom2)
    rho = _mesolve_final(h1, rho, t_pi, c_ops_atom1)
    return rho


def _propagate_process_outputs(
    omega: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float,
) -> dict[tuple[int, int], np.ndarray]:
    outputs: dict[tuple[int, int], np.ndarray] = {}
    for row in range(_COMPUTATIONAL_DIM):
        for col in range(_COMPUTATIONAL_DIM):
            basis_operator = np.zeros((_COMPUTATIONAL_DIM, _COMPUTATIONAL_DIM), dtype=complex)
            basis_operator[row, col] = 1.0
            rho0 = _embed_computational_operator(basis_operator)
            outputs[(row, col)] = _project_to_computational(
                _run_three_pulse_sequence(rho0, omega, gamma_e, delta_p, omega1_over_omega2)
            )
    return outputs


def _zero_scattering_result(
    omega: float,
    gamma_e: float,
    delta_p: float,
    omega1_over_omega2: float,
) -> ScatteringGateResult:
    zero = run_decay_gate(omega=omega, gamma=0.0)
    omega1, omega2 = single_photon_rabi_frequencies(omega, delta_p, omega1_over_omega2)
    return ScatteringGateResult(
        omega=omega,
        gamma_e=gamma_e,
        delta_p=delta_p,
        omega1_over_omega2=omega1_over_omega2,
        omega1=omega1,
        omega2=omega2,
        scattering_rate_per_us=0.0,
        t_pi=zero.t_pi,
        average_gate_fidelity=zero.average_gate_fidelity,
        raw_average_gate_fidelity=zero.raw_average_gate_fidelity,
        kraus_ops=zero.kraus_ops,
        corrected_kraus_ops=zero.corrected_kraus_ops,
        choi_matrix=zero.choi_matrix,
        final_density_matrices=zero.final_density_matrices,
    )


def run_scattering_gate(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    gamma_e: float = 0.0,
    delta_p: float = 2.0 * np.pi * 1000.0,
    omega1_over_omega2: float = 1.0,
) -> ScatteringGateResult:
    """Run the isolated intermediate-state-scattering channel.

    ``omega`` is the effective two-photon angular Rabi frequency. ``gamma_e`` is
    the intermediate-state decay rate and ``delta_p`` is the one-photon detuning,
    both in rad/us.  Scattering collapse operators are applied only to the atom
    being actively driven in each pulse segment.
    """

    omega = _check_omega(omega)
    gamma_e = _check_gamma(gamma_e)
    delta_p = _check_delta_p(delta_p)
    omega1_over_omega2 = _check_beam_ratio(omega1_over_omega2)
    omega1, omega2 = single_photon_rabi_frequencies(omega, delta_p, omega1_over_omega2)
    rate = scattering_rate_per_us(omega, gamma_e, delta_p, omega1_over_omega2)

    if rate == 0.0:
        return _zero_scattering_result(omega, gamma_e, delta_p, omega1_over_omega2)

    outputs = _propagate_process_outputs(omega, gamma_e, delta_p, omega1_over_omega2)
    choi = _choi_from_process_outputs(outputs)
    kraus_ops = _kraus_from_choi(choi)
    corrected_kraus_ops = tuple(LOCAL_Z_PRODUCT @ op for op in kraus_ops)

    raw_fidelity = pedersen_fidelity_kraus(kraus_ops, CZ_TARGET)
    corrected_fidelity = pedersen_fidelity_kraus(corrected_kraus_ops, CZ_TARGET)
    final_density_matrices = {
        label: outputs[(index, index)]
        for index, label in enumerate(COMPUTATIONAL_BASIS)
    }

    return ScatteringGateResult(
        omega=omega,
        gamma_e=gamma_e,
        delta_p=delta_p,
        omega1_over_omega2=omega1_over_omega2,
        omega1=omega1,
        omega2=omega2,
        scattering_rate_per_us=rate,
        t_pi=np.pi / omega,
        average_gate_fidelity=corrected_fidelity,
        raw_average_gate_fidelity=raw_fidelity,
        kraus_ops=kraus_ops,
        corrected_kraus_ops=corrected_kraus_ops,
        choi_matrix=choi,
        final_density_matrices=final_density_matrices,
    )
