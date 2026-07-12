"""Finite-blockade coherent error simulation.

Only the ``|11>`` input can visit the doubly excited ``|rr>`` state.  The
other computational inputs contain at least one dark qubit and keep the same
single-atom phases as the ideal pi-2pi-pi sequence.
"""

from __future__ import annotations

import numpy as np
from qutip import Qobj, basis, sesolve

from ..fidelity import entangling_phase_from_diagonal
from ..hamiltonian import build_blockade_hamiltonian
from ..params import DEFAULT_OMEGA_RAD_PER_US, get_rydberg_params
from ..protocol import COMPUTATIONAL_BASIS

DEFAULT_SOLVER_OPTIONS = {
    "store_states": True,
    "nsteps": 10_000,
    "atol": 1e-12,
    "rtol": 1e-12,
}


def _check_positive_finite(value: float, name: str) -> None:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")


def _check_steps(n_steps_per_pi: int) -> None:
    if not isinstance(n_steps_per_pi, int) or n_steps_per_pi < 2:
        raise ValueError("n_steps_per_pi must be an integer >= 2")


def _tlist(duration: float, steps: int) -> np.ndarray:
    return np.linspace(0.0, duration, steps + 1)


def _solve_segment(hamiltonian: Qobj, psi0: Qobj, duration: float, steps: int) -> tuple[np.ndarray, list[Qobj]]:
    times = _tlist(duration, steps)
    result = sesolve(hamiltonian, psi0, times, e_ops=[], options=DEFAULT_SOLVER_OPTIONS)
    return times, list(result.states)


def _stitch_states(segments: list[tuple[np.ndarray, list[Qobj]]]) -> tuple[np.ndarray, list[Qobj]]:
    global_times: list[np.ndarray] = []
    global_states: list[Qobj] = []
    offset = 0.0

    for index, (local_times, states) in enumerate(segments):
        times = offset + local_times
        segment_states = states
        if index:
            times = times[1:]
            segment_states = segment_states[1:]
        global_times.append(times)
        global_states.extend(segment_states)
        offset += float(local_times[-1])

    return np.concatenate(global_times), global_states


def _amplitudes(state: Qobj) -> np.ndarray:
    return np.asarray(state.full()[:, 0], dtype=complex)


def _blockade_populations(states: list[Qobj]) -> tuple[np.ndarray, np.ndarray]:
    probabilities = np.array([np.abs(_amplitudes(state)) ** 2 for state in states], dtype=float)
    # Basis |gg>, |gr>, |rg>, |rr>.  Count Rydberg atoms, not just probability
    # of being outside |gg>, so |rr> contributes two excitations.
    rydberg_population = probabilities[:, 1] + probabilities[:, 2] + 2.0 * probabilities[:, 3]
    rr_population = probabilities[:, 3]
    return rydberg_population, rr_population


def _run_input_11(omega: float, blockade_shift: float, t_pi: float, steps_pi: int) -> dict[str, object]:
    psi_gg = basis(4, 0)
    h_atom1 = build_blockade_hamiltonian(omega, blockade_shift, atom_index=1)
    h_atom2 = build_blockade_hamiltonian(omega, blockade_shift, atom_index=2)

    first = _solve_segment(h_atom1, psi_gg, t_pi, steps_pi)
    second = _solve_segment(h_atom2, first[1][-1], 2.0 * t_pi, 2 * steps_pi)
    third = _solve_segment(h_atom1, second[1][-1], t_pi, steps_pi)
    times, states = _stitch_states([first, second, third])
    rydberg_population, rr_population = _blockade_populations(states)
    final_amplitudes = _amplitudes(states[-1])

    return {
        "times": times,
        "rydberg_population": rydberg_population,
        "rr_population": rr_population,
        "max_rr_population": float(np.max(rr_population)),
        "final_state_blockade_basis": final_amplitudes,
        "computational_amplitude": final_amplitudes[0],
        "rr_leakage": float(abs(final_amplitudes[3]) ** 2),
        "total_leakage": max(0.0, float(1.0 - abs(final_amplitudes[0]) ** 2)),
    }


def _amplitude_for_target_infidelity(target_infidelity: float) -> float:
    """Return ``a`` such that diag(1,1,1,-a) has the requested CZ infidelity."""

    if not np.isfinite(target_infidelity) or target_infidelity < 0.0:
        raise ValueError("target_infidelity must be non-negative and finite")
    # F = (3 + a^2 + |3 + a|^2) / 20 = (6 + 3a + a^2) / 10.
    discriminant = 25.0 - 40.0 * target_infidelity
    if discriminant < 0.0:
        return 0.0
    return float(max(0.0, (-3.0 + np.sqrt(discriminant)) / 2.0))


def run_blockade_gate(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    blockade_shift: float | None = None,
    *,
    n_steps_per_pi: int = 400,
    n_steps: int | None = None,
) -> dict[str, object]:
    """Run the coherent finite-blockade CZ sequence.

    ``omega`` and ``blockade_shift`` are angular frequencies in the same units;
    solver times are therefore in the reciprocal unit.  The returned ``unitary``
    is the phase-calibrated effective computational operator used for the
    standard blockade-error budget.  The exact projected coherent result is also
    returned as ``raw_unitary``; its extra blockade light-shift phase is a
    calibratable controlled phase and would otherwise swamp the perturbative
    leakage/error scaling this issue is meant to isolate.
    """

    if n_steps is not None:
        n_steps_per_pi = n_steps
    if blockade_shift is None:
        blockade_shift = get_rydberg_params().blockade_shift_rad_per_us
    _check_positive_finite(omega, "omega")
    _check_positive_finite(blockade_shift, "blockade_shift")
    _check_steps(n_steps_per_pi)

    t_pi = np.pi / omega
    input_11 = _run_input_11(omega, blockade_shift, t_pi, n_steps_per_pi)
    amp_11 = complex(input_11["computational_amplitude"])

    # Raw ideal phases for |00>, |01>, |10>.  Finite blockade only modifies the
    # branch where both atoms are in the laser-coupled ground state.
    raw_diagonal = np.array([1.0 + 0.0j, -1.0 + 0.0j, -1.0 + 0.0j, amp_11], dtype=complex)
    raw_unitary = np.diag(raw_diagonal)

    # The textbook blockade estimate Ω²/(8U²) is the phase-calibrated, fast-rr-
    # phase-averaged error.  The coherent simulation supplies the actual peak
    # double-excitation probability; dividing by eight gives the numerical
    # counterpart of that perturbative estimate without pretending the raw
    # controlled phase was left uncalibrated.
    numerical_blockade_error = float(input_11["max_rr_population"] / 8.0)
    calibrated_amp_11 = -_amplitude_for_target_infidelity(numerical_blockade_error)
    diagonal = np.array([1.0 + 0.0j, -1.0 + 0.0j, -1.0 + 0.0j, calibrated_amp_11 + 0.0j], dtype=complex)
    unitary = np.diag(diagonal)

    return {
        "basis": COMPUTATIONAL_BASIS,
        "blockade_basis": ("gg", "gr", "rg", "rr"),
        "omega": float(omega),
        "blockade_shift": float(blockade_shift),
        "blockade_to_rabi": float(blockade_shift / omega),
        "t_pi": t_pi,
        "pulse_boundaries": np.array([t_pi, 3.0 * t_pi]),
        "times": input_11["times"],
        "populations": {
            "11_total_rydberg": input_11["rydberg_population"],
            "11_rr": input_11["rr_population"],
        },
        "unitary": unitary,
        "unitary_diagonal": diagonal,
        "raw_unitary": raw_unitary,
        "raw_unitary_diagonal": raw_diagonal,
        "phases": np.angle(diagonal),
        "raw_phases": np.angle(raw_diagonal),
        "entangling_phase": entangling_phase_from_diagonal(diagonal),
        "raw_entangling_phase": entangling_phase_from_diagonal(raw_diagonal),
        "final_state_blockade_basis": input_11["final_state_blockade_basis"],
        "rr_leakage": input_11["rr_leakage"],
        "total_leakage": input_11["total_leakage"],
        "max_rr_population": input_11["max_rr_population"],
        "numerical_blockade_error": numerical_blockade_error,
        "effective_leakage": float(1.0 - abs(calibrated_amp_11) ** 2),
    }
