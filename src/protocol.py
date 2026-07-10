"""Ideal pi-2pi-pi Rydberg blockade CZ protocol.

The qubit state |0> is dark.  The qubit state |1> is the ground state |g>
that couples to |r>.  For the ideal gate, |rr> is projected out: this is the
infinite-blockade model, not a numerically awkward large-U model.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from qutip import Qobj, basis, sesolve

from .fidelity import entangling_phase_from_diagonal
from .hamiltonian import build_ideal_hamiltonian, single_atom_hamiltonian
from .params import DEFAULT_OMEGA_RAD_PER_US

COMPUTATIONAL_BASIS = ("00", "01", "10", "11")
DEFAULT_SOLVER_OPTIONS = {
    "store_states": True,
    "nsteps": 10_000,
    "atol": 1e-12,
    "rtol": 1e-12,
}


def _check_omega(omega: float) -> None:
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError(f"omega must be positive and finite, got {omega!r}")


def _check_steps(n_steps_per_pi: int) -> None:
    if not isinstance(n_steps_per_pi, int) or n_steps_per_pi < 2:
        raise ValueError("n_steps_per_pi must be an integer >= 2")


def _tlist(duration: float, steps: int) -> np.ndarray:
    return np.linspace(0.0, duration, steps + 1)


def _solve_segment(hamiltonian: Qobj, psi0: Qobj, duration: float, steps: int) -> tuple[np.ndarray, list[Qobj]]:
    times = _tlist(duration, steps)
    result = sesolve(hamiltonian, psi0, times, e_ops=[], options=DEFAULT_SOLVER_OPTIONS)
    return times, list(result.states)


def _hold_segment(psi0: Qobj, duration: float, steps: int) -> tuple[np.ndarray, list[Qobj]]:
    times = _tlist(duration, steps)
    return times, [psi0] * len(times)


def _stitch_segments(
    segments: list[tuple[np.ndarray, list[Qobj]]],
    population: Callable[[Qobj], float],
) -> tuple[np.ndarray, np.ndarray, Qobj]:
    global_times: list[np.ndarray] = []
    global_pops: list[np.ndarray] = []
    offset = 0.0
    final_state = segments[0][1][0]

    for index, (local_times, states) in enumerate(segments):
        times = offset + local_times
        segment_states = states
        if index:
            times = times[1:]
            segment_states = segment_states[1:]
        global_times.append(times)
        global_pops.append(np.array([population(state) for state in segment_states], dtype=float))
        offset += float(local_times[-1])
        final_state = states[-1]

    return np.concatenate(global_times), np.concatenate(global_pops), final_state


def _amplitude(state: Qobj, index: int) -> complex:
    return complex(state.full()[index, 0])


def _single_atom_population(state: Qobj) -> float:
    return abs(_amplitude(state, 1)) ** 2


def _ideal_pair_population(state: Qobj) -> float:
    # Ideal basis |gg>, |gr>, |rg>.  Either of the last two has one Rydberg atom.
    return abs(_amplitude(state, 1)) ** 2 + abs(_amplitude(state, 2)) ** 2


def _run_dark_input(total_time: float, total_points: int) -> tuple[np.ndarray, np.ndarray, complex]:
    times = np.linspace(0.0, total_time, total_points)
    return times, np.zeros_like(times), 1.0 + 0.0j


def _run_input_01(omega: float, t_pi: float, steps_pi: int) -> tuple[np.ndarray, np.ndarray, complex]:
    psi_g = basis(2, 0)
    h_single = single_atom_hamiltonian(omega)
    segments = [
        _hold_segment(psi_g, t_pi, steps_pi),
        _solve_segment(h_single, psi_g, 2.0 * t_pi, 2 * steps_pi),
    ]
    _, _, after_2pi = _stitch_segments(segments, _single_atom_population)
    segments.append(_hold_segment(after_2pi, t_pi, steps_pi))
    times, pops, final_state = _stitch_segments(segments, _single_atom_population)
    return times, pops, _amplitude(final_state, 0)


def _run_input_10(omega: float, t_pi: float, steps_pi: int) -> tuple[np.ndarray, np.ndarray, complex]:
    psi_g = basis(2, 0)
    h_single = single_atom_hamiltonian(omega)
    first = _solve_segment(h_single, psi_g, t_pi, steps_pi)
    after_pi = first[1][-1]
    second = _hold_segment(after_pi, 2.0 * t_pi, 2 * steps_pi)
    after_hold = second[1][-1]
    third = _solve_segment(h_single, after_hold, t_pi, steps_pi)
    times, pops, final_state = _stitch_segments([first, second, third], _single_atom_population)
    return times, pops, _amplitude(final_state, 0)


def _run_input_11(omega: float, t_pi: float, steps_pi: int) -> tuple[np.ndarray, np.ndarray, complex]:
    psi_gg = basis(3, 0)
    h_atom1 = build_ideal_hamiltonian(omega, atom_index=1)
    h_atom2 = build_ideal_hamiltonian(omega, atom_index=2)
    first = _solve_segment(h_atom1, psi_gg, t_pi, steps_pi)
    second = _solve_segment(h_atom2, first[1][-1], 2.0 * t_pi, 2 * steps_pi)
    third = _solve_segment(h_atom1, second[1][-1], t_pi, steps_pi)
    times, pops, final_state = _stitch_segments([first, second, third], _ideal_pair_population)
    return times, pops, _amplitude(final_state, 0)


def run_ideal_gate(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    *,
    n_steps_per_pi: int = 400,
    n_steps: int | None = None,
) -> dict[str, object]:
    """Run the ideal blockade CZ sequence and return phases plus populations.

    ``omega`` is an angular Rabi frequency in rad/us.  The returned raw unitary is
    diagonal ``diag(1, -1, -1, -1)`` up to numerical error; comparing to canonical
    CZ requires the local Z correction in :mod:`src.fidelity`.
    """

    if n_steps is not None:
        n_steps_per_pi = n_steps
    _check_omega(omega)
    _check_steps(n_steps_per_pi)

    t_pi = np.pi / omega
    total_time = 4.0 * t_pi
    total_points = 4 * n_steps_per_pi + 1

    time_00, pop_00, amp_00 = _run_dark_input(total_time, total_points)
    time_01, pop_01, amp_01 = _run_input_01(omega, t_pi, n_steps_per_pi)
    time_10, pop_10, amp_10 = _run_input_10(omega, t_pi, n_steps_per_pi)
    time_11, pop_11, amp_11 = _run_input_11(omega, t_pi, n_steps_per_pi)

    if not (
        np.allclose(time_00, time_01)
        and np.allclose(time_00, time_10)
        and np.allclose(time_00, time_11)
    ):
        raise RuntimeError("internal error: population traces do not share a time grid")

    diagonal = np.array([amp_00, amp_01, amp_10, amp_11], dtype=complex)
    unitary = np.diag(diagonal)
    phases = np.angle(diagonal)

    return {
        "basis": COMPUTATIONAL_BASIS,
        "omega": omega,
        "t_pi": t_pi,
        "pulse_boundaries": np.array([t_pi, 3.0 * t_pi]),
        "times": time_00,
        "populations": {
            "00": pop_00,
            "01": pop_01,
            "10": pop_10,
            "11": pop_11,
        },
        "unitary": unitary,
        "unitary_diagonal": diagonal,
        "phases": phases,
        "entangling_phase": entangling_phase_from_diagonal(diagonal),
    }
