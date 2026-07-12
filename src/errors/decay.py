"""Rydberg-decay Lindblad simulation for the pi-2pi-pi CZ gate.

The model keeps the dark state |0>, the laser-coupled ground state |g>, and the
Rydberg state |r> for each atom, but projects out |rr> to isolate decay in the
same infinite-blockade limit used by the ideal gate.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import numpy as np
from qutip import Qobj, mesolve

from ..fidelity import CZ_TARGET, LOCAL_Z_PRODUCT, pedersen_fidelity_kraus
from ..params import DEFAULT_OMEGA_RAD_PER_US

DECAY_BASIS: Final[tuple[str, ...]] = ("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg")
"""Two-atom basis with the doubly excited Rydberg state projected out."""

COMPUTATIONAL_BASIS: Final[tuple[str, ...]] = ("00", "01", "10", "11")
COMPUTATIONAL_TO_DECAY: Final[dict[str, str]] = {
    "00": "00",
    "01": "0g",
    "10": "g0",
    "11": "gg",
}

_BASIS_INDEX: Final[dict[str, int]] = {label: index for index, label in enumerate(DECAY_BASIS)}
_COMPUTATIONAL_INDICES: Final[tuple[int, ...]] = tuple(_BASIS_INDEX[COMPUTATIONAL_TO_DECAY[label]] for label in COMPUTATIONAL_BASIS)

DEFAULT_MESOLVE_OPTIONS: Final[dict[str, float | int | bool]] = {
    "store_states": True,
    "nsteps": 10_000,
    "atol": 1e-11,
    "rtol": 1e-11,
}


@dataclass(frozen=True)
class DecayGateResult:
    """Numerical result for the isolated Rydberg-decay channel."""

    omega: float
    gamma: float
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


@dataclass(frozen=True)
class DecayPopulationTrace:
    """Diagnostic Rydberg population trace for one computational input."""

    label: str
    times: np.ndarray
    rydberg_population: np.ndarray
    final_density_matrix: np.ndarray


def _check_omega(omega: float) -> float:
    omega = float(omega)
    if not np.isfinite(omega) or omega <= 0.0:
        raise ValueError(f"omega must be positive and finite, got {omega!r}")
    return omega


def _check_gamma(gamma: float) -> float:
    gamma = float(gamma)
    if not np.isfinite(gamma) or gamma < 0.0:
        raise ValueError(f"gamma must be non-negative and finite, got {gamma!r}")
    return gamma


def _check_steps(n_steps_per_pi: int) -> int:
    if not isinstance(n_steps_per_pi, int) or n_steps_per_pi < 2:
        raise ValueError("n_steps_per_pi must be an integer >= 2")
    return n_steps_per_pi


def _qobj(data: np.ndarray) -> Qobj:
    return Qobj(np.asarray(data, dtype=complex))


def _hamiltonian(omega: float, atom_index: int) -> Qobj:
    """Return the infinite-blockade Hamiltonian in :data:`DECAY_BASIS`."""

    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")

    data = np.zeros((len(DECAY_BASIS), len(DECAY_BASIS)), dtype=complex)
    if atom_index == 1:
        couplings = (("g0", "r0"), ("gg", "rg"))
    else:
        couplings = (("0g", "0r"), ("gg", "gr"))

    for left, right in couplings:
        i = _BASIS_INDEX[left]
        j = _BASIS_INDEX[right]
        data[i, j] = data[j, i] = omega / 2.0
    return _qobj(data)


def _collapse_operators(gamma: float) -> list[Qobj]:
    """Return ``sqrt(gamma) |g><r|`` collapse operators for both atoms."""

    if gamma == 0.0:
        return []

    size = len(DECAY_BASIS)
    amplitude = np.sqrt(gamma)
    atom1 = np.zeros((size, size), dtype=complex)
    atom2 = np.zeros((size, size), dtype=complex)

    for destination, source in (("g0", "r0"), ("gg", "rg")):
        atom1[_BASIS_INDEX[destination], _BASIS_INDEX[source]] = amplitude
    for destination, source in (("0g", "0r"), ("gg", "gr")):
        atom2[_BASIS_INDEX[destination], _BASIS_INDEX[source]] = amplitude

    return [_qobj(atom1), _qobj(atom2)]


def _time_grid(duration: float, points: int) -> np.ndarray:
    return np.linspace(0.0, duration, points + 1)


def _mesolve_states(hamiltonian: Qobj, rho0: Qobj, duration: float, points: int, c_ops: list[Qobj]) -> list[Qobj]:
    result = mesolve(
        hamiltonian,
        rho0,
        _time_grid(duration, points),
        c_ops=c_ops,
        e_ops=[],
        options=DEFAULT_MESOLVE_OPTIONS,
    )
    return list(result.states)


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


def _embed_computational_operator(operator: np.ndarray) -> Qobj:
    matrix = np.asarray(operator, dtype=complex)
    if matrix.shape != (4, 4):
        raise ValueError(f"computational operator must be 4x4, got {matrix.shape}")

    embedded = np.zeros((len(DECAY_BASIS), len(DECAY_BASIS)), dtype=complex)
    for row, decay_row in enumerate(_COMPUTATIONAL_INDICES):
        for col, decay_col in enumerate(_COMPUTATIONAL_INDICES):
            embedded[decay_row, decay_col] = matrix[row, col]
    return _qobj(embedded)


def _project_to_computational(operator: Qobj) -> np.ndarray:
    full = np.asarray(operator.full(), dtype=complex)
    projected = full[np.ix_(_COMPUTATIONAL_INDICES, _COMPUTATIONAL_INDICES)]
    return np.asarray(projected, dtype=complex)


def _run_three_pulse_sequence(rho0: Qobj, omega: float, gamma: float) -> Qobj:
    t_pi = np.pi / omega
    h1 = _hamiltonian(omega, atom_index=1)
    h2 = _hamiltonian(omega, atom_index=2)
    c_ops = _collapse_operators(gamma)

    rho = _mesolve_final(h1, rho0, t_pi, c_ops)
    rho = _mesolve_final(h2, rho, 2.0 * t_pi, c_ops)
    rho = _mesolve_final(h1, rho, t_pi, c_ops)
    return rho


def _propagate_process_outputs(omega: float, gamma: float) -> dict[tuple[int, int], np.ndarray]:
    outputs: dict[tuple[int, int], np.ndarray] = {}
    for row in range(4):
        for col in range(4):
            basis_operator = np.zeros((4, 4), dtype=complex)
            basis_operator[row, col] = 1.0
            rho0 = _embed_computational_operator(basis_operator)
            outputs[(row, col)] = _project_to_computational(_run_three_pulse_sequence(rho0, omega, gamma))
    return outputs


def _choi_from_process_outputs(outputs: dict[tuple[int, int], np.ndarray]) -> np.ndarray:
    """Build the unnormalized Choi matrix from propagated operator outputs.

    The convention is input index first and output index second:
    ``J[(i,a),(j,b)] = E(|i><j|)[a,b]``.  The inverse in
    :func:`_kraus_from_choi` transposes the reshaped eigenvectors back to normal
    operator form.
    """

    dim = 4
    choi = np.zeros((dim * dim, dim * dim), dtype=complex)
    for input_row in range(dim):
        for input_col in range(dim):
            block = np.asarray(outputs[(input_row, input_col)], dtype=complex)
            if block.shape != (dim, dim):
                raise ValueError(f"process block must be 4x4, got {block.shape}")
            row_slice = slice(input_row * dim, (input_row + 1) * dim)
            col_slice = slice(input_col * dim, (input_col + 1) * dim)
            choi[row_slice, col_slice] = block

    return 0.5 * (choi + choi.conj().T)


def _kraus_from_choi(choi: np.ndarray, *, tolerance: float = 1e-12) -> tuple[np.ndarray, ...]:
    """Convert an unnormalized Choi matrix to Kraus operators."""

    matrix = np.asarray(choi, dtype=complex)
    if matrix.shape != (16, 16):
        raise ValueError(f"Choi matrix must be 16x16 for two qubits, got {matrix.shape}")

    eigenvalues, eigenvectors = np.linalg.eigh(0.5 * (matrix + matrix.conj().T))
    kraus_ops: list[np.ndarray] = []
    for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors.T):
        if eigenvalue <= tolerance:
            continue
        input_output = np.sqrt(float(eigenvalue)) * eigenvector.reshape((4, 4), order="C")
        kraus_ops.append(input_output.T.copy())

    if not kraus_ops:
        raise RuntimeError("Choi reconstruction produced no Kraus operators")
    return tuple(kraus_ops)


def run_decay_gate(
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    gamma: float = 0.0,
) -> DecayGateResult:
    """Run the isolated Rydberg-decay channel and compute Pedersen fidelity.

    ``omega`` is the angular Rabi frequency in rad/us and ``gamma`` is the
    Rydberg decay rate in us^-1.  The reported ``average_gate_fidelity`` applies
    the same local-Z frame correction used by the ideal gate before comparing to
    canonical CZ.
    """

    omega = _check_omega(omega)
    gamma = _check_gamma(gamma)
    if gamma == 0.0:
        return _zero_decay_result(omega)

    outputs = _propagate_process_outputs(omega, gamma)
    choi = _choi_from_process_outputs(outputs)
    kraus_ops = _kraus_from_choi(choi)
    corrected_kraus_ops = tuple(LOCAL_Z_PRODUCT @ op for op in kraus_ops)

    raw_fidelity = pedersen_fidelity_kraus(kraus_ops, CZ_TARGET)
    corrected_fidelity = pedersen_fidelity_kraus(corrected_kraus_ops, CZ_TARGET)

    final_density_matrices = {
        label: outputs[(index, index)]
        for index, label in enumerate(COMPUTATIONAL_BASIS)
    }

    return DecayGateResult(
        omega=omega,
        gamma=gamma,
        t_pi=np.pi / omega,
        average_gate_fidelity=corrected_fidelity,
        raw_average_gate_fidelity=raw_fidelity,
        kraus_ops=kraus_ops,
        corrected_kraus_ops=corrected_kraus_ops,
        choi_matrix=choi,
        final_density_matrices=final_density_matrices,
    )


def _zero_decay_result(omega: float) -> DecayGateResult:
    """Return the exact boundary result for the no-decay channel."""

    raw_gate = CZ_TARGET @ LOCAL_Z_PRODUCT
    corrected_gate = LOCAL_Z_PRODUCT @ raw_gate
    outputs: dict[tuple[int, int], np.ndarray] = {}
    for row in range(4):
        for col in range(4):
            basis_operator = np.zeros((4, 4), dtype=complex)
            basis_operator[row, col] = 1.0
            outputs[(row, col)] = raw_gate @ basis_operator @ raw_gate.conj().T
    choi = _choi_from_process_outputs(outputs)
    kraus_ops = (raw_gate,)
    corrected_kraus_ops = (corrected_gate,)
    final_density_matrices = {
        label: outputs[(index, index)]
        for index, label in enumerate(COMPUTATIONAL_BASIS)
    }
    return DecayGateResult(
        omega=omega,
        gamma=0.0,
        t_pi=np.pi / omega,
        average_gate_fidelity=pedersen_fidelity_kraus(corrected_kraus_ops, CZ_TARGET),
        raw_average_gate_fidelity=pedersen_fidelity_kraus(kraus_ops, CZ_TARGET),
        kraus_ops=kraus_ops,
        corrected_kraus_ops=corrected_kraus_ops,
        choi_matrix=choi,
        final_density_matrices=final_density_matrices,
    )


def run_decay_population_trace(
    label: str,
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    gamma: float = 0.0,
    *,
    n_steps_per_pi: int = 120,
) -> DecayPopulationTrace:
    """Return a diagnostic total-Rydberg-population trace for one basis input."""

    omega = _check_omega(omega)
    gamma = _check_gamma(gamma)
    n_steps_per_pi = _check_steps(n_steps_per_pi)
    if label not in COMPUTATIONAL_BASIS:
        raise ValueError(f"label must be one of {COMPUTATIONAL_BASIS}, got {label!r}")

    basis_index = COMPUTATIONAL_BASIS.index(label)
    rho_matrix = np.zeros((4, 4), dtype=complex)
    rho_matrix[basis_index, basis_index] = 1.0
    rho = _embed_computational_operator(rho_matrix)

    t_pi = np.pi / omega
    h1 = _hamiltonian(omega, atom_index=1)
    h2 = _hamiltonian(omega, atom_index=2)
    c_ops = _collapse_operators(gamma)

    segments = [
        (h1, t_pi, n_steps_per_pi),
        (h2, 2.0 * t_pi, 2 * n_steps_per_pi),
        (h1, t_pi, n_steps_per_pi),
    ]

    times_parts: list[np.ndarray] = []
    population_parts: list[np.ndarray] = []
    offset = 0.0
    for segment_index, (hamiltonian, duration, points) in enumerate(segments):
        states = _mesolve_states(hamiltonian, rho, duration, points, c_ops)
        local_times = _time_grid(duration, points)
        segment_times = offset + local_times
        populations = np.array([_total_rydberg_population(state) for state in states], dtype=float)
        if segment_index:
            segment_times = segment_times[1:]
            populations = populations[1:]
        times_parts.append(segment_times)
        population_parts.append(populations)
        offset += duration
        rho = states[-1]

    return DecayPopulationTrace(
        label=label,
        times=np.concatenate(times_parts),
        rydberg_population=np.concatenate(population_parts),
        final_density_matrix=_project_to_computational(rho),
    )


def _total_rydberg_population(rho: Qobj) -> float:
    diagonal = np.real(np.diag(rho.full()))
    total = 0.0
    for label in ("0r", "gr", "r0", "rg"):
        total += diagonal[_BASIS_INDEX[label]]
    return float(np.clip(total, 0.0, 1.0))
