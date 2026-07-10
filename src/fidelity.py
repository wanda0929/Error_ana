"""Pedersen average gate fidelity utilities.

The bare pi-2pi-pi blockade protocol produces CZ with known local Z phases.  We
keep the raw Pedersen formula separate from the optional local-phase correction
so tests can catch both the physics and the bookkeeping.
"""

from __future__ import annotations

import numpy as np

CZ_TARGET = np.diag([1.0, 1.0, 1.0, -1.0]).astype(complex)
LOCAL_Z_PRODUCT = np.diag([1.0, -1.0, -1.0, 1.0]).astype(complex)


def as_operator_matrix(operator: object) -> np.ndarray:
    """Return ``operator`` as a square complex matrix.

    A one-dimensional input is treated as a diagonal operator.  This is handy for
    protocol results, which naturally report only the computational-basis phases.
    """

    if hasattr(operator, "full"):
        operator = operator.full()
    matrix = np.asarray(operator, dtype=complex)
    if matrix.ndim == 1:
        matrix = np.diag(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError(f"operator must be a square matrix or diagonal vector, got {matrix.shape}")
    return matrix


def _unit_interval_with_roundoff(value: float, tol: float = 1e-12) -> float:
    if -tol <= value < 0.0:
        return 0.0
    if 1.0 < value <= 1.0 + tol:
        return 1.0
    return value


def pedersen_fidelity(
    actual: object,
    target: object = CZ_TARGET,
    *,
    correct_local_z: bool = False,
) -> float:
    """Average gate fidelity between two operators using Pedersen's formula.

    ``correct_local_z=True`` applies the local Z phase correction appropriate for
    diagonal two-qubit gates before comparing against ``target``.
    """

    actual_matrix = as_operator_matrix(actual)
    target_matrix = as_operator_matrix(target)
    if actual_matrix.shape != target_matrix.shape:
        raise ValueError(f"shape mismatch: actual {actual_matrix.shape}, target {target_matrix.shape}")
    if correct_local_z:
        actual_matrix = apply_local_z_correction(actual_matrix)

    d = target_matrix.shape[0]
    m = target_matrix.conj().T @ actual_matrix
    fidelity = (np.trace(m.conj().T @ m).real + abs(np.trace(m)) ** 2) / (d * (d + 1))
    return _unit_interval_with_roundoff(float(np.real_if_close(fidelity)))


def pedersen_fidelity_kraus(kraus_ops: list[object] | tuple[object, ...], target: object = CZ_TARGET) -> float:
    """Average gate fidelity for an open-system map represented by Kraus operators."""

    target_matrix = as_operator_matrix(target)
    d = target_matrix.shape[0]
    term_norm = 0.0
    term_trace = 0.0
    for op in kraus_ops:
        k = as_operator_matrix(op)
        if k.shape != target_matrix.shape:
            raise ValueError(f"Kraus shape mismatch: {k.shape} vs target {target_matrix.shape}")
        m = target_matrix.conj().T @ k
        term_norm += np.trace(m.conj().T @ m).real
        term_trace += abs(np.trace(m)) ** 2
    fidelity = (term_norm + term_trace) / (d * (d + 1))
    return _unit_interval_with_roundoff(float(np.real_if_close(fidelity)))


def entangling_phase_from_diagonal(diagonal: object) -> float:
    """Return phi_11 - phi_10 - phi_01 + phi_00 wrapped to (-pi, pi]."""

    diag = np.diag(as_operator_matrix(diagonal))
    if diag.shape != (4,):
        raise ValueError("entangling phase is defined here for a two-qubit diagonal gate")
    phases = np.angle(diag)
    return wrap_phase(phases[3] - phases[2] - phases[1] + phases[0])


def wrap_phase(angle: float) -> float:
    """Wrap a phase to (-pi, pi]."""

    wrapped = (float(angle) + np.pi) % (2.0 * np.pi) - np.pi
    if wrapped <= -np.pi + 1e-15:
        return float(np.pi)
    return float(wrapped)


def optimal_local_z_correction(actual: object, *, min_diagonal_magnitude: float = 1e-9) -> np.ndarray:
    """Return the diagonal local phase correction for a two-qubit diagonal gate.

    For diagonal phases ``phi_00, phi_01, phi_10, phi_11``, local Z rotations can
    remove the first three phases.  The only invariant left is the entangling
    phase ``phi_11 - phi_10 - phi_01 + phi_00``.
    """

    matrix = as_operator_matrix(actual)
    if matrix.shape != (4, 4):
        raise ValueError("local Z correction is implemented for 4x4 two-qubit operators")
    diag = np.diag(matrix)
    if np.any(np.abs(diag) < min_diagonal_magnitude):
        raise ValueError("cannot infer local phases from near-zero diagonal entries")

    phi00, phi01, phi10, _phi11 = np.angle(diag)
    correction_diagonal = np.array(
        [
            np.exp(-1j * phi00),
            np.exp(-1j * phi01),
            np.exp(-1j * phi10),
            np.exp(1j * (phi00 - phi01 - phi10)),
        ],
        dtype=complex,
    )
    return np.diag(correction_diagonal)


def apply_local_z_correction(actual: object) -> np.ndarray:
    """Apply the inferred local Z phase correction to ``actual``."""

    matrix = as_operator_matrix(actual)
    return optimal_local_z_correction(matrix) @ matrix
