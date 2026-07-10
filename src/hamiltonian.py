"""Hamiltonian builders for the Rb-Rb Rydberg CZ protocol.

Units are deliberately caller-owned: if ``omega`` is in rad/us, times must be
in us and blockade shifts must also be in rad/us.  The basis order is explicit
because silent basis swaps are how you get pretty plots of wrong physics.
"""

from __future__ import annotations

from typing import Final

import numpy as np
from qutip import Qobj

IDEAL_BASIS: Final[tuple[str, ...]] = ("gg", "gr", "rg")
"""Infinite-blockade basis.  The doubly excited state ``rr`` is projected out."""

BLOCKADE_BASIS: Final[tuple[str, ...]] = ("gg", "gr", "rg", "rr")
"""Finite-blockade basis used once ``rr`` leakage is intentionally modeled."""

SINGLE_ATOM_BASIS: Final[tuple[str, ...]] = ("g", "r")


def _check_atom_index(atom_index: int) -> None:
    if atom_index not in (1, 2):
        raise ValueError(f"atom_index must be 1 or 2, got {atom_index!r}")


def _check_rabi_frequency(omega: float) -> None:
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError(f"omega must be a positive finite Rabi frequency, got {omega!r}")


def single_atom_hamiltonian(omega: float, detuning: float = 0.0) -> Qobj:
    """Return the resonant two-level Hamiltonian in basis ``|g>, |r>``.

    ``detuning`` follows the project convention H_rr = -Delta.
    """

    _check_rabi_frequency(omega)
    data = np.array(
        [[0.0, omega / 2.0], [omega / 2.0, -float(detuning)]],
        dtype=complex,
    )
    return Qobj(data)


def build_ideal_hamiltonian(omega: float, atom_index: int, detuning: float = 0.0) -> Qobj:
    """Return the infinite-blockade Hamiltonian in ``IDEAL_BASIS``.

    The basis is ``|gg>, |gr>, |rg>``.  Couplings that would enter ``|rr>`` are
    absent by construction; this is the ideal U -> infinity model, not a large-U
    approximation.
    """

    _check_rabi_frequency(omega)
    _check_atom_index(atom_index)

    data = np.zeros((3, 3), dtype=complex)
    if atom_index == 1:
        # Pulse atom 1: |gg> <-> |rg|.  |gr> -> |rr> is projected out.
        data[0, 2] = data[2, 0] = omega / 2.0
        data[2, 2] = -float(detuning)
    else:
        # Pulse atom 2: |gg> <-> |gr|.  |rg> -> |rr> is projected out.
        data[0, 1] = data[1, 0] = omega / 2.0
        data[1, 1] = -float(detuning)
    return Qobj(data)


def build_blockade_hamiltonian(
    omega: float,
    blockade_shift: float,
    atom_index: int,
    detuning: float = 0.0,
) -> Qobj:
    """Return the finite-blockade Hamiltonian in ``BLOCKADE_BASIS``.

    Basis order: ``|gg>, |gr>, |rg>, |rr>``.  The ``rr`` state sits at energy
    ``blockade_shift``.  This builder is present now so later finite-blockade
    work has one obvious seam; the ideal protocol intentionally does not use it.
    """

    _check_rabi_frequency(omega)
    _check_atom_index(atom_index)
    if not np.isfinite(blockade_shift):
        raise ValueError("blockade_shift must be finite")

    data = np.zeros((4, 4), dtype=complex)
    if atom_index == 1:
        # Atom 1 flips in |gg> <-> |rg> and |gr> <-> |rr>.
        couplings = ((0, 2), (1, 3))
        detuned_states = (2, 3)
    else:
        # Atom 2 flips in |gg> <-> |gr> and |rg> <-> |rr>.
        couplings = ((0, 1), (2, 3))
        detuned_states = (1, 3)

    for i, j in couplings:
        data[i, j] = data[j, i] = omega / 2.0
    for i in detuned_states:
        data[i, i] += -float(detuning)
    data[3, 3] += float(blockade_shift)
    return Qobj(data)
