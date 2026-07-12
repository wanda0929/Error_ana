"""Closed-form error estimates for the Rydberg CZ error budget."""

from __future__ import annotations

import numpy as np


def _check_positive_finite(value: float, name: str) -> None:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")


def epsilon_blockade(omega: float, blockade_shift: float) -> float:
    """Return the perturbative finite-blockade infidelity ``Ω²/(8U²)``.

    Parameters use the same angular-frequency units.  The formula only depends
    on the ratio, so ``omega`` and ``blockade_shift`` may be rad/us, rad/s, or
    any consistent unit.
    """

    _check_positive_finite(omega, "omega")
    _check_positive_finite(blockade_shift, "blockade_shift")
    return float(omega**2 / (8.0 * blockade_shift**2))
