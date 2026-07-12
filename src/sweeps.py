"""Parameter sweeps for Rydberg CZ error channels."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from .analytical import epsilon_blockade
from .errors.blockade import run_blockade_gate
from .fidelity import CZ_TARGET, pedersen_fidelity
from .params import get_rydberg_params


def blockade_ratios(num: int = 30, minimum: float = 5.0, maximum: float = 500.0) -> np.ndarray:
    """Return log-spaced blockade ratios ``U/Ω`` for the standard sweep."""

    if not isinstance(num, int) or num < 2:
        raise ValueError("num must be an integer >= 2")
    if not np.isfinite(minimum) or not np.isfinite(maximum) or minimum <= 0 or maximum <= minimum:
        raise ValueError("minimum and maximum must be positive finite values with maximum > minimum")
    return np.geomspace(float(minimum), float(maximum), int(num))


def sweep_blockade(
    omega: float | None = None,
    blockade_shifts: Iterable[float] | None = None,
    *,
    ratios: Iterable[float] | None = None,
    n_steps_per_pi: int = 160,
) -> list[dict[str, float]]:
    """Sweep finite-blockade strength and return numerical/analytical errors.

    Provide either absolute ``blockade_shifts`` or dimensionless ``ratios``.  If
    neither is provided, the standard 30-point ``U/Ω`` sweep from 5 to 500 is
    used.
    """

    if omega is None:
        omega = get_rydberg_params().omega_rad_per_us
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError(f"omega must be positive and finite, got {omega!r}")

    if blockade_shifts is not None and ratios is not None:
        raise ValueError("provide either blockade_shifts or ratios, not both")
    if ratios is None:
        if blockade_shifts is None:
            ratio_values = blockade_ratios()
            shift_values = omega * ratio_values
        else:
            shift_values = np.asarray(list(blockade_shifts), dtype=float)
            ratio_values = shift_values / omega
    else:
        ratio_values = np.asarray(list(ratios), dtype=float)
        shift_values = omega * ratio_values

    if shift_values.ndim != 1 or len(shift_values) == 0:
        raise ValueError("sweep requires at least one blockade value")
    if np.any(~np.isfinite(shift_values)) or np.any(shift_values <= 0):
        raise ValueError("all blockade shifts must be positive and finite")
    if np.any(~np.isfinite(ratio_values)) or np.any(ratio_values <= 0):
        raise ValueError("all blockade ratios must be positive and finite")

    rows: list[dict[str, float]] = []
    for ratio, blockade_shift in zip(ratio_values, shift_values):
        result = run_blockade_gate(omega, blockade_shift, n_steps_per_pi=n_steps_per_pi)
        fidelity = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
        analytical_error = epsilon_blockade(omega, blockade_shift)
        rows.append(
            {
                "blockade_to_rabi": float(ratio),
                "blockade_shift": float(blockade_shift),
                "fidelity": float(fidelity),
                "numerical_error": float(1.0 - fidelity),
                "analytical_error": float(analytical_error),
                "analytical_fidelity": float(1.0 - analytical_error),
                "rr_leakage": float(result["rr_leakage"]),
                "total_leakage": float(result["total_leakage"]),
            }
        )
    return rows
