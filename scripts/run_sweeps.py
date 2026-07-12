#!/usr/bin/env python3
"""Run parameter sweeps for implemented Rydberg CZ error channels."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.analytical import epsilon_decay
from src.params import get_rydberg_params
from src.sweeps import sweep_decay, write_decay_sweep_csv


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"


def main() -> None:
    params = get_rydberg_params()
    rows = sweep_decay(num_points=25, decades=2.0)
    output_path = write_decay_sweep_csv(rows, DECAY_SWEEP_CSV)
    baseline_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)

    print(params.summary())
    print(f"Decay sweep points: {len(rows)}")
    print(
        "Gamma range: "
        f"{rows[0].gamma_per_us:.4e} to {rows[-1].gamma_per_us:.4e} us^-1 "
        f"({rows[-1].gamma_per_us / rows[0].gamma_per_us:.1f}x)"
    )
    print(f"Baseline gamma: {params.rydberg_decay_rate_per_us:.4e} us^-1")
    print(f"Baseline analytical decay error: {baseline_error:.4e}")
    print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
