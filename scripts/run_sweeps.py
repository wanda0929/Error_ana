#!/usr/bin/env python3
"""Run finite-blockade parameter sweeps and save tabular data."""

from __future__ import annotations

import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.params import get_rydberg_params
from src.sweeps import sweep_blockade


def write_blockade_sweep(output_path: Path, *, n_steps_per_pi: int = 160) -> list[dict[str, float]]:
    params = get_rydberg_params()
    rows = sweep_blockade(params.omega_rad_per_us, n_steps_per_pi=n_steps_per_pi)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "blockade_to_rabi",
        "blockade_shift",
        "fidelity",
        "numerical_error",
        "analytical_error",
        "analytical_fidelity",
        "rr_leakage",
        "total_leakage",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main() -> None:
    output_path = ROOT / "data" / "blockade_sweep.csv"
    rows = write_blockade_sweep(output_path)
    baseline = get_rydberg_params().blockade_to_rabi
    print(f"Wrote {len(rows)} blockade sweep points to {output_path.relative_to(ROOT)}")
    print(f"Baseline U/Omega = {baseline:.2f}")


if __name__ == "__main__":
    main()
