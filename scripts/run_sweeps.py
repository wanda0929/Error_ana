#!/usr/bin/env python3
"""Run parameter sweeps for implemented Rydberg CZ error channels."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.analytical import epsilon_blockade, epsilon_decay
from src.params import get_rydberg_params
from src.sweeps import sweep_blockade, sweep_decay, write_blockade_sweep_csv, write_decay_sweep_csv


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
BLOCKADE_SWEEP_CSV = ROOT / "figures" / "blockade_sweep.csv"


def main() -> None:
    params = get_rydberg_params()

    decay_rows = sweep_decay(num_points=25, decades=2.0)
    decay_path = write_decay_sweep_csv(decay_rows, DECAY_SWEEP_CSV)
    baseline_decay_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)

    blockade_rows = sweep_blockade(params.omega_rad_per_us, n_steps_per_pi=160)
    blockade_path = write_blockade_sweep_csv(blockade_rows, BLOCKADE_SWEEP_CSV)
    baseline_blockade_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)

    print(params.summary())
    print(f"Decay sweep points: {len(decay_rows)}")
    print(
        "Gamma range: "
        f"{decay_rows[0].gamma_per_us:.4e} to {decay_rows[-1].gamma_per_us:.4e} us^-1 "
        f"({decay_rows[-1].gamma_per_us / decay_rows[0].gamma_per_us:.1f}x)"
    )
    print(f"Baseline gamma: {params.rydberg_decay_rate_per_us:.4e} us^-1")
    print(f"Baseline analytical decay error: {baseline_decay_error:.4e}")
    print(f"Saved {decay_path.relative_to(ROOT)}")

    print(f"Blockade sweep points: {len(blockade_rows)}")
    print(
        "U/Omega range: "
        f"{blockade_rows[0].blockade_to_rabi:.2f} to {blockade_rows[-1].blockade_to_rabi:.2f}"
    )
    print(f"Baseline U/Omega: {params.blockade_to_rabi:.2f}")
    print(f"Baseline analytical blockade error: {baseline_blockade_error:.4e}")
    print(f"Saved {blockade_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
