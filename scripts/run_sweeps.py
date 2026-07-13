#!/usr/bin/env python3
"""Run parameter sweeps for implemented Rydberg CZ error channels."""

from __future__ import annotations

import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.analytical import epsilon_blockade, epsilon_decay, epsilon_scattering
from src.params import get_rydberg_params
from src.sweeps import (
    sweep_blockade,
    sweep_decay,
    sweep_scattering,
    write_blockade_sweep_csv,
    write_decay_sweep_csv,
    write_scattering_sweep_csv,
)


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
BLOCKADE_SWEEP_CSV = ROOT / "figures" / "blockade_sweep.csv"
SCATTERING_SWEEP_CSV = ROOT / "figures" / "scattering_sweep.csv"
BASELINE_INTERMEDIATE_DETUNING_MHZ = 1000.0


def main() -> None:
    params = get_rydberg_params()

    decay_rows = sweep_decay(num_points=25, decades=2.0)
    decay_path = write_decay_sweep_csv(decay_rows, DECAY_SWEEP_CSV)
    baseline_decay_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)

    blockade_rows = sweep_blockade(params.omega_rad_per_us, n_steps_per_pi=160)
    blockade_path = write_blockade_sweep_csv(blockade_rows, BLOCKADE_SWEEP_CSV)
    baseline_blockade_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)

    scattering_rows = sweep_scattering(num_points=25)
    scattering_path = write_scattering_sweep_csv(scattering_rows, SCATTERING_SWEEP_CSV)
    baseline_delta_p = 2.0 * math.pi * BASELINE_INTERMEDIATE_DETUNING_MHZ
    baseline_scattering_error = epsilon_scattering(params.intermediate_decay_rate_per_us, baseline_delta_p)

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

    print(f"Scattering sweep points: {len(scattering_rows)}")
    print(
        "Intermediate detuning range: "
        f"{scattering_rows[0].delta_p_mhz / 1000.0:.3g} to "
        f"{scattering_rows[-1].delta_p_mhz / 1000.0:.3g} GHz"
    )
    print(f"Baseline intermediate detuning: {BASELINE_INTERMEDIATE_DETUNING_MHZ / 1000.0:.1f} GHz")
    print(f"Baseline analytical scattering error: {baseline_scattering_error:.4e}")
    print(f"Saved {scattering_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
