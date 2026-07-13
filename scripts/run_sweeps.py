#!/usr/bin/env python3
"""Run parameter sweeps for implemented Rydberg CZ error channels."""

from __future__ import annotations

import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.analytical import epsilon_amplitude, epsilon_blockade, epsilon_decay, epsilon_doppler, epsilon_scattering
from src.errors.amplitude import DEFAULT_SIGMA_OMEGA
from src.params import DEFAULT_TEMPERATURE_K, K_EFF_RAD_PER_UM, RB87_MASS_KG, get_rydberg_params
from src.sweeps import (
    combined_error_budget_rows,
    sweep_amplitude,
    sweep_blockade,
    sweep_decay,
    sweep_doppler,
    sweep_scattering,
    write_amplitude_sweep_csv,
    write_blockade_sweep_csv,
    write_combined_budget_csv,
    write_decay_sweep_csv,
    write_doppler_sweep_csv,
    write_scattering_sweep_csv,
)


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
BLOCKADE_SWEEP_CSV = ROOT / "figures" / "blockade_sweep.csv"
DOPPLER_SWEEP_CSV = ROOT / "figures" / "doppler_sweep.csv"
SCATTERING_SWEEP_CSV = ROOT / "figures" / "scattering_sweep.csv"
AMPLITUDE_SWEEP_CSV = ROOT / "figures" / "amplitude_sweep.csv"
COMBINED_BUDGET_CSV = ROOT / "figures" / "combined_error_budget.csv"
BASELINE_INTERMEDIATE_DETUNING_MHZ = 1000.0


def main() -> None:
    params = get_rydberg_params()

    decay_rows = sweep_decay(num_points=25, decades=2.0)
    decay_path = write_decay_sweep_csv(decay_rows, DECAY_SWEEP_CSV)
    baseline_decay_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)

    blockade_rows = sweep_blockade(params.omega_rad_per_us, n_steps_per_pi=160)
    blockade_path = write_blockade_sweep_csv(blockade_rows, BLOCKADE_SWEEP_CSV)
    baseline_blockade_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)

    doppler_rows = sweep_doppler(num_points=25, n_samples=500)
    doppler_path = write_doppler_sweep_csv(doppler_rows, DOPPLER_SWEEP_CSV)
    baseline_doppler_error = epsilon_doppler(
        K_EFF_RAD_PER_UM,
        DEFAULT_TEMPERATURE_K,
        RB87_MASS_KG,
        params.omega_rad_per_us,
    )

    scattering_rows = sweep_scattering(num_points=25)
    scattering_path = write_scattering_sweep_csv(scattering_rows, SCATTERING_SWEEP_CSV)
    baseline_delta_p = 2.0 * math.pi * BASELINE_INTERMEDIATE_DETUNING_MHZ
    baseline_scattering_error = epsilon_scattering(params.intermediate_decay_rate_per_us, baseline_delta_p)

    amplitude_rows = sweep_amplitude(num_points=25, n_samples=500)
    amplitude_path = write_amplitude_sweep_csv(amplitude_rows, AMPLITUDE_SWEEP_CSV)
    baseline_amplitude_error = epsilon_amplitude(DEFAULT_SIGMA_OMEGA)

    combined_rows = combined_error_budget_rows(n_samples=24, seed=2024, individual_n_samples=300)
    combined_path = write_combined_budget_csv(combined_rows, COMBINED_BUDGET_CSV)
    additive_row = next(row for row in combined_rows if row.source == "Total (additive)")
    combined_row = next(row for row in combined_rows if row.source == "Total (combined)")
    relative_gap = abs(additive_row.numerical_error - combined_row.numerical_error) / combined_row.numerical_error

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

    print(f"Doppler sweep points: {len(doppler_rows)}")
    print(
        "Temperature range: "
        f"{doppler_rows[0].temperature_uK:.3g} to {doppler_rows[-1].temperature_uK:.3g} uK "
        f"with N={doppler_rows[0].n_samples} shots/point"
    )
    print(f"Baseline temperature: {DEFAULT_TEMPERATURE_K * 1e6:.1f} uK")
    print(f"Baseline analytical Doppler error: {baseline_doppler_error:.4e}")
    print(f"Saved {doppler_path.relative_to(ROOT)}")

    print(f"Scattering sweep points: {len(scattering_rows)}")
    print(
        "Intermediate detuning range: "
        f"{scattering_rows[0].delta_p_mhz / 1000.0:.3g} to "
        f"{scattering_rows[-1].delta_p_mhz / 1000.0:.3g} GHz"
    )
    print(f"Baseline intermediate detuning: {BASELINE_INTERMEDIATE_DETUNING_MHZ / 1000.0:.1f} GHz")
    print(f"Baseline analytical scattering error: {baseline_scattering_error:.4e}")
    print(f"Saved {scattering_path.relative_to(ROOT)}")

    print(f"Amplitude-noise sweep points: {len(amplitude_rows)}")
    print(
        "Fractional Rabi-noise range: "
        f"{amplitude_rows[0].sigma_percent:.3g}% to {amplitude_rows[-1].sigma_percent:.3g}% "
        f"with N={amplitude_rows[0].n_samples} shots/point"
    )
    print(f"Evered-like fractional Rabi noise: {DEFAULT_SIGMA_OMEGA * 100.0:.1f}%")
    print(f"Baseline analytical amplitude-noise error: {baseline_amplitude_error:.4e}")
    print(f"Saved {amplitude_path.relative_to(ROOT)}")

    print("Combined baseline budget:")
    print(f"Additive numerical sum: {additive_row.numerical_error:.4e}")
    print(f"Full combined numerical error: {combined_row.numerical_error:.4e}")
    print(f"Relative additive-vs-combined gap: {relative_gap:.1%}")
    print(f"Saved {combined_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
