#!/usr/bin/env python3
"""Generate publication figures for implemented Rydberg CZ error channels."""

from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.params import DEFAULT_TEMPERATURE_K, get_rydberg_params
from src.sweeps import (
    read_decay_sweep_csv,
    read_doppler_sweep_csv,
    sweep_decay,
    sweep_doppler,
    write_decay_sweep_csv,
    write_doppler_sweep_csv,
)


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
DOPPLER_SWEEP_CSV = ROOT / "figures" / "doppler_sweep.csv"
DECAY_FIGURE = ROOT / "figures" / "fidelity_vs_decay_rate.png"
DOPPLER_FIGURE = ROOT / "figures" / "fidelity_vs_temperature.png"


def _load_or_create_decay_rows():
    if DECAY_SWEEP_CSV.exists():
        return read_decay_sweep_csv(DECAY_SWEEP_CSV)
    rows = sweep_decay(num_points=25, decades=2.0)
    write_decay_sweep_csv(rows, DECAY_SWEEP_CSV)
    return rows


def _load_or_create_doppler_rows():
    if DOPPLER_SWEEP_CSV.exists():
        return read_doppler_sweep_csv(DOPPLER_SWEEP_CSV)
    rows = sweep_doppler(num_points=25, n_samples=500)
    write_doppler_sweep_csv(rows, DOPPLER_SWEEP_CSV)
    return rows


def plot_decay_sweep(output_path: Path = DECAY_FIGURE) -> Path:
    """Plot numerical and analytical fidelity versus decay rate."""

    params = get_rydberg_params()
    rows = _load_or_create_decay_rows()
    gamma = np.array([row.gamma_per_us for row in rows], dtype=float)
    numerical_fidelity = np.array([row.numerical_fidelity for row in rows], dtype=float)
    analytical_fidelity = np.array([row.analytical_fidelity for row in rows], dtype=float)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)

    ax.plot(gamma, analytical_fidelity, color="#5e81ac", lw=2.2, label=r"Analytical $1-\Gamma T_R$")
    ax.scatter(
        gamma,
        numerical_fidelity,
        s=42,
        color="#bf616a",
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
        label="Lindblad simulation",
    )
    ax.axvline(
        params.rydberg_decay_rate_per_us,
        color="#2e3440",
        lw=1.2,
        ls="--",
        alpha=0.72,
        label="Project ARC baseline",
    )

    baseline_text = (
        rf"$\tau={params.rydberg_lifetime_us:.0f}\,\mu s$" "\n"
        rf"$\Omega/2\pi={params.omega_mhz:.1f}\,MHz$"
    )
    ax.text(
        params.rydberg_decay_rate_per_us * 1.08,
        min(0.9996, max(analytical_fidelity) - 0.00025),
        baseline_text,
        ha="left",
        va="top",
        fontsize=9,
        color="#2e3440",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#d8dee9", "alpha": 0.92},
    )

    ax.set_xscale("log")
    ax.set_title("Rydberg decay error in the blockade CZ gate")
    ax.set_xlabel(r"Rydberg decay rate $\Gamma$ ($\mu s^{-1}$)")
    ax.set_ylabel(r"Average gate fidelity $F_{avg}$")
    ax.set_ylim(min(analytical_fidelity.min(), numerical_fidelity.min()) - 0.0007, 1.00015)
    ax.grid(True, which="major", alpha=0.25)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(frameon=False, loc="lower left")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def plot_doppler_sweep(output_path: Path = DOPPLER_FIGURE) -> Path:
    """Plot numerical and analytical fidelity versus atom temperature."""

    rows = _load_or_create_doppler_rows()
    temperature_uK = np.array([row.temperature_uK for row in rows], dtype=float)
    numerical_fidelity = np.array([row.numerical_fidelity for row in rows], dtype=float)
    numerical_se = np.array([row.numerical_standard_error for row in rows], dtype=float)
    analytical_fidelity = np.array([row.analytical_fidelity for row in rows], dtype=float)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)

    ax.plot(
        temperature_uK,
        analytical_fidelity,
        color="#5e81ac",
        lw=2.2,
        label=r"Analytical $1-\epsilon_D$",
    )
    ax.errorbar(
        temperature_uK,
        numerical_fidelity,
        yerr=numerical_se,
        fmt="o",
        ms=5.4,
        color="#bf616a",
        ecolor="#d08770",
        elinewidth=1.0,
        capsize=2.5,
        markeredgecolor="white",
        markeredgewidth=0.7,
        zorder=3,
        label="Monte Carlo coherent simulation",
    )
    ax.axvline(
        DEFAULT_TEMPERATURE_K * 1e6,
        color="#2e3440",
        lw=1.2,
        ls="--",
        alpha=0.72,
        label="Evered-like 10 µK point",
    )

    baseline_index = int(np.argmin(np.abs(temperature_uK - DEFAULT_TEMPERATURE_K * 1e6)))
    baseline_text = (
        rf"$T={DEFAULT_TEMPERATURE_K * 1e6:.0f}\,\mu K$" "\n"
        rf"$1-F\approx{1.0 - numerical_fidelity[baseline_index]:.1e}$"
    )
    ax.text(
        DEFAULT_TEMPERATURE_K * 1e6 * 1.08,
        min(0.99998, max(analytical_fidelity) - 0.00002),
        baseline_text,
        ha="left",
        va="top",
        fontsize=9,
        color="#2e3440",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#d8dee9", "alpha": 0.92},
    )

    ax.set_xscale("log")
    ax.set_title("Doppler dephasing in the blockade CZ gate")
    ax.set_xlabel(r"Atom temperature $T$ ($\mu K$)")
    ax.set_ylabel(r"Average gate fidelity $F_{avg}$")
    ax.set_ylim(min(analytical_fidelity.min(), numerical_fidelity.min()) - 0.00025, 1.00004)
    ax.grid(True, which="major", alpha=0.25)
    ax.grid(True, which="minor", alpha=0.10)
    ax.legend(frameon=False, loc="lower left")
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def main() -> None:
    decay_output_path = plot_decay_sweep()
    doppler_output_path = plot_doppler_sweep()
    print(f"Saved {decay_output_path.relative_to(ROOT)}")
    print(f"Saved {doppler_output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
