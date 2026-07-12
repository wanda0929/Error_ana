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

from src.params import get_rydberg_params
from src.sweeps import read_decay_sweep_csv, sweep_decay, write_decay_sweep_csv


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
DECAY_FIGURE = ROOT / "figures" / "fidelity_vs_decay_rate.png"


def _load_or_create_decay_rows():
    if DECAY_SWEEP_CSV.exists():
        return read_decay_sweep_csv(DECAY_SWEEP_CSV)
    rows = sweep_decay(num_points=25, decades=2.0)
    write_decay_sweep_csv(rows, DECAY_SWEEP_CSV)
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


def main() -> None:
    output_path = plot_decay_sweep()
    print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
