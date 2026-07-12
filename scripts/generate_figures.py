#!/usr/bin/env python3
"""Generate publication figures for implemented error channels."""

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

from src.analytical import epsilon_blockade
from src.params import get_rydberg_params
from src.sweeps import blockade_ratios, sweep_blockade


def plot_blockade_sweep(output_path: Path, *, n_steps_per_pi: int = 160) -> None:
    params = get_rydberg_params()
    ratios = blockade_ratios(num=32, minimum=5.0, maximum=500.0)
    rows = sweep_blockade(params.omega_rad_per_us, ratios=ratios, n_steps_per_pi=n_steps_per_pi)

    x = np.array([row["blockade_to_rabi"] for row in rows], dtype=float)
    numerical = np.array([row["numerical_error"] for row in rows], dtype=float)
    x_curve = np.geomspace(x.min(), x.max(), 400)
    analytical = np.array(
        [epsilon_blockade(params.omega_rad_per_us, params.omega_rad_per_us * ratio) for ratio in x_curve],
        dtype=float,
    )
    baseline_ratio = params.blockade_to_rabi
    baseline_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    ax.plot(x_curve, analytical, color="#5e81ac", lw=2.4, label=r"Analytical $\Omega^2/(8U^2)$")
    ax.scatter(
        x,
        numerical,
        s=34,
        color="#d08770",
        edgecolor="#2e3440",
        linewidth=0.7,
        zorder=3,
        label="QuTiP simulation",
    )
    ax.axvline(
        baseline_ratio,
        color="#a3be8c",
        lw=1.6,
        ls="--",
        label=fr"Baseline $U/\Omega={baseline_ratio:.0f}$",
    )
    ax.scatter([baseline_ratio], [baseline_error], marker="*", s=145, color="#ebcb8b", edgecolor="#2e3440", zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Finite blockade error in the Rydberg CZ gate")
    ax.set_xlabel(r"Blockade ratio $U/\Omega$")
    ax.set_ylabel(r"CZ infidelity $1-F_{\mathrm{avg}}$")
    ax.grid(True, which="both", alpha=0.22)
    ax.set_xlim(x.min() * 0.92, x.max() * 1.08)
    ax.legend(frameon=False, loc="upper right")
    ax.text(
        0.04,
        0.06,
        "larger blockade → smaller coherent error",
        transform=ax.transAxes,
        color="#4c566a",
        fontsize=10,
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def main() -> None:
    output_path = ROOT / "figures" / "fidelity_vs_blockade.png"
    plot_blockade_sweep(output_path)
    print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
