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

from src.analytical import epsilon_blockade, epsilon_scattering
from src.params import get_rydberg_params
from src.sweeps import (
    blockade_ratios,
    read_blockade_sweep_csv,
    read_decay_sweep_csv,
    read_scattering_sweep_csv,
    scattering_detunings_mhz,
    sweep_blockade,
    sweep_decay,
    sweep_scattering,
    write_blockade_sweep_csv,
    write_decay_sweep_csv,
    write_scattering_sweep_csv,
)


DECAY_SWEEP_CSV = ROOT / "figures" / "decay_sweep.csv"
DECAY_FIGURE = ROOT / "figures" / "fidelity_vs_decay_rate.png"
BLOCKADE_SWEEP_CSV = ROOT / "figures" / "blockade_sweep.csv"
BLOCKADE_FIGURE = ROOT / "figures" / "fidelity_vs_blockade.png"
SCATTERING_SWEEP_CSV = ROOT / "figures" / "scattering_sweep.csv"
SCATTERING_FIGURE = ROOT / "figures" / "fidelity_vs_detuning.png"
BASELINE_INTERMEDIATE_DETUNING_MHZ = 1000.0


def _load_or_create_decay_rows():
    if DECAY_SWEEP_CSV.exists():
        return read_decay_sweep_csv(DECAY_SWEEP_CSV)
    rows = sweep_decay(num_points=25, decades=2.0)
    write_decay_sweep_csv(rows, DECAY_SWEEP_CSV)
    return rows


def _load_or_create_blockade_rows():
    if BLOCKADE_SWEEP_CSV.exists():
        return read_blockade_sweep_csv(BLOCKADE_SWEEP_CSV)
    rows = sweep_blockade(n_steps_per_pi=160)
    write_blockade_sweep_csv(rows, BLOCKADE_SWEEP_CSV)
    return rows


def _load_or_create_scattering_rows():
    if SCATTERING_SWEEP_CSV.exists():
        return read_scattering_sweep_csv(SCATTERING_SWEEP_CSV)
    rows = sweep_scattering(num_points=25)
    write_scattering_sweep_csv(rows, SCATTERING_SWEEP_CSV)
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


def plot_blockade_sweep(output_path: Path = BLOCKADE_FIGURE, *, n_steps_per_pi: int = 160) -> Path:
    """Plot numerical and analytical infidelity versus blockade ratio."""

    params = get_rydberg_params()
    rows = _load_or_create_blockade_rows()
    x = np.array([row.blockade_to_rabi for row in rows], dtype=float)
    numerical = np.array([row.numerical_error for row in rows], dtype=float)

    if x.size < 2:
        ratios = blockade_ratios(num=32, minimum=5.0, maximum=500.0)
        rows = sweep_blockade(params.omega_rad_per_us, ratios=ratios, n_steps_per_pi=n_steps_per_pi)
        write_blockade_sweep_csv(rows, BLOCKADE_SWEEP_CSV)
        x = np.array([row.blockade_to_rabi for row in rows], dtype=float)
        numerical = np.array([row.numerical_error for row in rows], dtype=float)

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
    return output_path


def plot_scattering_sweep(output_path: Path = SCATTERING_FIGURE) -> Path:
    """Plot numerical and analytical infidelity versus intermediate detuning."""

    params = get_rydberg_params()
    rows = _load_or_create_scattering_rows()
    x_ghz = np.array([row.delta_p_mhz / 1000.0 for row in rows], dtype=float)
    numerical = np.array([row.numerical_error for row in rows], dtype=float)
    analytical_points = np.array([row.analytical_error for row in rows], dtype=float)

    if x_ghz.size < 2:
        detunings = scattering_detunings_mhz(num=25)
        rows = sweep_scattering(detunings_mhz=detunings)
        write_scattering_sweep_csv(rows, SCATTERING_SWEEP_CSV)
        x_ghz = np.array([row.delta_p_mhz / 1000.0 for row in rows], dtype=float)
        numerical = np.array([row.numerical_error for row in rows], dtype=float)
        analytical_points = np.array([row.analytical_error for row in rows], dtype=float)

    x_curve_mhz = np.geomspace(x_ghz.min() * 1000.0, x_ghz.max() * 1000.0, 400)
    analytical_curve = np.array(
        [
            epsilon_scattering(params.intermediate_decay_rate_per_us, 2.0 * np.pi * detuning_mhz)
            for detuning_mhz in x_curve_mhz
        ],
        dtype=float,
    )
    baseline_error = epsilon_scattering(
        params.intermediate_decay_rate_per_us,
        2.0 * np.pi * BASELINE_INTERMEDIATE_DETUNING_MHZ,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    ax.plot(
        x_curve_mhz / 1000.0,
        analytical_curve,
        color="#5e81ac",
        lw=2.4,
        label=r"Analytical $\propto \Gamma_e/\Delta_p$",
    )
    ax.scatter(
        x_ghz,
        numerical,
        s=38,
        color="#b48ead",
        edgecolor="#2e3440",
        linewidth=0.7,
        zorder=3,
        label="Lindblad simulation",
    )
    ax.scatter(
        x_ghz,
        analytical_points,
        s=18,
        color="#5e81ac",
        alpha=0.45,
        label="Analytical sweep points",
    )
    ax.axvline(
        BASELINE_INTERMEDIATE_DETUNING_MHZ / 1000.0,
        color="#a3be8c",
        lw=1.6,
        ls="--",
        label=r"Evered-like $\Delta_p/2\pi=1\,GHz$",
    )
    ax.scatter(
        [BASELINE_INTERMEDIATE_DETUNING_MHZ / 1000.0],
        [baseline_error],
        marker="*",
        s=145,
        color="#ebcb8b",
        edgecolor="#2e3440",
        zorder=4,
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("Intermediate-state scattering in the Rydberg CZ gate")
    ax.set_xlabel(r"Intermediate detuning $\Delta_p/2\pi$ (GHz)")
    ax.set_ylabel(r"CZ infidelity $1-F_{\mathrm{avg}}$")
    ax.grid(True, which="both", alpha=0.22)
    ax.set_xlim(x_ghz.min() * 0.92, x_ghz.max() * 1.08)
    ax.legend(frameon=False, loc="upper right")
    ax.text(
        0.04,
        0.06,
        "larger detuning → less virtual 5P scattering",
        transform=ax.transAxes,
        color="#4c566a",
        fontsize=10,
    )
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def main() -> None:
    outputs = [plot_decay_sweep(), plot_blockade_sweep(), plot_scattering_sweep()]
    for output_path in outputs:
        print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
