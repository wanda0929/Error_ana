#!/usr/bin/env python3
"""Run the ideal Rydberg CZ simulation and save the population-dynamics figure."""

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

from src.fidelity import CZ_TARGET, pedersen_fidelity
from src.params import get_rydberg_params
from src.protocol import run_ideal_gate


def plot_population_dynamics(result: dict[str, object], output_path: Path) -> None:
    times = np.asarray(result["times"])
    populations = result["populations"]
    boundaries = np.asarray(result["pulse_boundaries"])

    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
    colors = {
        "00": "#4c566a",
        "01": "#5e81ac",
        "10": "#d08770",
        "11": "#a3be8c",
    }
    for label in ("00", "01", "10", "11"):
        ax.plot(times, populations[label], lw=2.2, label=fr"$|{label}\rangle$", color=colors[label])

    for boundary in boundaries:
        ax.axvline(boundary, color="#2e3440", lw=1.0, ls="--", alpha=0.55)

    ax.text(boundaries[0] / 2, 1.04, r"$\pi_1$", ha="center", va="bottom")
    ax.text((boundaries[0] + boundaries[1]) / 2, 1.04, r"$2\pi_2$", ha="center", va="bottom")
    ax.text((boundaries[1] + times[-1]) / 2, 1.04, r"$\pi_1$", ha="center", va="bottom")

    ax.set_title("Ideal Rydberg blockade CZ population dynamics")
    ax.set_xlabel("Time (µs)")
    ax.set_ylabel("Total Rydberg population")
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(-0.03, 1.12)
    ax.grid(True, alpha=0.22)
    ax.legend(frameon=False, ncols=4, loc="upper center", bbox_to_anchor=(0.5, -0.16))
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def main() -> None:
    params = get_rydberg_params()
    result = run_ideal_gate(params.omega_rad_per_us)
    raw_fidelity = pedersen_fidelity(result["unitary"], CZ_TARGET)
    corrected_fidelity = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
    output_path = ROOT / "figures" / "population_dynamics.png"
    plot_population_dynamics(result, output_path)

    print(params.summary())
    print(f"Raw diagonal unitary: {np.diag(result['unitary'])}")
    print(f"Raw fidelity vs CZ: {raw_fidelity:.12f}")
    print(f"Local-Z-corrected fidelity vs CZ: {corrected_fidelity:.12f}")
    print(f"Entangling phase: {result['entangling_phase']:.12f} rad")
    print(f"Saved {output_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
