"""Parameter sweeps for Rydberg CZ error channels."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .analytical import epsilon_decay_from_gamma, epsilon_doppler
from .errors.decay import run_decay_gate
from .errors.doppler import run_doppler_gate_mc
from .params import (
    DEFAULT_OMEGA_RAD_PER_US,
    K_EFF_RAD_PER_UM,
    RB87_MASS_KG,
    get_rydberg_params,
)


@dataclass(frozen=True)
class DecaySweepRow:
    """One row of the isolated Rydberg-decay sweep."""

    gamma_per_us: float
    lifetime_us: float
    numerical_fidelity: float
    numerical_error: float
    analytical_error: float
    analytical_fidelity: float


@dataclass(frozen=True)
class DopplerSweepRow:
    """One row of the isolated Doppler-temperature sweep."""

    temperature_K: float
    temperature_uK: float
    numerical_fidelity: float
    numerical_error: float
    numerical_std: float
    numerical_standard_error: float
    analytical_error: float
    analytical_fidelity: float
    n_samples: int


def _check_positive_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be positive and finite, got {value!r}")
    return value


def _check_nonnegative_finite(name: str, value: float) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be non-negative and finite, got {value!r}")
    return value


def _check_num_points(num_points: int) -> int:
    if not isinstance(num_points, int) or num_points < 2:
        raise ValueError("num_points must be an integer >= 2")
    return num_points


def decay_gamma_grid(
    baseline_gamma: float,
    *,
    num_points: int = 25,
    decades: float = 2.0,
) -> np.ndarray:
    """Return a logarithmic gamma grid centered on ``baseline_gamma``.

    ``decades`` is the total span in orders of magnitude.  The default covers
    baseline/10 through baseline*10, exactly two decades.
    """

    baseline_gamma = _check_positive_finite("baseline_gamma", baseline_gamma)
    num_points = _check_num_points(num_points)
    decades = _check_positive_finite("decades", decades)
    half_span = decades / 2.0
    return np.logspace(
        np.log10(baseline_gamma) - half_span,
        np.log10(baseline_gamma) + half_span,
        num_points,
    )


def doppler_temperature_grid(
    *,
    min_temperature_K: float = 0.1e-6,
    max_temperature_K: float = 100.0e-6,
    num_points: int = 25,
) -> np.ndarray:
    """Return a logarithmic temperature grid for Doppler sweeps."""

    min_temperature_K = _check_positive_finite("min_temperature_K", min_temperature_K)
    max_temperature_K = _check_positive_finite("max_temperature_K", max_temperature_K)
    num_points = _check_num_points(num_points)
    if max_temperature_K <= min_temperature_K:
        raise ValueError("max_temperature_K must be greater than min_temperature_K")
    return np.logspace(np.log10(min_temperature_K), np.log10(max_temperature_K), num_points)


def sweep_decay(
    *,
    omega: float | None = None,
    gammas: Iterable[float] | None = None,
    baseline_gamma: float | None = None,
    num_points: int = 25,
    decades: float = 2.0,
) -> list[DecaySweepRow]:
    """Sweep Rydberg decay rate and compare simulation with the formula."""

    params = get_rydberg_params()
    if omega is None:
        omega = params.omega_rad_per_us
    omega = _check_positive_finite("omega", omega)

    if gammas is None:
        if baseline_gamma is None:
            baseline_gamma = params.rydberg_decay_rate_per_us
        gamma_values = decay_gamma_grid(baseline_gamma, num_points=num_points, decades=decades)
    else:
        gamma_values = np.array(list(gammas), dtype=float)
        if gamma_values.size < 2:
            raise ValueError("gammas must contain at least two points")

    rows: list[DecaySweepRow] = []
    for gamma in gamma_values:
        gamma = _check_positive_finite("gamma", float(gamma))
        result = run_decay_gate(omega=omega, gamma=gamma)
        analytical_error = epsilon_decay_from_gamma(omega, gamma)
        rows.append(
            DecaySweepRow(
                gamma_per_us=gamma,
                lifetime_us=1.0 / gamma,
                numerical_fidelity=result.average_gate_fidelity,
                numerical_error=1.0 - result.average_gate_fidelity,
                analytical_error=analytical_error,
                analytical_fidelity=1.0 - analytical_error,
            )
        )
    return rows


def write_decay_sweep_csv(rows: Iterable[DecaySweepRow], path: str | Path) -> Path:
    """Write decay sweep rows to CSV and return the path."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        raise ValueError("cannot write an empty decay sweep")

    fieldnames = list(asdict(rows[0]).keys())
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return output_path


def read_decay_sweep_csv(path: str | Path) -> list[DecaySweepRow]:
    """Read decay sweep rows from a CSV produced by :func:`write_decay_sweep_csv`."""

    input_path = Path(path)
    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [DecaySweepRow(**{key: float(value) for key, value in row.items()}) for row in reader]
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows


def sweep_doppler(
    *,
    omega: float = DEFAULT_OMEGA_RAD_PER_US,
    k_eff_rad_per_um: float = K_EFF_RAD_PER_UM,
    mass_kg: float = RB87_MASS_KG,
    temperatures_K: Iterable[float] | None = None,
    num_points: int = 25,
    min_temperature_K: float = 0.1e-6,
    max_temperature_K: float = 100.0e-6,
    n_samples: int = 500,
    seed: int | None = 1234,
) -> list[DopplerSweepRow]:
    """Sweep atom temperature and compare Doppler MC with the formula."""

    omega = _check_positive_finite("omega", omega)
    k_eff_rad_per_um = _check_nonnegative_finite("k_eff_rad_per_um", k_eff_rad_per_um)
    mass_kg = _check_positive_finite("mass_kg", mass_kg)
    if not isinstance(n_samples, int) or n_samples < 1:
        raise ValueError("n_samples must be an integer >= 1")

    if temperatures_K is None:
        temperature_values = doppler_temperature_grid(
            min_temperature_K=min_temperature_K,
            max_temperature_K=max_temperature_K,
            num_points=num_points,
        )
    else:
        temperature_values = np.array(list(temperatures_K), dtype=float)
        if temperature_values.size < 2:
            raise ValueError("temperatures_K must contain at least two points")

    rows: list[DopplerSweepRow] = []
    for index, temperature in enumerate(temperature_values):
        temperature = _check_nonnegative_finite("temperature", float(temperature))
        row_seed = None if seed is None else int(seed) + index
        result = run_doppler_gate_mc(
            omega=omega,
            k_eff_rad_per_um=k_eff_rad_per_um,
            temperature_K=temperature,
            mass_kg=mass_kg,
            n_samples=n_samples,
            seed=row_seed,
        )
        analytical_error = epsilon_doppler(k_eff_rad_per_um, temperature, mass_kg, omega)
        rows.append(
            DopplerSweepRow(
                temperature_K=temperature,
                temperature_uK=temperature * 1e6,
                numerical_fidelity=result.average_fidelity,
                numerical_error=1.0 - result.average_fidelity,
                numerical_std=result.std_fidelity,
                numerical_standard_error=result.standard_error,
                analytical_error=analytical_error,
                analytical_fidelity=1.0 - analytical_error,
                n_samples=result.n_samples,
            )
        )
    return rows


def write_doppler_sweep_csv(rows: Iterable[DopplerSweepRow], path: str | Path) -> Path:
    """Write Doppler sweep rows to CSV and return the path."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        raise ValueError("cannot write an empty Doppler sweep")

    fieldnames = list(asdict(rows[0]).keys())
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return output_path


def read_doppler_sweep_csv(path: str | Path) -> list[DopplerSweepRow]:
    """Read Doppler sweep rows from a CSV produced by :func:`write_doppler_sweep_csv`."""

    input_path = Path(path)
    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = []
        for row in reader:
            values = {key: float(value) for key, value in row.items()}
            values["n_samples"] = int(values["n_samples"])
            rows.append(DopplerSweepRow(**values))
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows
