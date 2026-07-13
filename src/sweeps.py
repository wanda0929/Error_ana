"""Parameter sweeps for Rydberg CZ error channels."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from .analytical import epsilon_blockade, epsilon_decay_from_gamma, epsilon_doppler, epsilon_scattering
from .errors.blockade import run_blockade_gate
from .errors.decay import run_decay_gate
from .errors.doppler import run_doppler_gate_mc
from .errors.scattering import run_scattering_gate
from .fidelity import CZ_TARGET, pedersen_fidelity
from .params import DEFAULT_OMEGA_RAD_PER_US, K_EFF_RAD_PER_UM, RB87_MASS_KG, get_rydberg_params


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
class BlockadeSweepRow:
    """One row of the isolated finite-blockade sweep."""

    blockade_to_rabi: float
    blockade_shift: float
    numerical_fidelity: float
    numerical_error: float
    analytical_error: float
    analytical_fidelity: float
    rr_leakage: float
    total_leakage: float
    max_rr_population: float


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


@dataclass(frozen=True)
class ScatteringSweepRow:
    """One row of the isolated intermediate-state-scattering sweep."""

    delta_p_mhz: float
    delta_p_rad_per_us: float
    omega1_over_omega2: float
    numerical_fidelity: float
    numerical_error: float
    analytical_error: float
    analytical_fidelity: float


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


def blockade_ratios(num: int = 30, minimum: float = 5.0, maximum: float = 500.0) -> np.ndarray:
    """Return log-spaced blockade ratios ``U/Ω`` for the standard sweep."""

    num = _check_num_points(num)
    minimum = _check_positive_finite("minimum", minimum)
    maximum = _check_positive_finite("maximum", maximum)
    if maximum <= minimum:
        raise ValueError("maximum must be greater than minimum")
    return np.geomspace(minimum, maximum, num)


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


def scattering_detunings_mhz(
    num: int = 25,
    minimum_mhz: float = 100.0,
    maximum_mhz: float = 100_000.0,
) -> np.ndarray:
    """Return log-spaced intermediate detunings in cycles/us (MHz)."""

    num = _check_num_points(num)
    minimum_mhz = _check_positive_finite("minimum_mhz", minimum_mhz)
    maximum_mhz = _check_positive_finite("maximum_mhz", maximum_mhz)
    if maximum_mhz <= minimum_mhz:
        raise ValueError("maximum_mhz must be greater than minimum_mhz")
    return np.geomspace(minimum_mhz, maximum_mhz, num)


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


def sweep_blockade(
    omega: float | None = None,
    blockade_shifts: Iterable[float] | None = None,
    *,
    ratios: Iterable[float] | None = None,
    n_steps_per_pi: int = 160,
) -> list[BlockadeSweepRow]:
    """Sweep finite-blockade strength and compare simulation with the formula.

    Provide either absolute ``blockade_shifts`` or dimensionless ``ratios``.  If
    neither is provided, the standard 30-point ``U/Ω`` sweep from 5 to 500 is
    used.
    """

    if omega is None:
        omega = get_rydberg_params().omega_rad_per_us
    omega = _check_positive_finite("omega", omega)

    if blockade_shifts is not None and ratios is not None:
        raise ValueError("provide either blockade_shifts or ratios, not both")
    if ratios is None:
        if blockade_shifts is None:
            ratio_values = blockade_ratios()
            shift_values = omega * ratio_values
        else:
            shift_values = np.asarray(list(blockade_shifts), dtype=float)
            ratio_values = shift_values / omega
    else:
        ratio_values = np.asarray(list(ratios), dtype=float)
        shift_values = omega * ratio_values

    if shift_values.ndim != 1 or shift_values.size == 0:
        raise ValueError("sweep requires at least one blockade value")
    if np.any(~np.isfinite(shift_values)) or np.any(shift_values <= 0):
        raise ValueError("all blockade shifts must be positive and finite")
    if np.any(~np.isfinite(ratio_values)) or np.any(ratio_values <= 0):
        raise ValueError("all blockade ratios must be positive and finite")

    rows: list[BlockadeSweepRow] = []
    for ratio, blockade_shift in zip(ratio_values, shift_values):
        result = run_blockade_gate(omega, float(blockade_shift), n_steps_per_pi=n_steps_per_pi)
        fidelity = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
        analytical_error = epsilon_blockade(omega, float(blockade_shift))
        rows.append(
            BlockadeSweepRow(
                blockade_to_rabi=float(ratio),
                blockade_shift=float(blockade_shift),
                numerical_fidelity=float(fidelity),
                numerical_error=float(1.0 - fidelity),
                analytical_error=float(analytical_error),
                analytical_fidelity=float(1.0 - analytical_error),
                rr_leakage=float(result["rr_leakage"]),
                total_leakage=float(result["total_leakage"]),
                max_rr_population=float(result["max_rr_population"]),
            )
        )
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


def sweep_scattering(
    *,
    omega: float | None = None,
    gamma_e: float | None = None,
    detunings_mhz: Iterable[float] | None = None,
    num_points: int = 25,
    omega1_over_omega2: float = 1.0,
) -> list[ScatteringSweepRow]:
    """Sweep intermediate detuning and compare simulation with the formula."""

    params = get_rydberg_params()
    if omega is None:
        omega = params.omega_rad_per_us
    if gamma_e is None:
        gamma_e = params.intermediate_decay_rate_per_us
    omega = _check_positive_finite("omega", omega)
    gamma_e = _check_nonnegative_finite("gamma_e", gamma_e)
    omega1_over_omega2 = _check_positive_finite("omega1_over_omega2", omega1_over_omega2)

    if detunings_mhz is None:
        detuning_values = scattering_detunings_mhz(num=num_points)
    else:
        detuning_values = np.asarray(list(detunings_mhz), dtype=float)
        if detuning_values.ndim != 1 or detuning_values.size < 2:
            raise ValueError("detunings_mhz must contain at least two points")
    if np.any(~np.isfinite(detuning_values)) or np.any(detuning_values <= 0.0):
        raise ValueError("all detunings_mhz values must be positive and finite")

    rows: list[ScatteringSweepRow] = []
    for delta_p_mhz in detuning_values:
        delta_p = 2.0 * np.pi * float(delta_p_mhz)
        result = run_scattering_gate(
            omega=omega,
            gamma_e=gamma_e,
            delta_p=delta_p,
            omega1_over_omega2=omega1_over_omega2,
        )
        analytical_error = epsilon_scattering(gamma_e, delta_p, omega1_over_omega2)
        rows.append(
            ScatteringSweepRow(
                delta_p_mhz=float(delta_p_mhz),
                delta_p_rad_per_us=float(delta_p),
                omega1_over_omega2=float(omega1_over_omega2),
                numerical_fidelity=float(result.average_gate_fidelity),
                numerical_error=float(1.0 - result.average_gate_fidelity),
                analytical_error=float(analytical_error),
                analytical_fidelity=float(1.0 - analytical_error),
            )
        )
    return rows


def _write_dataclass_csv(rows: Iterable[object], path: str | Path, *, empty_message: str) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    if not rows:
        raise ValueError(empty_message)

    fieldnames = list(asdict(rows[0]).keys())
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))
    return output_path


def write_decay_sweep_csv(rows: Iterable[DecaySweepRow], path: str | Path) -> Path:
    """Write decay sweep rows to CSV and return the path."""

    return _write_dataclass_csv(rows, path, empty_message="cannot write an empty decay sweep")


def write_blockade_sweep_csv(rows: Iterable[BlockadeSweepRow], path: str | Path) -> Path:
    """Write blockade sweep rows to CSV and return the path."""

    return _write_dataclass_csv(rows, path, empty_message="cannot write an empty blockade sweep")


def write_doppler_sweep_csv(rows: Iterable[DopplerSweepRow], path: str | Path) -> Path:
    """Write Doppler sweep rows to CSV and return the path."""

    return _write_dataclass_csv(rows, path, empty_message="cannot write an empty Doppler sweep")


def write_scattering_sweep_csv(rows: Iterable[ScatteringSweepRow], path: str | Path) -> Path:
    """Write scattering sweep rows to CSV and return the path."""

    return _write_dataclass_csv(rows, path, empty_message="cannot write an empty scattering sweep")


def read_decay_sweep_csv(path: str | Path) -> list[DecaySweepRow]:
    """Read decay sweep rows from a CSV produced by :func:`write_decay_sweep_csv`."""

    input_path = Path(path)
    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [DecaySweepRow(**{key: float(value) for key, value in row.items()}) for row in reader]
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows


def read_blockade_sweep_csv(path: str | Path) -> list[BlockadeSweepRow]:
    """Read blockade sweep rows from a CSV produced by :func:`write_blockade_sweep_csv`."""

    input_path = Path(path)
    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [BlockadeSweepRow(**{key: float(value) for key, value in row.items()}) for row in reader]
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows


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


def read_scattering_sweep_csv(path: str | Path) -> list[ScatteringSweepRow]:
    """Read scattering sweep rows from a CSV produced by :func:`write_scattering_sweep_csv`."""

    input_path = Path(path)
    with input_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [ScatteringSweepRow(**{key: float(value) for key, value in row.items()}) for row in reader]
    if not rows:
        raise ValueError(f"no rows found in {input_path}")
    return rows
