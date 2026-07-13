from pathlib import Path

import numpy as np
import pytest

from src.analytical import DOPPLER_PROTOCOL_COEFFICIENT, epsilon_doppler
from src.errors.doppler import projected_gate_for_detunings, run_doppler_gate_mc
from src.fidelity import CZ_TARGET, pedersen_fidelity
from src.params import (
    BOLTZMANN_J_PER_K,
    DEFAULT_OMEGA_RAD_PER_US,
    DEFAULT_TEMPERATURE_K,
    K_EFF_RAD_PER_UM,
    RB87_MASS_KG,
    thermal_velocity_rms_um_per_us,
)
from src.sweeps import doppler_temperature_grid, sweep_doppler


def test_epsilon_doppler_formula_current_baseline():
    eps = epsilon_doppler(K_EFF_RAD_PER_UM, DEFAULT_TEMPERATURE_K, RB87_MASS_KG, DEFAULT_OMEGA_RAD_PER_US)
    expected = (
        DOPPLER_PROTOCOL_COEFFICIENT
        * K_EFF_RAD_PER_UM**2
        * BOLTZMANN_J_PER_K
        * DEFAULT_TEMPERATURE_K
        / (RB87_MASS_KG * DEFAULT_OMEGA_RAD_PER_US**2)
    )
    assert np.isclose(eps, expected)
    assert np.isclose(eps, 9.45e-5, rtol=0.03)


def test_thermal_velocity_units_match_solver_units():
    sigma_v = thermal_velocity_rms_um_per_us(DEFAULT_TEMPERATURE_K, RB87_MASS_KG)
    assert np.isclose(sigma_v, 0.0309, rtol=0.02)
    assert np.isclose(K_EFF_RAD_PER_UM * sigma_v / DEFAULT_OMEGA_RAD_PER_US, 6.2e-3, rtol=0.03)


def test_zero_temperature_is_ideal_after_local_z_correction():
    result = run_doppler_gate_mc(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        k_eff_rad_per_um=K_EFF_RAD_PER_UM,
        temperature_K=0.0,
        mass_kg=RB87_MASS_KG,
        n_samples=3,
        seed=123,
    )
    assert result.average_fidelity == 1.0
    assert result.infidelity < 1e-10
    assert np.all(result.per_shot_fidelities == 1.0)


def test_zero_wavevector_is_ideal_even_at_finite_temperature():
    result = run_doppler_gate_mc(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        k_eff_rad_per_um=0.0,
        temperature_K=50e-6,
        mass_kg=RB87_MASS_KG,
        n_samples=4,
        seed=123,
    )
    assert result.average_fidelity == 1.0


def test_resonant_single_shot_matches_ideal_gate():
    gate = projected_gate_for_detunings(DEFAULT_OMEGA_RAD_PER_US, delta1=0.0, delta2=0.0)
    np.testing.assert_allclose(np.diag(gate), [1.0, -1.0, -1.0, -1.0], atol=1e-10)
    assert pedersen_fidelity(gate, CZ_TARGET, correct_local_z=True) == 1.0


def test_doppler_baseline_matches_analytical_within_factor_2():
    result = run_doppler_gate_mc(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        k_eff_rad_per_um=K_EFF_RAD_PER_UM,
        temperature_K=DEFAULT_TEMPERATURE_K,
        mass_kg=RB87_MASS_KG,
        n_samples=100,
        seed=42,
    )
    numerical_error = result.infidelity
    analytical_error = epsilon_doppler(K_EFF_RAD_PER_UM, DEFAULT_TEMPERATURE_K, RB87_MASS_KG, DEFAULT_OMEGA_RAD_PER_US)
    assert numerical_error > 0.0
    assert 0.5 < numerical_error / analytical_error < 2.0


def test_doppler_error_increases_with_temperature():
    low = run_doppler_gate_mc(
        temperature_K=1e-6,
        n_samples=50,
        seed=7,
    )
    high = run_doppler_gate_mc(
        temperature_K=50e-6,
        n_samples=50,
        seed=7,
    )
    assert high.average_fidelity < low.average_fidelity
    assert high.infidelity > 10.0 * low.infidelity


def test_doppler_sweep_shape_range_and_error_growth():
    rows = sweep_doppler(num_points=20, n_samples=2, seed=11)
    assert len(rows) == 20
    assert np.isclose(rows[0].temperature_uK, 0.1)
    assert np.isclose(rows[-1].temperature_uK, 100.0)
    assert rows[-1].analytical_error > rows[0].analytical_error
    assert all(row.n_samples == 2 for row in rows)

    grid = doppler_temperature_grid(num_points=20)
    assert len(grid) == 20
    assert np.all(np.diff(grid) > 0.0)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"omega": 0.0}, "omega"),
        ({"k_eff_rad_per_um": -1.0}, "k_eff"),
        ({"temperature_K": -1e-6}, "temperature"),
        ({"mass_kg": 0.0}, "mass"),
        ({"n_samples": 0}, "n_samples"),
    ],
)
def test_doppler_rejects_invalid_parameters(kwargs, match):
    with pytest.raises(ValueError, match=match):
        run_doppler_gate_mc(**kwargs)


def test_html_contains_doppler_section_and_figure_reference():
    html = (Path(__file__).resolve().parents[1] / "site" / "index.html").read_text(encoding="utf-8")
    assert 'id="sec-error-doppler"' in html
    assert "Doppler Dephasing" in html
    assert "fidelity_vs_temperature.png" in html
    assert "showDopplerFallback" in html
