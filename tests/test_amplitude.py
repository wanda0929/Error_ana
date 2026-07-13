from pathlib import Path
import math

import numpy as np
import pytest

from src.analytical import AMPLITUDE_PROTOCOL_SENSITIVITY, epsilon_amplitude
from src.errors.amplitude import (
    DEFAULT_SIGMA_OMEGA,
    projected_gate_for_rabi_scale,
    run_amplitude_noise_gate,
)
from src.fidelity import CZ_TARGET, LOCAL_Z_PRODUCT, pedersen_fidelity
from src.params import DEFAULT_OMEGA_RAD_PER_US
from src.sweeps import amplitude_sigma_grid, sweep_amplitude


def test_epsilon_amplitude_formula_current_baseline():
    expected = (math.pi**2 / 2.0) * DEFAULT_SIGMA_OMEGA**2
    assert np.isclose(AMPLITUDE_PROTOCOL_SENSITIVITY, math.pi**2 / 2.0)
    assert np.isclose(epsilon_amplitude(DEFAULT_SIGMA_OMEGA), expected)
    assert np.isclose(expected, 1.97e-3, rtol=0.01)


def test_resonant_rabi_scale_matches_ideal_gate():
    gate = projected_gate_for_rabi_scale(DEFAULT_OMEGA_RAD_PER_US, rabi_scale=1.0)
    np.testing.assert_allclose(np.diag(gate), [1.0, -1.0, -1.0, -1.0], atol=1e-10)
    assert pedersen_fidelity(LOCAL_Z_PRODUCT @ gate, CZ_TARGET) == 1.0


def test_small_deterministic_rabi_error_matches_sensitivity():
    eps = 1e-3
    gate = projected_gate_for_rabi_scale(DEFAULT_OMEGA_RAD_PER_US, rabi_scale=1.0 + eps)
    numerical = 1.0 - pedersen_fidelity(LOCAL_Z_PRODUCT @ gate, CZ_TARGET)
    assert np.isclose(numerical, epsilon_amplitude(eps), rtol=2e-5)


def test_amplitude_zero_noise_is_ideal_after_local_z_correction():
    result = run_amplitude_noise_gate(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        sigma_omega=0.0,
        n_samples=5,
        seed=123,
    )
    assert result.average_fidelity == 1.0
    assert result.infidelity < 1e-10
    assert np.all(result.per_shot_fidelities == 1.0)
    assert np.all(result.fractional_rabi_errors == 0.0)


def test_amplitude_baseline_matches_analytical_within_factor_two():
    result = run_amplitude_noise_gate(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        sigma_omega=DEFAULT_SIGMA_OMEGA,
        n_samples=120,
        seed=42,
    )
    numerical_error = result.infidelity
    analytical_error = epsilon_amplitude(DEFAULT_SIGMA_OMEGA)
    assert numerical_error > 0.0
    assert 0.5 < numerical_error / analytical_error < 2.0
    assert result.standard_error > 0.0


def test_amplitude_error_scales_quadratically():
    low = run_amplitude_noise_gate(sigma_omega=0.02, n_samples=120, seed=7)
    high = run_amplitude_noise_gate(sigma_omega=0.04, n_samples=120, seed=7)
    ratio = high.infidelity / low.infidelity
    assert 3.5 < ratio < 4.5


def test_amplitude_sweep_shape_range_and_error_growth():
    rows = sweep_amplitude(num_points=20, n_samples=2, seed=11)
    assert len(rows) == 20
    assert np.isclose(rows[0].sigma_omega, 0.001)
    assert np.isclose(rows[-1].sigma_omega, 0.20)
    assert rows[-1].analytical_error > rows[0].analytical_error
    assert rows[-1].numerical_error > rows[0].numerical_error
    assert all(row.n_samples == 2 for row in rows)

    grid = amplitude_sigma_grid(num_points=20)
    assert len(grid) == 20
    assert np.all(np.diff(grid) > 0.0)


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda: epsilon_amplitude(-0.01), "sigma_omega"),
        (lambda: epsilon_amplitude(0.01, sensitivity=0.0), "sensitivity"),
        (lambda: projected_gate_for_rabi_scale(0.0, 1.0), "omega"),
        (lambda: projected_gate_for_rabi_scale(DEFAULT_OMEGA_RAD_PER_US, math.inf), "rabi_scale"),
        (lambda: run_amplitude_noise_gate(sigma_omega=-0.01), "sigma_omega"),
        (lambda: run_amplitude_noise_gate(n_samples=0), "n_samples"),
    ],
)
def test_amplitude_rejects_invalid_parameters(call, match):
    with pytest.raises(ValueError, match=match):
        call()


def test_html_contains_amplitude_section_and_figure_reference():
    html = (Path(__file__).resolve().parents[1] / "site" / "index.html").read_text(encoding="utf-8")
    assert 'id="sec-error-amplitude"' in html
    assert "Amplitude Noise" in html
    assert "fidelity_vs_amplitude_noise.png" in html
    assert "showAmplitudeFallback" in html
