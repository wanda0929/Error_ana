from pathlib import Path

import numpy as np
import pytest

from src.analytical import epsilon_decay, epsilon_decay_from_gamma, rydberg_decay_exposure_time
from src.errors.decay import COMPUTATIONAL_BASIS, run_decay_gate
from src.params import DEFAULT_OMEGA_RAD_PER_US, get_rydberg_params
from src.sweeps import sweep_decay


def test_epsilon_decay_current_baseline():
    params = get_rydberg_params()
    eps = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)
    assert np.isclose(eps, 5.85e-4, rtol=0.03)
    assert np.isclose(eps, epsilon_decay_from_gamma(params.omega_rad_per_us, params.rydberg_decay_rate_per_us))
    assert np.isclose(rydberg_decay_exposure_time(params.omega_rad_per_us), 7.0 * np.pi / (4.0 * params.omega_rad_per_us))


def test_decay_zero_gamma_is_ideal_after_local_z_correction():
    result = run_decay_gate(omega=DEFAULT_OMEGA_RAD_PER_US, gamma=0.0)
    assert result.average_gate_fidelity == 1.0
    assert 1.0 - result.average_gate_fidelity < 1e-10
    assert np.isclose(result.raw_average_gate_fidelity, 0.2)
    assert set(result.final_density_matrices) == set(COMPUTATIONAL_BASIS)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"omega": 0.0, "gamma": 0.0}, "omega"),
        ({"omega": DEFAULT_OMEGA_RAD_PER_US, "gamma": -1e-3}, "gamma"),
    ],
)
def test_decay_rejects_invalid_parameters(kwargs, match):
    with pytest.raises(ValueError, match=match):
        run_decay_gate(**kwargs)


def test_decay_baseline_matches_analytical_within_50_percent():
    params = get_rydberg_params()
    result = run_decay_gate(params.omega_rad_per_us, params.rydberg_decay_rate_per_us)
    numerical_error = 1.0 - result.average_gate_fidelity
    analytical_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)
    assert numerical_error > 0.0
    assert abs(numerical_error - analytical_error) / analytical_error < 0.5


def test_large_decay_rate_lowers_fidelity_significantly():
    result = run_decay_gate(omega=DEFAULT_OMEGA_RAD_PER_US, gamma=0.05)
    assert result.average_gate_fidelity < 0.995
    assert result.infidelity > 10 * epsilon_decay_from_gamma(DEFAULT_OMEGA_RAD_PER_US, 0.002)


def test_decay_sweep_shape_range_and_monotonic_error():
    params = get_rydberg_params()
    rows = sweep_decay(baseline_gamma=params.rydberg_decay_rate_per_us, num_points=20, decades=2.0)
    assert len(rows) == 20
    assert rows[-1].gamma_per_us / rows[0].gamma_per_us >= 99.9
    assert all(row.numerical_fidelity <= 1.0 for row in rows)
    assert rows[-1].numerical_error > rows[0].numerical_error
    assert rows[-1].analytical_error > rows[0].analytical_error


def test_html_contains_decay_section_and_figure_reference():
    html = (Path(__file__).resolve().parents[1] / "site" / "index.html").read_text(encoding="utf-8")
    assert 'id="sec-error-decay"' in html
    assert "Rydberg Decay" in html
    assert "fidelity_vs_decay_rate.png" in html
    assert "showDecayFallback" in html
