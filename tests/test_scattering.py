from pathlib import Path

import numpy as np
import pytest

from src.analytical import epsilon_scattering
from src.errors.scattering import (
    run_scattering_gate,
    scattering_rate_per_us,
    single_photon_rabi_frequencies,
)
from src.params import DEFAULT_OMEGA_RAD_PER_US, get_rydberg_params
from src.sweeps import scattering_detunings_mhz, sweep_scattering


def _delta_p(detuning_mhz: float) -> float:
    return 2.0 * np.pi * detuning_mhz


def test_epsilon_scattering_baseline_and_ratio_symmetry():
    params = get_rydberg_params()
    baseline = epsilon_scattering(params.intermediate_decay_rate_per_us, _delta_p(1000.0))
    assert np.isclose(baseline, 1.67e-2, rtol=0.04)
    assert np.isclose(
        epsilon_scattering(params.intermediate_decay_rate_per_us, _delta_p(1000.0), omega1_over_omega2=3.0),
        epsilon_scattering(params.intermediate_decay_rate_per_us, _delta_p(1000.0), omega1_over_omega2=1.0 / 3.0),
    )


def test_single_photon_rabi_frequencies_reproduce_effective_omega_and_rate_symmetry():
    delta_p = _delta_p(1000.0)
    for ratio in (1.0, 3.0, 1.0 / 3.0):
        omega1, omega2 = single_photon_rabi_frequencies(DEFAULT_OMEGA_RAD_PER_US, delta_p, ratio)
        assert np.isclose(omega1 / omega2, ratio)
        assert np.isclose(omega1 * omega2 / (2.0 * delta_p), DEFAULT_OMEGA_RAD_PER_US)

    gamma_e = get_rydberg_params().intermediate_decay_rate_per_us
    assert np.isclose(
        scattering_rate_per_us(DEFAULT_OMEGA_RAD_PER_US, gamma_e, delta_p, 3.0),
        scattering_rate_per_us(DEFAULT_OMEGA_RAD_PER_US, gamma_e, delta_p, 1.0 / 3.0),
    )


def test_scattering_zero_gamma_is_ideal_after_local_z_correction():
    result = run_scattering_gate(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        gamma_e=0.0,
        delta_p=_delta_p(1000.0),
    )
    assert result.average_gate_fidelity == 1.0
    assert result.infidelity < 1e-10
    assert np.isclose(result.raw_average_gate_fidelity, 0.2)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"omega": 0.0, "gamma_e": 0.0, "delta_p": _delta_p(1000.0)}, "omega"),
        ({"omega": DEFAULT_OMEGA_RAD_PER_US, "gamma_e": -1.0, "delta_p": _delta_p(1000.0)}, "gamma"),
        ({"omega": DEFAULT_OMEGA_RAD_PER_US, "gamma_e": 1.0, "delta_p": 0.0}, "delta_p"),
        ({"omega": DEFAULT_OMEGA_RAD_PER_US, "gamma_e": 1.0, "delta_p": _delta_p(1000.0), "omega1_over_omega2": 0.0}, "omega1_over_omega2"),
    ],
)
def test_scattering_rejects_invalid_parameters(kwargs, match):
    with pytest.raises(ValueError, match=match):
        run_scattering_gate(**kwargs)


def test_scattering_baseline_matches_analytical_within_factor_two():
    params = get_rydberg_params()
    result = run_scattering_gate(
        omega=params.omega_rad_per_us,
        gamma_e=params.intermediate_decay_rate_per_us,
        delta_p=_delta_p(1000.0),
    )
    analytical = epsilon_scattering(params.intermediate_decay_rate_per_us, _delta_p(1000.0))
    numerical = result.infidelity
    assert numerical > 0.0
    assert 0.5 < numerical / analytical < 2.0


def test_large_intermediate_detuning_suppresses_scattering():
    params = get_rydberg_params()
    result = run_scattering_gate(
        omega=params.omega_rad_per_us,
        gamma_e=params.intermediate_decay_rate_per_us,
        delta_p=_delta_p(100_000.0),
    )
    assert result.average_gate_fidelity > 0.999
    assert result.infidelity < epsilon_scattering(params.intermediate_decay_rate_per_us, _delta_p(1000.0)) / 50.0


def test_scattering_sweep_shape_range_and_monotonic_error():
    rows = sweep_scattering(detunings_mhz=[1000.0, 10_000.0, 100_000.0])
    assert len(rows) == 3
    assert rows[-1].delta_p_mhz / rows[0].delta_p_mhz == 100.0
    assert all(row.numerical_fidelity <= 1.0 for row in rows)
    assert rows[0].numerical_error > rows[-1].numerical_error
    assert rows[0].analytical_error > rows[-1].analytical_error
    assert len(scattering_detunings_mhz(num=20)) == 20


def test_html_contains_scattering_section_and_figure_reference():
    html = (Path(__file__).resolve().parents[1] / "site" / "index.html").read_text(encoding="utf-8")
    assert 'id="sec-error-scattering"' in html
    assert "Intermediate-State Scattering" in html
    assert "fidelity_vs_detuning.png" in html
    assert "showScatteringFallback" in html
