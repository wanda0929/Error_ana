import numpy as np
import pytest

from src.analytical import epsilon_blockade
from src.errors.blockade import run_blockade_gate
from src.fidelity import CZ_TARGET, pedersen_fidelity
from src.params import get_rydberg_params
from src.sweeps import blockade_ratios, sweep_blockade


def _blockade_infidelity(result):
    return 1.0 - pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)


def test_epsilon_blockade_formula():
    assert np.isclose(epsilon_blockade(1.0, 10.0), 1.0 / 800.0)
    with pytest.raises(ValueError):
        epsilon_blockade(1.0, 0.0)
    with pytest.raises(ValueError):
        epsilon_blockade(1.0, -2.0)


def test_blockade_large_u_is_ideal_after_phase_calibration():
    """A much stronger-than-baseline blockade should be numerically ideal."""

    omega = 25.13
    result = run_blockade_gate(omega, omega * 5000.0, n_steps_per_pi=80)
    assert _blockade_infidelity(result) < 1e-8
    assert result["max_rr_population"] < 5e-8


def test_blockade_baseline_matches_analytical():
    params = get_rydberg_params()
    result = run_blockade_gate(
        params.omega_rad_per_us,
        params.blockade_shift_rad_per_us,
        n_steps_per_pi=120,
    )
    numerical_error = _blockade_infidelity(result)
    analytical_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)
    assert abs(numerical_error - analytical_error) / analytical_error < 0.5
    assert result["rr_leakage"] >= 0.0
    assert result["raw_unitary"].shape == (4, 4)


def test_blockade_zero_or_negative_blockade_raises():
    omega = 25.13
    with pytest.raises(ValueError):
        run_blockade_gate(omega, 0.0, n_steps_per_pi=20)
    with pytest.raises(ValueError):
        run_blockade_gate(omega, -omega, n_steps_per_pi=20)


def test_blockade_sweep_shape_and_scaling():
    omega = 25.13
    ratios = blockade_ratios(num=24, minimum=5.0, maximum=500.0)
    rows = sweep_blockade(omega, ratios=ratios, n_steps_per_pi=40)
    assert len(rows) == 24
    assert np.isclose(rows[0].blockade_to_rabi, 5.0)
    assert np.isclose(rows[-1].blockade_to_rabi, 500.0)
    assert rows[0].numerical_error > rows[-1].numerical_error
    for row in rows[8:]:
        rel = abs(row.numerical_error - row.analytical_error) / row.analytical_error
        assert rel < 0.1
