import numpy as np
import pytest

from src.params import get_rydberg_params


def test_arc_parameters_for_rb87_70s_are_in_expected_range():
    params = get_rydberg_params()
    assert np.isclose(params.rydberg_lifetime_us, 373.9, rtol=0.01)
    assert 850.0 < params.c6_abs_ghz_um6 < 890.0
    assert params.distance_um == 3.0
    assert params.blockade_to_rabi > 250.0
    assert np.isclose(params.intermediate_linewidth_mhz, 6.07, rtol=0.03)


def test_invalid_quantum_numbers_fail_before_arc_lookup():
    with pytest.raises(ValueError, match="l must satisfy"):
        get_rydberg_params(n=3, l=3)
