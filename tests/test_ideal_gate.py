import numpy as np

from src.fidelity import CZ_TARGET, pedersen_fidelity
from src.params import DEFAULT_OMEGA_RAD_PER_US
from src.protocol import run_ideal_gate


def test_ideal_gate_fidelity_with_local_z_correction():
    result = run_ideal_gate(omega=DEFAULT_OMEGA_RAD_PER_US, n_steps_per_pi=120)
    fidelity = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
    assert 1.0 - fidelity < 1e-10


def test_ideal_gate_raw_phases_and_entangling_phase():
    result = run_ideal_gate(omega=DEFAULT_OMEGA_RAD_PER_US, n_steps_per_pi=120)
    np.testing.assert_allclose(result["unitary_diagonal"], [1, -1, -1, -1], atol=1e-10)
    assert abs(abs(result["entangling_phase"]) - np.pi) < 1e-8


def test_population_traces_start_and_end_in_computational_space():
    result = run_ideal_gate(omega=DEFAULT_OMEGA_RAD_PER_US, n_steps_per_pi=80)
    times = result["times"]
    assert np.isclose(times[0], 0.0)
    assert np.isclose(times[-1], 4.0 * result["t_pi"])
    for population in result["populations"].values():
        assert np.isclose(population[0], 0.0, atol=1e-12)
        assert np.isclose(population[-1], 0.0, atol=1e-10)
        assert np.all(population >= -1e-12)
        assert np.all(population <= 1.0 + 1e-10)
