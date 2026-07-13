from pathlib import Path

import numpy as np
import pytest

from src.analytical import epsilon_total_additive
from src.errors.combined import (
    DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
    COMBINED_BASIS,
    evaluate_error_budget,
    run_combined_gate,
)
from src.errors.decay import run_decay_gate
from src.params import DEFAULT_OMEGA_RAD_PER_US, RB87_MASS_KG, get_rydberg_params


def test_combined_basis_contains_dark_and_rr_states():
    assert COMBINED_BASIS == ("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg", "rr")


def test_epsilon_total_additive_includes_all_five_channels():
    params = get_rydberg_params()
    total, components = epsilon_total_additive(
        omega=params.omega_rad_per_us,
        gamma=params.rydberg_decay_rate_per_us,
        blockade_shift=params.blockade_shift_rad_per_us,
        gamma_e=params.intermediate_decay_rate_per_us,
        delta_p=DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
        sigma_omega=0.02,
        k_eff_rad_per_um=0.0,
        temperature_K=0.0,
        mass_kg=RB87_MASS_KG,
    )
    assert set(components) == {"decay", "blockade", "doppler", "scattering", "amplitude"}
    assert np.isclose(total, sum(components.values()))
    assert components["decay"] > components["blockade"]
    assert components["scattering"] > components["amplitude"]


def test_combined_all_zero_errors_is_ideal():
    result = run_combined_gate(
        omega=DEFAULT_OMEGA_RAD_PER_US,
        gamma=0.0,
        blockade_shift=float("inf"),
        gamma_e=0.0,
        delta_p=DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
        sigma_omega=0.0,
        k_eff_rad_per_um=0.0,
        temperature_K=0.0,
        mass_kg=RB87_MASS_KG,
        n_samples=1,
        seed=42,
    )
    assert 1.0 - result.average_fidelity < 1e-8
    assert result.additive_sum == 0.0
    assert np.all(result.per_shot_fidelities <= 1.0 + 1e-12)


def test_combined_decay_only_matches_isolated_decay():
    params = get_rydberg_params()
    isolated = run_decay_gate(params.omega_rad_per_us, params.rydberg_decay_rate_per_us)
    combined = run_combined_gate(
        omega=params.omega_rad_per_us,
        gamma=params.rydberg_decay_rate_per_us,
        blockade_shift=float("inf"),
        gamma_e=0.0,
        delta_p=DEFAULT_INTERMEDIATE_DETUNING_RAD_PER_US,
        sigma_omega=0.0,
        k_eff_rad_per_um=0.0,
        temperature_K=0.0,
        mass_kg=RB87_MASS_KG,
        n_samples=1,
        seed=42,
    )
    assert combined.infidelity > 0.0
    assert abs(combined.infidelity - isolated.infidelity) / isolated.infidelity < 0.10


def test_combined_budget_additive_approximation_at_baseline():
    budget = evaluate_error_budget(n_samples=2, seed=1, individual_n_samples=1)
    numerical = budget.combined_numerical_error
    additive = budget.additive_numerical_sum
    assert 0.005 < numerical < 0.03
    assert abs(additive - numerical) / numerical < 0.20


def test_html_contains_combined_and_summary_sections():
    html = (Path(__file__).resolve().parents[1] / "site" / "index.html").read_text(encoding="utf-8")
    assert 'id="sec-combined-budget"' in html
    assert 'id="sec-summary-table"' in html
    assert "error_budget_combined.png" in html
    assert "showCombinedFallback" in html
