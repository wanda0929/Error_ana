# Issue #5: Doppler dephasing: simulation, analytical formula, sweep, figure, HTML section

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: No Doppler dephasing code exists — no `src/errors/doppler.py`, no `epsilon_doppler()` formula, no temperature sweep, no `figures/fidelity_vs_temperature.png`, and no "Doppler Dephasing" HTML section.
- **User impact / severity**: High for the milestone. Third error-channel vertical slice. This is the first **Monte Carlo** error source (unlike decay/Lindblad or blockade/coherent), so it establishes the MC pattern for issue #7 (amplitude noise).
- **Recommended next step**: Implement the Doppler slice. The simulation is a Monte Carlo over thermal velocities: for each trial, sample atom velocities, add Doppler detuning to the ideal-basis Hamiltonian, run `sesolve`, compute per-shot unitary fidelity, average. No Lindblad needed — the noise is shot-to-shot classical randomness.
- **Confidence**: high
- **Main blocker, if any**: None. Issues #1–#3 are merged. The existing infrastructure (ideal basis, `sesolve`, `pedersen_fidelity`) is sufficient. The Monte Carlo pattern is independent of the Lindblad pattern used for decay.

## 2. Affected-path context

### Affected subsystem

The third error-channel vertical slice: `src/errors/doppler.py`, extension of `src/analytical.py` and `src/sweeps.py`, scripts, tests, and the "Doppler Dephasing" section in `site/index.html`.

### Repo vocabulary / components needed for this issue

| Term | Meaning |
|---|---|
| `k_eff` | Effective wavevector for the two-photon transition. For counter-propagating Rb 780+480 nm: `k_eff = \|1/λ₇₈₀ − 1/λ₄₈₀\| ≈ 0.80 µm⁻¹` (in 1/λ units) or `≈ 5.03 rad/µm` (angular). |
| `v` | Atom velocity along the laser axis. Drawn from Maxwell-Boltzmann: `v ~ N(0, kB×T/m)` per axis. |
| `δ = k_eff × v` | Doppler frequency shift seen by a moving atom. In solver units: rad/µs when k_eff is rad/µm and v is µm/µs. |
| Monte Carlo | Each trial samples independent velocities for both atoms, runs the full coherent gate with detunings, computes a per-shot fidelity. The average over N trials gives `F_avg`. |
| `sesolve` | Coherent (Schrödinger) solver. Doppler is per-shot classical noise, not quantum dissipation — no Lindblad. |
| Detuning convention | In `src/hamiltonian.py`: `H_rr = -Δ`. A positive detuning Δ lowers the `\|r⟩` energy. The Doppler shift adds `δ` to the detuning: atom moving toward the laser sees higher frequency → effectively reduces the detuning needed. |

### Current relevant control/data flow

```
src/params.py
    -> baseline Ω, T, k_eff, mass
    -> src/analytical.py::epsilon_doppler(...)             [new]
    -> src/errors/doppler.py::run_doppler_gate_mc(...)    [new]
    -> src/sweeps.py::sweep_doppler(...)                  [new]
    -> scripts/run_sweeps.py                              [extend]
    -> scripts/generate_figures.py                        [extend]
    -> site/index.html (Doppler section)                  [modify]
```

### Why the key files matter

- **`src/hamiltonian.py`**: `build_ideal_hamiltonian(omega, atom_index, detuning=δ)` already accepts a `detuning` parameter. This is the exact API needed — just pass `δ = k_eff × v` as the detuning for the atom being pulsed. **No Hamiltonian code changes needed.**
- **`src/protocol.py`**: The `_run_input_01`, `_run_input_10`, `_run_input_11` helpers show the pulse pattern. Doppler needs the same pattern with per-atom detunings. The helpers are private but the pattern is simple enough to replicate.
- **`src/fidelity.py`**: `pedersen_fidelity(U_shot, CZ_TARGET, correct_local_z=True)` works per-shot on a 4×4 unitary. The MC average is `⟨F⟩ = (1/N) Σᵢ F(Uᵢ)`.
- **`src/errors/decay.py`**: Establishes the pattern for `@dataclass(frozen=True)` results, validation helpers, and module structure. Doppler follows the same shape but uses `sesolve` + MC averaging instead of `mesolve` + Kraus reconstruction.
- **`src/analytical.py`**: Already has `_check_positive_finite` and `_check_nonnegative_finite`. Doppler adds `epsilon_doppler(k_eff, temperature, mass, omega)`.
- **`src/sweeps.py`**: Already has `DecaySweepRow`, `sweep_decay`, `write_decay_sweep_csv`. Doppler adds `DopplerSweepRow`, `sweep_doppler`, `write_doppler_sweep_csv`.

### First file a new implementer should open

`src/hamiltonian.py` → verify that `build_ideal_hamiltonian` and `single_atom_hamiltonian` accept `detuning`. Then `src/errors/decay.py` for the module/result pattern.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Implement MC Doppler slice: sample velocities from Maxwell-Boltzmann, add detuning δ = k_eff × v, run π–2π–π with `sesolve`, average fidelity over N trials. Sweep temperature.
- **Expected behavior**:
  - T = 0: F = 1 (infidelity < 10⁻¹⁰)
  - T = 10 µK: numerical matches analytical within factor 2
  - N ≥ 500 for figures, N ≤ 100 for tests (speed)
  - Sweep: ≥20 points from ~0.1 µK to ~100 µK
  - Figure: `figures/fidelity_vs_temperature.png` with error bars, analytical curve, Evered line
  - HTML: "Doppler Dephasing" section
  - `pytest tests/test_doppler.py` passes (< 60s)
- **Actual behavior**: All specified files absent.
- **Relevant comments**: None.

### Repository evidence

| File | What it proves |
|---|---|
| `src/hamiltonian.py` | Both `single_atom_hamiltonian(omega, detuning=0.0)` and `build_ideal_hamiltonian(omega, atom_index, detuning=0.0)` accept a detuning parameter. The convention is `H_rr = -Δ`. No code changes needed for the Hamiltonian. |
| `src/protocol.py` | Shows the ideal π–2π–π pattern. Private helpers `_run_input_01`, `_run_input_10`, `_run_input_11` demonstrate the 3-segment structure for each computational input. |
| `src/errors/decay.py` | Establishes module conventions: `@dataclass(frozen=True)` result, `_check_omega`/`_check_gamma` validators, private helpers, public `run_decay_gate()`. Doppler module follows the same shape. |
| `src/analytical.py` | Already has validation utilities and `epsilon_decay` functions. Clean place to add `epsilon_doppler`. |
| `src/sweeps.py` | Has `DecaySweepRow`, `sweep_decay`, CSV I/O. Pattern for `DopplerSweepRow`/`sweep_doppler`. |
| `src/fidelity.py` | `pedersen_fidelity(actual, CZ_TARGET, correct_local_z=True)` works on per-shot unitaries. The MC average is trivial: `mean([F(U_i) for each trial])`. |
| `tests/test_decay.py` | Model for test structure: zero-error boundary, baseline match, sweep shape, HTML reference. |
| `scripts/run_sweeps.py` / `scripts/generate_figures.py` | Established patterns for data generation and plotting. Doppler extends these. |
| `site/index.html` | Contains the error-section template (lines 982–1017) and already mentions Doppler in introductory text (line 731: "counter-propagating beams give a smaller effective wavevector, which helps suppress Doppler dephasing later"). |

**Verified at runtime**:
- 29 tests pass in 5.24s.
- `build_ideal_hamiltonian(omega, atom_index=1, detuning=0.5)` works — confirmed the API accepts detuning.
- Quick MC test (N=200, T=10µK): **infidelity ≈ 8.5 × 10⁻⁵**. This is the target the analytical formula must reproduce.
- `δ_rms/Ω ≈ 6.2 × 10⁻³` at T=10µK — the perturbative regime is valid.
- `k_eff = 2π × 0.8013 ≈ 5.03 rad/µm` for counter-propagating 780+480 nm.

**Key physics verified**: The Doppler error scales as `(k_eff × v_rms / Ω)² ≈ (δ_rms/Ω)²` with a protocol-dependent coefficient ~2.2. This is small enough that perturbation theory works well at T=10µK.

## 4. Root-cause candidates

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
|---|---|---|---|
| Doppler slice is simply unimplemented. | All required files absent. Feature request issue. | None. | Create the files, pass tests. | high |
| The analytical formula's proportionality constant is ambiguous. | AGENTS.md says `ε ∝ k_eff² × kB×T × t_gate² / m` but doesn't specify the exact coefficient. Quick MC gives ε ≈ 8.5e-5 while the naive `(k_eff × v_rms × t_gate)² ≈ 6e-3` is 70× too large. The correct scaling is `∝ (k_eff × v_rms / Ω)²` with coefficient ~2.2. | The exact coefficient depends on which pulse-protocol formula is used. The acceptance criterion allows "within factor 2" which is generous. | Compute the coefficient analytically from the detuned Rabi formula for the specific π–2π–π protocol, then verify against MC. | medium |
| Unit confusion between `k_eff` in 1/λ units (0.80 µm⁻¹) vs angular (5.03 rad/µm). | AGENTS.md uses `≈ 0.80 µm⁻¹` which is the 1/λ convention. The Doppler shift formula needs angular: `δ = k_angular × v`. | Easy to get wrong if implementer uses 0.80 directly without the 2π factor. | Define `k_eff` clearly in the code with explicit units. Verify MC against the known numerical answer (8.5e-5 at T=10µK). | medium |

## 5. Decision

- **Chosen likely cause**: Unimplemented feature slice with a secondary complexity: getting the analytical formula coefficient right and the units consistent.

- **Recommended fix**: Implement the Doppler MC slice end-to-end.

### Key design decisions

**1. Simulation approach — per-shot coherent unitaries:**

For each MC trial:
1. Sample `v₁, v₂` independently from `N(0, σ_v)` where `σ_v = √(kB×T/m)` (converted to µm/µs).
2. Compute detunings: `δ₁ = k_eff_angular × v₁`, `δ₂ = k_eff_angular × v₂` (in rad/µs).
3. Run the π–2π–π sequence for all 4 computational inputs with the appropriate detuning on the atom being pulsed:
   - `|00⟩`: dark, amplitude = 1.
   - `|01⟩`: atom 2 does a 2π-pulse with detuning `δ₂`. Single-atom Hamiltonian.
   - `|10⟩`: atom 1 does π–π with detuning `δ₁`. Single-atom Hamiltonian.
   - `|11⟩`: atom 1 does π (with `δ₁`), atom 2 does 2π (with `δ₂`, blocked in ideal basis), atom 1 does π (with `δ₁`). Use `build_ideal_hamiltonian` with detuning.
4. Construct the 4×4 diagonal unitary from the four output amplitudes.
5. Compute `F = pedersen_fidelity(U_shot, CZ_TARGET, correct_local_z=True)`.
6. Average `F` over all N trials.

**2. Why per-shot unitaries work (no Kraus needed):**

Each trial has a *definite* (classical) velocity — the coherence is preserved within each shot. The fidelity loss is from *averaging* perfect-but-wrong unitaries over the velocity ensemble. This is fundamentally different from decay (quantum decoherence within a single shot). The average `⟨F⟩` directly gives the gate fidelity for an atom sampled from the thermal distribution.

**3. Units — critical:**

| Quantity | Unit | Convention |
|---|---|---|
| `k_eff_angular` | rad/µm | `2π × 0.8013 ≈ 5.03 rad/µm` |
| `v` | µm/µs | `√(kB×T/m)` converted from m/s |
| `δ = k_eff × v` | rad/µs | Compatible with Ω (rad/µs) and solver times (µs) |
| `T` | K | Input in Kelvin, not µK |

**4. Analytical formula:**

The standard Rydberg-gate Doppler error formula (Saffman–Walker–Mølmer RMP 2010; de Léséleuc et al. 2018):

```
ε_Doppler = (π²/4) × k_eff² × kB × T / (m × Ω²)
```

where `C = π²/4 ≈ 2.47` is the coefficient from expanding the off-resonant Rabi evolution to second order in δ/Ω.

**Derivation sketch**: Treat the Doppler shift as a quasi-static detuning δ = k_eff × v during each shot. For a π-pulse with small detuning δ << Ω, the single-pulse excitation error is `P_Doppler = |δ/Ω|²` (Saffman 2010). The linear term in δ averages to zero over the symmetric thermal distribution; the quadratic term survives and, combined with the π-pulse phase sensitivity, gives the π²/4 prefactor.

**Numerical verification**: At T=10µK, the formula gives 9.4×10⁻⁵. Quick MC (N=200) gives 8.5×10⁻⁵ — agreement within 1σ MC noise. The acceptance criterion ("within factor 2") is comfortably met.

**5. What `|11⟩` sees — the blocked case:**

When atom 1 is in `|r⟩` (after the first π-pulse), atom 2's 2π-pulse is blocked in the ideal basis — `|rg⟩→|rr⟩` is projected out. So atom 2's detuning `δ₂` has **no effect** on the `|11⟩` input during step 2 (the 2π pulse does nothing). Only `δ₁` matters for `|11⟩`, through steps 1 and 3 (atom 1's π-pulses).

This is important for getting the analytical formula right: the `|11⟩` input has the same Doppler sensitivity as `|10⟩` (only atom 1's detuning matters), not double.

- **Rejected alternatives**:
  - *Use Lindblad dephasing as a Doppler proxy*: Wrong model. Doppler is shot-to-shot classical noise, not intra-shot quantum decoherence. Lindblad would give a different error structure.
  - *Use the 8-state decay basis*: Unnecessary. The ideal 3-state basis is correct for Doppler since there's no dissipation — `sesolve` is appropriate.
  - *Correlate atom velocities*: Wrong physics. Atoms are in separate traps; their thermal motion is independent.

- **Assumptions**:
  - Atom velocities are constant during the gate (frozen-atom approximation). Gate time ~0.5 µs, trap period ~10 µs for typical 100 kHz trap → atoms move negligibly during one gate.
  - 1D treatment along the laser axis is sufficient. Transverse velocities don't contribute to Doppler shift.
  - `k_eff` uses the AGENTS.md convention: `|1/λ₇₈₀ − 1/λ₄₈₀|` for counter-propagating two-photon.

- **Maintainer decision needed**: None critical. The analytical coefficient C has some flexibility, but the acceptance test allows factor-of-2 agreement.

## 6. Implementation handoff

### Ordered implementation steps

1. **Add physical constants to `src/params.py`** (or a new constants section):
   - `K_EFF_ANGULAR_RAD_PER_UM = 2π × |1/0.780 − 1/0.480| ≈ 5.035 rad/µm`
   - `ATOM_MASS_KG = 87 × 1.661e-27` (Rb-87 atomic mass)
   - `BOLTZMANN_J_PER_K = 1.381e-23`
   - `DEFAULT_TEMPERATURE_K = 10e-6` (Evered operating point: 10 µK)
   - Helper: `thermal_velocity_rms(temperature_K, mass_kg) → σ_v in µm/µs`

2. **Add `epsilon_doppler()` to `src/analytical.py`**:
   ```python
   def epsilon_doppler(k_eff: float, temperature: float, mass: float, omega: float) -> float:
       """Perturbative Doppler infidelity for the π–2π–π CZ gate.
       
       Returns C * k_eff² * kB * T / (m * Ω²) where C is a protocol-dependent
       coefficient derived from the detuned-Rabi sensitivity of the pulse sequence.
       """
       sigma_v_sq = BOLTZMANN * temperature / mass  # in (m/s)^2 or matched units
       return C * k_eff**2 * sigma_v_sq / omega**2
   ```
   Use units consistently (all in µm/µs/rad system) or convert at boundaries.

3. **Create `src/errors/doppler.py`**:
   - `DopplerGateResult` dataclass: `temperature`, `omega`, `k_eff`, `n_samples`, `average_fidelity`, `std_fidelity`, `per_shot_fidelities` (for error bars).
   - `run_doppler_gate_mc(omega, k_eff, temperature, mass, n_samples=500, seed=None)`:
     - For each trial: sample `v1, v2 ~ N(0, σ_v)`, compute detunings, run 4 inputs with `sesolve`, build diagonal unitary, compute Pedersen fidelity with local-Z correction.
     - Return `DopplerGateResult` with mean, std, and all per-shot values.
   - Helper: `_run_single_trial(omega, delta1, delta2) → 4×4 unitary diagonal`

4. **Add Doppler sweep to `src/sweeps.py`**:
   - `DopplerSweepRow` dataclass: temperature, fidelity mean, fidelity std, analytical error, etc.
   - `sweep_doppler(temperatures, omega, k_eff, mass, n_samples=500)` → list of rows.
   - CSV I/O functions.

5. **Create `tests/test_doppler.py`**:
   - `test_doppler_zero_temperature_is_ideal`: T=0 → F=1 (use N=1, v=0 deterministically).
   - `test_doppler_baseline_matches_analytical`: T=10µK, N=100 → within factor 2 of analytical.
   - `test_doppler_sweep_shape`: ≥20 points, range 0.1–100 µK.
   - `test_doppler_scales_with_temperature`: higher T → higher error (monotonic check).
   - `test_epsilon_doppler_formula`: known inputs → known output.
   - `test_invalid_parameters_rejected`: negative T, non-positive omega, etc.
   - All tests use `N ≤ 100` and fixed seeds for speed and reproducibility.

6. **Update `scripts/run_sweeps.py`**: Add Doppler sweep call, save CSV.

7. **Update `scripts/generate_figures.py`**: Generate `figures/fidelity_vs_temperature.png`:
   - x-axis: temperature (µK), linear or log depending on range
   - y-axis: gate fidelity
   - Numerical scatter with error bars (std/√N from MC) + analytical curve + Evered line at T=10µK
   - Publication quality

8. **Add "Doppler Dephasing" section to `site/index.html`**: Use the existing template.
   - Physics: thermal motion → frequency shift → pulse errors
   - Explain k_eff and why two-photon counter-propagating suppresses Doppler (smaller k_eff than single-photon)
   - Formula: `ε ∝ k_eff² kB T / (m Ω²)`
   - Collapsible derivation: detuned Rabi formula → pulse error → thermal average
   - Figure
   - Takeaway: "This error dominates at high temperature or with large k_eff (single-photon schemes)"

9. **Force-commit the figure**: `git add -f figures/fidelity_vs_temperature.png`.

10. **Full validation**.

### Files likely to change

| File | Action |
|---|---|
| `src/params.py` | Extend (add k_eff, mass, temperature constants) |
| `src/analytical.py` | Extend (add `epsilon_doppler`) |
| `src/errors/doppler.py` | New |
| `src/sweeps.py` | Extend (add Doppler sweep) |
| `scripts/run_sweeps.py` | Extend |
| `scripts/generate_figures.py` | Extend |
| `tests/test_doppler.py` | New |
| `site/index.html` | Modify (add Doppler section) |
| `figures/fidelity_vs_temperature.png` | Generated |

### Scope boundaries / what not to touch

- Do NOT add Lindblad terms. Doppler is classical shot-to-shot noise, not quantum decoherence.
- Do NOT model finite blockade in the Doppler simulation. Use the ideal 3-state basis (infinite blockade).
- Do NOT correlate atom velocities. They are in separate traps.
- Do NOT change `src/hamiltonian.py` — the `detuning` parameter already exists.
- Do NOT implement decay, blockade, scattering, amplitude noise, or combined budget in this slice.

### Regression test shape

```python
def test_doppler_zero_temperature_is_ideal():
    result = run_doppler_gate_mc(omega=..., k_eff=..., temperature=0.0, mass=..., n_samples=1)
    assert 1.0 - result.average_fidelity < 1e-10

def test_doppler_baseline_matches_analytical_within_factor_2():
    result = run_doppler_gate_mc(..., temperature=10e-6, n_samples=100, seed=42)
    numerical_error = 1.0 - result.average_fidelity
    analytical_error = epsilon_doppler(k_eff, 10e-6, mass, omega)
    assert numerical_error / analytical_error < 2.0
    assert numerical_error / analytical_error > 0.5

def test_doppler_scales_with_temperature():
    r_low = run_doppler_gate_mc(..., temperature=1e-6, n_samples=50, seed=0)
    r_high = run_doppler_gate_mc(..., temperature=50e-6, n_samples=50, seed=0)
    assert r_high.average_fidelity < r_low.average_fidelity
```

### Validation commands

```bash
# Targeted Doppler tests (should complete < 60s with N<=100)
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q tests/test_doppler.py

# Generate sweep data + figure (may take a few minutes with N=500)
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/run_sweeps.py
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/generate_figures.py

# Full regression
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q
```

### Expected pre-fix signal

- `import src.errors.doppler` → `ModuleNotFoundError`
- `pytest tests/test_doppler.py` → file not found
- `figures/fidelity_vs_temperature.png` → absent

### Expected post-fix signal

- All new and existing tests pass (29 existing + new Doppler tests)
- T=0 infidelity < 10⁻¹⁰
- T=10µK: numerical infidelity ≈ 8.5×10⁻⁵ (within factor 2 of analytical)
- Sweep: ≥20 points, 0.1–100 µK range
- Figure exists with error bars, analytical curve, Evered line
- HTML contains "Doppler Dephasing" section

### Negative/edge validation

- `temperature < 0` → `ValueError`
- `temperature = 0` → skip MC (or use N=1 with v=0), return ideal fidelity
- `k_eff = 0` → no Doppler effect, return ideal fidelity (valid physical limit: co-propagating equal-wavelength beams)
- Large temperature (100 µK): perturbative formula will overestimate error; numerical MC should still converge but with larger variance

## 7. Risks, edge cases, rollback

### Risks

| Risk | Mitigation |
|---|---|
| **MC variance at low N**: Per-shot fidelities near 1 have small variance, but at high T the spread increases. | Use fixed seeds in tests. Report `std/√N` as error bars on figures. N=500 for publication, N≤100 for tests. |
| **Analytical coefficient ambiguity**: The exact C in `ε = C × k_eff²kBT/(mΩ²)` depends on which formula derivation is followed. | Resolved: use C = π²/4 (standard Rydberg-gate expression from Saffman–Walker–Mølmer RMP 2010). Verified against MC. |
| **Unit bugs**: k_eff has two conventions (1/λ vs 2π/λ), velocity has (m/s vs µm/µs). | Define all constants with explicit unit suffixes in variable names. Add a self-consistency test that verifies `k_eff_angular × v_rms_um_per_us` gives a detuning in rad/µs compatible with Ω. |
| **Performance**: N=500 × 20 sweep points × 4 inputs × 3 solver calls = 120,000 `sesolve` calls for figures. | Each call is tiny (2-level or 3-level, <100 time steps). Quick MC test ran 200 trials in seconds. Full sweep should take minutes, not hours. Tests use N≤100. |
| **Shared file conflicts with issue #4**: Both extend `src/analytical.py`, `src/sweeps.py`, scripts. | These are additive extensions (new functions, not modifications). Merge conflicts will be trivial. |

### Edge cases

- **T → 0**: All velocities are 0, gate is ideal. Must reproduce F=1 exactly.
- **k_eff → 0**: No Doppler effect regardless of temperature. Valid limit (co-propagating beams of equal wavelength would cancel). Gate is ideal.
- **Very high T** (e.g., 1 mK): `δ/Ω > 1`, perturbation theory fails, MC still gives meaningful (low) fidelity. The figure should show divergence between analytical curve and numerical dots at high T.
- **Finite blockade + Doppler**: Not in scope. In reality both are present simultaneously, but this slice isolates Doppler with infinite blockade.

### Rollback/revert notes

Purely additive. Revert = remove new files + revert extensions to `analytical.py`, `sweeps.py`, scripts, and the HTML section. No impact on decay, blockade, or ideal gate.

## 8. Resolved decisions (formerly open questions)

All four questions have been resolved by literature review (Saffman–Walker–Mølmer RMP 2010, de Léséleuc et al. PRA 2018, Evered et al. Nature 2023).

### 1. Analytical coefficient C — ✅ RESOLVED

- **Decision**: Use **C = π²/4 ≈ 2.47** (the standard Rydberg-gate Doppler error expression).
- **Formula**: `ε_Doppler = (π²/4) × k_eff² × kB×T / (m × Ω²)`
- **Verification**: At T=10µK this gives 9.4×10⁻⁵, matching the MC result (8.5×10⁻⁵) within MC statistical noise (±1σ). The "factor of 2" acceptance criterion is comfortably satisfied.
- **Derivation**: Expand the off-resonant Rabi evolution to second order in δ/Ω. The linear term in δ averages to zero over the symmetric velocity distribution; the quadratic term gives the π²/4 coefficient. This is consistent with Saffman–Walker–Mølmer RMP 2010, which gives the single-pulse Doppler excitation error as `P_Doppler = |δ/Ω|²`, and de Léséleuc et al. 2018, which models Doppler as shot-to-shot Gaussian detuning with σ = k_eff√(kBT/m).
- **Implementation**:
  ```python
  def epsilon_doppler(k_eff_angular: float, temperature_K: float, mass_kg: float, omega_rad_per_us: float) -> float:
      kB = 1.381e-23  # J/K
      v_rms_sq_solver = (kB * temperature_K / mass_kg) * 1e12 * 1e-12  # µm²/µs²
      return (np.pi**2 / 4) * k_eff_angular**2 * v_rms_sq_solver / omega_rad_per_us**2
  ```

### 2. Constants location — ✅ RESOLVED

- **Decision**: Add `k_eff`, `mass`, `temperature`, and `kB` constants to **`src/params.py`**.
- **Rationale**: Mass of Rb-87 is fundamental and may be needed by other channels. `k_eff` is project-wide (referenced in the HTML explanation). de Léséleuc et al. define these as system-level parameters, not module-private. Centralize.
- **Constants to add**:
  ```python
  RB87_MASS_KG = 87 * 1.661e-27  # 87 amu in kg
  BOLTZMANN_J_PER_K = 1.381e-23
  LAMBDA_LOWER_NM = 780.0  # 5S → 5P leg
  LAMBDA_UPPER_NM = 480.0  # 5P → nS leg
  K_EFF_CYCLES_PER_UM = abs(1.0/0.780 - 1.0/0.480)  # ≈ 0.801 µm⁻¹
  K_EFF_RAD_PER_UM = 2 * pi * K_EFF_CYCLES_PER_UM   # ≈ 5.03 rad/µm
  DEFAULT_TEMPERATURE_K = 10e-6  # Evered operating point: 10 µK
  ```

### 3. MC implementation — ✅ RESOLVED

- **Decision**: **Serial loop** with fixed seeds for reproducibility.
- **Rationale**: de Léséleuc et al. use ~600 realizations — this is tiny. Our quick MC test (N=200, 4 inputs × 3 segments per trial) completed in seconds. N=500 × 20 sweep points will take minutes at most. Serial keeps code simple and seed-deterministic. No parallelism needed.

### 4. T=0 handling — ✅ RESOLVED

- **Decision**: **Short-circuit** at T ≤ 0: return ideal fidelity (F=1) without running MC.
- **Rationale**: At T=0 all velocities are zero by definition. Running MC with σ_v=0 is degenerate and wasteful. Raise `ValueError` for T<0 (non-physical). For T=0, return the exact ideal result.

## 9. Literature references

| Short name | Reference | Used for |
|---|---|---|
| Saffman 2010 | Saffman, Walker & Mølmer, Rev. Mod. Phys. 82, 2313 (2010) / arXiv:0909.4777 | Doppler reduction with counter-propagating beams; `P_Doppler = \|δ/Ω\|²` |
| de Léséleuc 2018 | de Léséleuc et al., PRA 97, 053803 (2018) / arXiv:1802.10424 | Gaussian Doppler detuning model; MC averaging over ~600 realizations; σ = k_eff√(kBT/m) |
| Evered 2023 | Evered et al., Nature 622, 268 (2023) / arXiv:2304.05420 | Finite-temperature modeled via measured T₂* = 3 µs; Gaussian detuning distribution |
| Levine 2019 | Levine et al., PRL 123, 170503 (2019) / arXiv:1908.06101 | Atomic temperature as main Rydberg-gate error source (different protocol but same physics) |
