# Issue #6: Intermediate-state scattering: simulation, analytical formula, sweep, figure, HTML section

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: No intermediate-state scattering error channel exists. No `src/errors/scattering.py`, no `epsilon_scattering()` in `src/analytical.py`, no scattering sweep, no `figures/fidelity_vs_detuning.png`, no "Intermediate-State Scattering" HTML section.
- **User impact / severity**: High for the milestone. This is the fourth error-channel vertical slice. Scattering is the dominant error source at small intermediate detuning and is the key physics motivation for choosing two-photon vs. single-photon excitation schemes.
- **Recommended next step**: Implement the scattering slice end-to-end following the established `decay.py` pattern (Lindblad mesolve, process-tomography extraction, Pedersen fidelity with local-Z correction).
- **Confidence**: high
- **Main blocker, if any**: Issues #1 and #2 are merged. The `params.py` already exposes `intermediate_linewidth_mhz` and `intermediate_decay_rate_per_us` from ARC. No missing infrastructure.

## 2. Affected-path context

### Affected subsystem

The fourth error-channel vertical slice: `src/errors/scattering.py`, `src/analytical.py` (add formula), `src/sweeps.py` (add sweep), scripts, tests, and the "Intermediate-State Scattering" section in `site/index.html`.

### Repo vocabulary / components needed to understand this issue

| Term | Meaning |
|---|---|
| `\|e⟩` = 5P₃/₂ | Intermediate excited state in the two-photon excitation path (780 nm leg). |
| `Γ_e` | Natural linewidth of 5P₃/₂ ≈ 2π × 6.07 MHz. Already in `params.py` as `intermediate_linewidth_mhz`. |
| `Δ_p` | Single-photon detuning from `\|e⟩`. Evered 2023 baseline: ~1 GHz. Sweep parameter. |
| `Ω_eff` | Effective two-photon Rabi frequency = `Ω₁ × Ω₂ / (2Δ_p)`. Project uses `omega_rad_per_us`. |
| `Ω₁` | Single-photon Rabi frequency on the 780 nm (lower/probe) leg. |
| `Ω₂` | Single-photon Rabi frequency on the 480 nm (upper/coupling) leg. |
| `q = Ω₁/Ω₂` | Beam ratio. Experiments report 3–12 for 780+480 nm (780 nm is always stronger due to ~600× larger dipole matrix element). |
| `γ_scatter` | Effective scattering rate while pulse is on. From |g⟩ side: `Γ_e × (Ω₁/(2Δ_p))²`. From |r⟩ side: `Γ_e × (Ω₂/(2Δ_p))²`. Both contribute. |
| Balanced beams | Special case `Ω₁ = Ω₂`. Minimizes scattering for fixed Ω_eff. Requires extreme 480 nm power (~300 mW into 3 µm waist). |
| `mesolve` | QuTiP Lindblad master equation solver. Used for all dissipative channels. |
| `DECAY_BASIS` | 8-state basis `(00, 0g, 0r, g0, gg, gr, r0, rg)` with `\|rr⟩` projected out. Reusable for scattering. |

### Current relevant control/data flow

```
src/params.py
    -> intermediate_linewidth_mhz, intermediate_decay_rate_per_us (from ARC)
    -> src/analytical.py::epsilon_scattering(...)           [new]
    -> src/errors/scattering.py::run_scattering_gate(...)  [new]
    -> src/sweeps.py scattering sweep                      [new]
    -> scripts/run_sweeps.py writes sweep CSV               [extend]
    -> scripts/generate_figures.py writes figure            [extend]
    -> site/index.html embeds figure + explanation          [modify]
```

### Why the key files matter

- **`src/errors/decay.py`** (249 lines): The direct implementation template. Scattering uses the same `DECAY_BASIS`, the same process-tomography via `_propagate_process_outputs`, and the same Choi → Kraus → Pedersen pipeline. The only difference is the collapse operator and that the scattering rate depends on `Δ_p`.
- **`src/params.py`**: Already computes `intermediate_linewidth_mhz = 6.0659` and `intermediate_decay_rate_per_us = 38.113` from ARC's `Rubidium87().getStateLifetime(5, 1, 1.5, ...)`. No new ARC calls needed.
- **`src/fidelity.py`**: Contains `pedersen_fidelity_kraus()` and `LOCAL_Z_PRODUCT` for frame correction. Scattering uses the Kraus path, same as decay.
- **`src/hamiltonian.py`**: No changes needed. Same infinite-blockade Hamiltonian from `decay.py`.
- **`site/index.html`**: HTML template at line ~1100 shows the section pattern. Scattering section would be section 05 or 06 depending on Doppler ordering.

### First file a new implementer should open

`src/errors/decay.py` — the scattering module is structurally identical except for the collapse operator definition and parameterization.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: End-to-end vertical slice for intermediate-state scattering, following the established pattern.
- **Reproduction steps**: N/A (feature request, not bug)
- **Expected behavior**: `src/errors/scattering.py` exists and produces fidelity results consistent with the analytical formula to within factor 2 at Δ_p = 1 GHz.
- **Actual behavior**: Module does not exist.
- **Acceptance criteria from issue**:
  - Δ_p = 100 GHz → F → 1
  - Δ_p = 1 GHz → numerical matches analytical within factor 2
  - Sweep: ≥20 points from ~100 MHz to ~100 GHz
  - Figure with analytical overlay + Evered line
  - HTML section added
  - Tests pass

### Repository evidence

- **`src/params.py` (lines 155-162)**: Already computes intermediate state lifetime from ARC and exposes:
  - `intermediate_lifetime_s = 2.624e-8 s`
  - `intermediate_decay_rate_per_us = 38.113 μs⁻¹`
  - `intermediate_linewidth_mhz = 6.066 MHz` (= Γ_e/2π)
  
  **Confirmed by runtime**: These values are computed from `Rubidium87().getStateLifetime(5, 1, 1.5, temperature=0, ...)`.

- **`src/errors/decay.py`**: Provides a complete template for Lindblad-based error channels:
  - 8-state `DECAY_BASIS` with `|rr⟩` projected out
  - `_hamiltonian(omega, atom_index)`: drives one atom at a time
  - `_collapse_operators(gamma)`: constructs `√γ |g⟩⟨r|` for decay
  - `_propagate_process_outputs()`: full process tomography (16 input states → Choi)
  - `_choi_from_process_outputs()` → `_kraus_from_choi()` → `pedersen_fidelity_kraus()`

- **`src/analytical.py`**: Currently has `epsilon_decay` and `epsilon_blockade`. No scattering formula.

- **`src/sweeps.py`**: Has `sweep_decay()` and `sweep_blockade()` with CSV I/O. Pattern is clear for adding `sweep_scattering()`.

- **`scripts/run_sweeps.py` and `scripts/generate_figures.py`**: Extend with scattering sweep and figure. Pattern established by decay/blockade.

- **Missing**: No `Δ_p` parameter in `RydbergParams`. The intermediate detuning is an experimental choice (not atomic physics), so it should be a function parameter, not a class field — same design as `gamma` in `run_decay_gate()`.

- **Existing tests**: `test_decay.py` shows the test pattern: zero-error boundary, baseline match to analytical, parameter rejection.

## 4. Root-cause candidates

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
|---|---|---|---|---|
| Module not yet implemented (planned feature) | No `src/errors/scattering.py` exists; issue explicitly requests it; milestone tracks it | None — this is the expected state | Check file existence | high |

Single candidate because this is a well-specified feature request, not a bug.

## 5. Decision

- **Chosen likely cause**: Feature not yet implemented. All prerequisites are in place.
- **Recommended fix**: Implement `src/errors/scattering.py` following the `decay.py` pattern with a dephasing Lindblad operator.
- **Rejected alternatives**:
  - *Population-loss model (decay from |r⟩ back to |g⟩ via |e⟩)*: Not correct. Scattering from the intermediate state during a two-photon drive causes decoherence (projects the atom's Rydberg/ground superposition), not population transfer between computational states. The standard model is effective dephasing.
  - *Full 3-level model (include |e⟩ explicitly)*: Overly complex. The adiabatic elimination is valid when Δ_p >> Γ_e, which is satisfied for the entire sweep range (~100 MHz to 100 GHz). Using an effective 2-level model with dephasing is both more efficient and physically clearer.
- **Assumptions**:
  - **Beam ratio parameterized** (default `q = Ω₁/Ω₂ = 1`, balanced). The scattering formula `ε ∝ (q + 1/q)/2` is symmetric and minimized at q=1. Experiments achieve ratios of 3:1 (Zhang 2010) to 12:1 (Gaëtan 2009) due to the ~600× dipole moment asymmetry between 5S→5P and 5P→70S, but the penalty is only 1.67× at ratio 3 — modest, and within the factor-of-2 test tolerance.
  - Adiabatic elimination is valid throughout the sweep range.
  - The collapse operator is dephasing: `L = √(γ_scatter) × |r⟩⟨r|` (pure dephasing in {|g⟩,|r⟩} basis). During a two-photon drive, |e⟩ gets virtual population from **both** |g⟩ (via Ω₁) and |r⟩ (via Ω₂). Both contribute to scattering.
- **Maintainer decision needed**: Choice of dephasing operator. The AGENTS.md says "effective dephasing Lindblad channel." Two natural choices:
  1. `L = √(γ_scatter) |r⟩⟨r|` — dephases the Rydberg state relative to ground. Standard for two-photon scattering.
  2. `L = √(γ_scatter) |g⟩⟨g|` — equivalent physics (up to a Lamb shift) since σ_z = 2|r⟩⟨r| - I.
  
  **Recommended default**: `L = √(γ_scatter) |r⟩⟨r|` per atom during the active pulse on that atom. This is the standard effective operator for off-resonant intermediate-state scattering (de Léséleuc 2018).

## 6. Implementation handoff

### Ordered implementation steps

1. **Add `epsilon_scattering()` to `src/analytical.py`**
   - General formula: `ε_scatter = (7π/8) × (Γ_e / Δ_p) × (q + 1/q) / 2` where `q = Ω₁/Ω₂`
   - Derivation: during a π-pulse, the atom time-averages 50% in |g⟩ and 50% in |r⟩. Virtual |e⟩ population comes from both sides: `(Ω₁/(2Δ_p))²` from |g⟩, `(Ω₂/(2Δ_p))²` from |r⟩. The time-averaged scattering rate is `Γ_e(Ω₁² + Ω₂²)/(8Δ_p²)`. Multiply by pulse duration `π/Ω_eff = 2πΔ_p/(Ω₁Ω₂)` to get scatter per π-pulse: `πΓ_e/(4Δ_p) × (q + 1/q)`. Sum over 7/2 effective π-pulses for the CZ gate → `(7π/8) × (Γ_e/Δ_p) × (q + 1/q)/2`.
   - For balanced beams (q=1): reduces to `7πΓ_e/(8Δ_p)` (independent of Ω, minimum)
   - At ratio q=3 (or 1/3): penalty factor = (3 + 1/3)/2 = 1.67×
   - Input: `gamma_e` (rad/μs), `delta_p` (rad/μs), `omega1_over_omega2` (default=1.0, balanced)
   - Validation: at Δ_p = 2π × 1000 MHz, balanced: ε ≈ 0.017; ratio=3: ε ≈ 0.028

2. **Create `src/errors/scattering.py`**
   - Reuse `DECAY_BASIS` and Hamiltonian infrastructure from `decay.py`
   - **Two** dephasing collapse operators per driven atom: one for virtual |e⟩ from |g⟩ side, one from |r⟩ side:
     - `L_g = √(γ_g) × |r⟩⟨r|` with `γ_g = Γ_e × (Ω₁/(2Δ_p))²`
     - `L_r = √(γ_r) × |r⟩⟨r|` with `γ_r = Γ_e × (Ω₂/(2Δ_p))²`
     - (Or equivalently, one operator with combined rate `γ_total = Γ_e(Ω₁² + Ω₂²)/(4Δ_p²)`)
   - Active only during the pulse on that specific atom. Atom 1's scattering: steps 1 and 3; atom 2's: step 2.
   - Interface: `run_scattering_gate(omega, gamma_e, delta_p, omega1_over_omega2=1.0) -> ScatteringGateResult`
   - Internally compute: `Ω₁ = sqrt(q × 2Ω_effΔ_p)`, `Ω₂ = sqrt(2Ω_effΔ_p/q)` where `q = omega1_over_omega2` and `Ω_eff = omega`
   - Then `γ_total = gamma_e × (Ω₁² + Ω₂²) / (4 × delta_p²)`
   - Use same Choi → Kraus → Pedersen pipeline as `decay.py`

3. **Add `sweep_scattering()` to `src/sweeps.py`**
   - Sweep `delta_p` logarithmically from 2π × 100 MHz to 2π × 100 GHz (≥20 points)
   - Return dataclass rows with: `delta_p_mhz`, `numerical_fidelity`, `numerical_error`, `analytical_error`, `analytical_fidelity`

4. **Add `ScatteringSweepRow` dataclass and CSV I/O** to `src/sweeps.py`

5. **Extend `scripts/run_sweeps.py`** — add scattering sweep call, print summary

6. **Extend `scripts/generate_figures.py`** — produce `figures/fidelity_vs_detuning.png`
   - x-axis: Δ_p in GHz (log scale)
   - y-axis: gate infidelity (log scale recommended given ~3 orders of magnitude range)
   - Analytical curve + numerical scatter + Evered vertical line at 1 GHz
   - Follow styling from `plot_blockade_sweep()`

7. **Create `tests/test_scattering.py`**
   - Zero-error boundary: delta_p → ∞ (use 2π × 100000 MHz) → F ≈ 1
   - Baseline match: delta_p = 2π × 1000 MHz → numerical ε matches analytical within factor 2
   - Parameter rejection: negative/zero delta_p, negative gamma_e

8. **Add HTML section** to `site/index.html`
   - Section 05 (or adjust numbering if Doppler lands first)
   - Physics: explain two-photon excitation, why intermediate state exists, what "scattering" means here
   - Collapsible derivation of the 7πΓ_e/(8Δ_p) formula
   - Embedded figure
   - Contrast with single-photon schemes (Yb) where this error is absent
   - Takeaway: "This error dominates when intermediate detuning is small"

### Files likely to change

| File | Change type |
|---|---|
| `src/errors/scattering.py` | **Create** — main simulation module |
| `src/errors/__init__.py` | Maybe add import (currently just a docstring) |
| `src/analytical.py` | Add `epsilon_scattering()` function |
| `src/sweeps.py` | Add `ScatteringSweepRow`, `sweep_scattering()`, CSV I/O |
| `scripts/run_sweeps.py` | Add scattering sweep call |
| `scripts/generate_figures.py` | Add scattering figure generation |
| `tests/test_scattering.py` | **Create** — boundary + match + rejection tests |
| `site/index.html` | Add scattering section (~80 lines) |

### Scope boundaries / what not to touch unless new evidence appears

- Do NOT modify `src/hamiltonian.py` — the existing Hamiltonians suffice.
- Do NOT modify `src/fidelity.py` — the Kraus fidelity path is already there.
- Do NOT add `delta_p` to `RydbergParams` — it's an experimental choice parameter, not an atomic constant. Pass it as a function argument.
- Do NOT implement a full 3-level model with explicit |e⟩ — adiabatic elimination is physically correct and computationally cleaner.

### Regression test shape

```python
def test_scattering_zero_error_boundary():
    """At very large Δ_p, scattering vanishes and fidelity → 1."""
    result = run_scattering_gate(omega=DEFAULT_OMEGA, gamma_e=38.1, delta_p=2*pi*100000)
    assert result.average_gate_fidelity > 1 - 1e-8

def test_scattering_baseline_matches_analytical():
    """At Δ_p = 1 GHz, numerical ε agrees with 7πΓ_e/(8Δ_p) within factor 2."""
    result = run_scattering_gate(omega=DEFAULT_OMEGA, gamma_e=38.1, delta_p=2*pi*1000)
    analytical = epsilon_scattering(gamma_e=38.1, delta_p=2*pi*1000)
    numerical = 1 - result.average_gate_fidelity
    assert 0.5 < numerical / analytical < 2.0

def test_scattering_rejects_invalid_params():
    with pytest.raises(ValueError):
        run_scattering_gate(omega=DEFAULT_OMEGA, gamma_e=-1, delta_p=2*pi*1000)
```

### Validation command(s)

```bash
pytest tests/test_scattering.py -v
python scripts/run_sweeps.py  # should include scattering sweep output
python scripts/generate_figures.py  # should produce fidelity_vs_detuning.png
```

### Expected pre-fix signal

- `ImportError: cannot import name 'run_scattering_gate' from 'src.errors.scattering'`
- `figures/fidelity_vs_detuning.png` does not exist

### Expected post-fix signal

- All three tests pass
- `figures/fidelity_vs_detuning.png` exists with correct content
- Analytical and numerical curves show 1/Δ_p scaling
- At Δ_p = 1 GHz, balanced: infidelity ≈ 1.7% (= 7πΓ_e/(8Δ_p)). At ratio 3: ≈ 2.8%.

### Negative/edge validation

- At Δ_p = 100 MHz (very close to resonance), formula may break down — adiabatic elimination requires Δ_p >> Γ_e ≈ 6 MHz. At 100 MHz this ratio is ~16, still OK but approaching the limit. If numerical diverges from analytical here, that's expected physics, not a bug.
- The balanced formula `7πΓ_e/(8Δ_p)` being Ω-independent is non-obvious — verify this holds in the simulation by checking two different Ω values give the same ε at the same Δ_p.
- **Ratio symmetry check**: run at q=3 and q=1/3 and verify identical infidelity. This tests that both scattering channels (from |g⟩ side and |r⟩ side) are correctly modeled.

## 7. Risks, edge cases, rollback

### Risks

- **Collapse operator choice matters**: If using `|g⟩⟨g|` instead of `|r⟩⟨r|`, the resulting dephasing is the same (both project onto the σ_z eigenbasis), but the Lamb shift differs. For the pure dephasing model, these are equivalent up to a unitary rotation that gets absorbed into the local-Z correction. Should still verify numerically.
- **Beam ratio**: The general formula `(7π/8)(Γ_e/Δ_p)(q+1/q)/2` includes the beam ratio `q = Ω₁/Ω₂`. The function `q + 1/q` is symmetric: ratios 3:1 and 1:3 give identical scattering. Balanced (q=1) is the minimum. Experimental ratios of 3 penalize by only 1.67×. The Lindblad simulation should reproduce this symmetry as a consistency check.
- **Scattering during non-driven periods**: If atom 1 is in |r⟩ during step 2, is there still scattering on atom 1? No — the scattering comes from virtual population in |e⟩, which only exists while the drive laser is on. Collapse ops should only be active during the relevant pulse.

### Edge cases

- Δ_p very small (< 10 × Γ_e): adiabatic elimination breaks down. Document this limit in the sweep but don't try to fix it — it's outside the physically relevant regime.
- Δ_p = 0: must be rejected (division by zero in γ_scatter).
- γ_e = 0: should return F = 1 (no scattering from a stable intermediate state).

### Rollback/revert notes

Module is self-contained. If scattering code is wrong, deleting `src/errors/scattering.py` and reverting the additions to `analytical.py`, `sweeps.py`, scripts, and HTML restores the previous state. No existing modules are modified in their core logic.

## 8. Open questions

### Q1: Collapse operator — dephasing on driven atom only or on both atoms throughout?

- **What was found**: The AGENTS.md says "effective dephasing Lindblad channel with rate γ_scatter = Γ_e × (Ω₁/(2Δ_p))² during the pulse." The key phrase is "during the pulse" — scattering only occurs while the drive is on.
- **Why it is ambiguous**: The decay channel in `decay.py` applies collapse operators throughout all three pulse steps to both atoms, because spontaneous emission from |r⟩ happens regardless of whether a pulse is on. Scattering from |e⟩ is different — it requires the drive laser to create virtual |e⟩ population.
- **Options**:
  - A) Collapse ops on both atoms at all times (simpler code, overestimates error by ~40%)
  - B) Collapse ops only on the driven atom during its pulse (physically correct, requires per-segment collapse operator lists)
- **Recommended default**: Option B. The per-segment structure already exists in `decay.py`'s `_run_three_pulse_sequence`. Simply pass different c_ops to each `_mesolve_final` call.

### Q2: Does the analytical formula use Γ_e as angular frequency or cycles?

- **What was found**: ARC returns `intermediate_decay_rate_per_us = 38.113 μs⁻¹` which is `Γ_e` in angular frequency (= 2π × 6.066 MHz). The formula `(7π/8) × (Γ_e/Δ_p) × (q + 1/q)/2` works when Γ_e and Δ_p are in the same angular units; `q` is dimensionless.
- **Why it is ambiguous**: Literature sometimes quotes Γ_e/2π in MHz and Δ_p in GHz, leading to factor-of-2π confusion.
- **Options**: Use angular units throughout (matching `params.py` convention)
- **Recommended default**: Express everything in rad/μs. `epsilon_scattering(gamma_e_rad_per_us, delta_p_rad_per_us, omega1_over_omega2=1.0)` with documentation stating units.

### Q3: HTML section numbering — 05 or 06?

- **What was found**: Current HTML has sections 01 (intro), 02 (ideal), 03 (decay), 04 (blockade). Doppler (issue #5) and scattering (this issue) are both pending.
- **Why it is ambiguous**: The AGENTS.md lists error sources as: decay, blockade, Doppler, scattering, amplitude. If Doppler lands first, scattering is section 06. If scattering lands first, it takes section 05.
- **Options**: Number it 05 and let Doppler be 06, or match AGENTS.md ordering.
- **Recommended default**: Follow AGENTS.md order — Doppler = 05, Scattering = 06. If Doppler hasn't landed, use section 06 anyway with a placeholder note. This avoids merge conflicts.
