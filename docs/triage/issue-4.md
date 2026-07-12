# Issue #4: Finite blockade: simulation, analytical formula, sweep, figure, HTML section

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: The repository has no finite-blockade error channel code. No `src/errors/blockade.py`, no `epsilon_blockade()` analytical formula, no blockade sweep, no `figures/fidelity_vs_blockade.png`, and no "Finite Blockade" HTML section. This is a feature-implementation issue, not a bug report.
- **User impact / severity**: High for the milestone. This is the second error-channel vertical slice (after decay, issue #3). It is the simplest error channel — **purely coherent**, no Lindblad dynamics — so it sets a clean baseline for the implementation pattern.
- **Recommended next step**: Implement the blockade slice end-to-end. The Hamiltonian infrastructure (`build_blockade_hamiltonian` in `src/hamiltonian.py`) already exists, so the core simulation is a thin wrapper around `sesolve` in the 4-state `BLOCKADE_BASIS`. The analytical formula is one line. The primary design question is how to extract a 4×4 computational unitary from the 4-state blockade-basis simulation where `|0⟩` is dark.
- **Confidence**: high
- **Main blocker, if any**: Issues #1 and #2 are merged. Issue #3 (decay) is open but is *not* a dependency — these error channels are designed for parallel development. The only shared infrastructure that issue #3 might create (`src/errors/__init__.py`, `src/analytical.py`, `src/sweeps.py`) is trivial to create here if #3 hasn't landed yet.

## 2. Affected-path context

### Affected subsystem

The second error-channel vertical slice: `src/errors/blockade.py`, `src/analytical.py`, `src/sweeps.py`, plotting/sweep scripts, blockade tests, and the "Finite Blockade" section in `site/index.html`.

### Repo vocabulary / components needed to understand this issue

| Term | Meaning |
|---|---|
| `\|0⟩` | Dark qubit state. Does not couple to the Rydberg laser. |
| `\|g⟩` | Qubit `\|1⟩`. Ground state that couples to `\|r⟩`. |
| `\|r⟩` | Auxiliary Rydberg state (`\|70S₁/₂⟩`). |
| `\|rr⟩` | Doubly-excited Rydberg state. Energy-shifted by U (blockade interaction). |
| `IDEAL_BASIS` | `(\|gg⟩, \|gr⟩, \|rg⟩)` — 3 states, `\|rr⟩` projected out. |
| `BLOCKADE_BASIS` | `(\|gg⟩, \|gr⟩, \|rg⟩, \|rr⟩)` — 4 states, `\|rr⟩` included at energy U. |
| `U` | Rydberg–Rydberg interaction strength. `U = C₆/R⁶`. In solver units: `blockade_shift_rad_per_us`. |
| `Ω` | Two-photon Rabi frequency. Solver units: `omega_rad_per_us`. |
| `U/Ω` | Blockade ratio. Baseline ≈ 296. Larger = better blockade. |
| `ε_block` | Blockade infidelity: `Ω²/(8U²) = 1/(8(U/Ω)²)`. |
| `sesolve` | QuTiP closed-system (Schrödinger equation) solver. No dissipation. |
| `LOCAL_Z_PRODUCT` | `diag(1, -1, -1, 1)`. Corrects the raw `diag(1, -1, -1, -1)` gate to canonical CZ. |

### Current relevant control/data flow

```
src/params.py
    -> baseline Ω, U (from ARC C₆ and R)
    -> src/analytical.py::epsilon_blockade(omega, U)         [new]
    -> src/errors/blockade.py::run_blockade_gate(...)        [new]
    -> src/sweeps.py blockade sweep                          [new]
    -> scripts/run_sweeps.py writes tabular sweep data       [new/extend]
    -> scripts/generate_figures.py writes fidelity_vs_blockade.png  [new/extend]
    -> site/index.html embeds figure + explanation           [modify]
```

### Why the key files matter

- **`src/hamiltonian.py`**: Already contains `build_blockade_hamiltonian(omega, blockade_shift, atom_index)` returning a 4×4 `Qobj` in `BLOCKADE_BASIS`. This is the core Hamiltonian builder — the blockade module wraps it, it does not rewrite it.
- **`src/protocol.py`**: Establishes the ideal π–2π–π pattern using `sesolve`. The blockade module follows the same 3-segment structure but uses `build_blockade_hamiltonian` instead of `build_ideal_hamiltonian`.
- **`src/fidelity.py`**: Contains `pedersen_fidelity(actual, target, correct_local_z=True)` and `apply_local_z_correction()`. The blockade error is coherent, so the result is a unitary — no Kraus operators needed. `pedersen_fidelity` works directly.
- **`site/index.html`**: Lines 981–1017 contain the `ERROR-SOURCE SECTION TEMPLATE` with slots for physics, formula, derivation, figure, and takeaway. The blockade section plugs into this template.
- **`scripts/run_ideal.py`**: Establishes script conventions (sys.path, Agg backend, save to `figures/`, print report).

### First file a new implementer should open

`src/hamiltonian.py` → `build_blockade_hamiltonian()`. Then `src/protocol.py` to see the pulse-sequence pattern. Then `src/fidelity.py` for the unitary fidelity path.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Implement end-to-end vertical slice for finite blockade: coherent simulation, analytical formula `ε_block = Ω²/(8U²)`, parameter sweep over U/Ω, figure, HTML section, tests.
- **Reproduction steps, if provided**: N/A — feature request.
- **Expected behavior**:
  - U/Ω = 1000 → infidelity < 10⁻⁸
  - Evered baseline U/Ω: numerical matches analytical within 50%
  - Sweep: ≥20 points, U/Ω from ~5 to ~500
  - Figure: `figures/fidelity_vs_blockade.png` with numerical dots, analytical curve, Evered line
  - HTML: "Finite Blockade" section with physics, formula, derivation, figure, takeaway
  - `pytest tests/test_blockade.py` passes
- **Actual behavior**: All specified files are absent.
- **Environment/version/config**: `qutip==5.0.4`, `ARC-Alkali-Rydberg-Calculator==3.9.0` in `requirements.txt`. Verified working conda env at `/opt/homebrew/Caskroom/miniconda/base/envs/qutip`.
- **Relevant comments or corrections**: No issue comments.

### Repository evidence

| File | What it proves |
|---|---|
| `src/hamiltonian.py` | `build_blockade_hamiltonian(omega, blockade_shift, atom_index)` exists and returns a 4×4 `Qobj` in `BLOCKADE_BASIS = ("gg", "gr", "rg", "rr")`. The `|rr⟩` state sits at energy `blockade_shift`. Both atom-1 and atom-2 pulse couplings include the `|gr⟩↔|rr⟩` and `|rg⟩↔|rr⟩` terms. |
| `tests/test_hamiltonian.py` | `test_blockade_hamiltonian_matrix_elements_atom2` verifies the 4×4 matrix elements for atom 2 with a known omega and blockade shift. |
| `src/protocol.py` | The ideal gate uses `sesolve` with 3 segments (π, 2π, π). Each segment evolves a state and stitches population traces. The raw unitary is `diag(1, -1, -1, -1)`. |
| `src/fidelity.py` | `pedersen_fidelity(actual, CZ_TARGET, correct_local_z=True)` works on unitary matrices. `LOCAL_Z_PRODUCT = diag(1, -1, -1, 1)` and `apply_local_z_correction()` are tested. |
| `src/params.py` | `get_rydberg_params()` returns `blockade_shift_rad_per_us ≈ 7435.74`, `blockade_to_rabi ≈ 295.86`, `omega_rad_per_us ≈ 25.13`. Verified at runtime. |
| `tests/test_params.py` | Asserts `blockade_to_rabi > 250.0` and C₆ in range `850–890 GHz·µm⁶`. |
| `site/index.html` | Contains `ERROR-SOURCE SECTION TEMPLATE FOR ISSUES #3-#7` at line 982 with the exact structure to use. |
| `scripts/run_ideal.py` | Establishes script conventions: add ROOT to sys.path, Agg backend, save PNG under `figures/`, constrained_layout, 220 dpi, publication-quality style. |
| `.gitignore` | `figures/*.png` is ignored. Per issue #3 triage decisions, publication figures should be force-committed with `git add -f`. |

**Files absent**: `src/errors/` directory, `src/analytical.py`, `src/sweeps.py`, `scripts/run_sweeps.py`, `scripts/generate_figures.py`, `tests/test_blockade.py`, `figures/fidelity_vs_blockade.png`.

**Existing tests**: 21 tests pass in 3.59s. None test blockade error or analytical formulas.

**Baseline numerical values verified at runtime**:
- `U/Ω ≈ 295.86` at `R = 3.0 µm`, `Ω/2π = 4 MHz`
- `ε_block = 1/(8 × 295.86²) ≈ 1.43 × 10⁻⁶` — very small at baseline

## 4. Root-cause candidates

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
|---|---|---|---|
| The blockade slice is simply unimplemented after the foundation work landed. | All required files absent. Issue explicitly requests implementation. | None — this is the direct state of the repo. | Create the files, make tests pass, generate figure and HTML section. | high |
| The blockade Hamiltonian builder might have a bug that makes the simulation silently wrong. | No direct evidence of a bug. `test_blockade_hamiltonian_matrix_elements_atom2` passes for atom 2, but atom 1 is not explicitly tested. | The builder follows the same pattern for both atoms. Matrix elements look correct on visual inspection. | Run a finite-blockade sim at very large U and verify it matches the ideal gate to high precision. This is already part of the acceptance criteria (`U/Ω = 1000 → infidelity < 10⁻⁸`). | low (bug unlikely) |
| The local-Z correction may not work correctly when the unitary is non-diagonal (finite blockade introduces small off-diagonal elements). | `apply_local_z_correction` reads phases from diagonal entries; if off-diagonal elements are large, the correction is not well-defined. | At moderate-to-strong blockade, the gate is nearly diagonal — off-diagonal leakage to `\|rr⟩` is small. The fidelity metric itself is invariant to local Z phases. | Test at U/Ω ≈ 10 where leakage is ~1% and verify the fidelity metric still gives a sensible answer. | medium |

## 5. Decision

- **Chosen likely cause**: Unimplemented feature slice. The Hamiltonian infrastructure is already in place.

- **Recommended fix**: Implement the blockade vertical slice end-to-end. This is the *simplest* error channel — purely coherent, `sesolve` only, no Monte Carlo, no Lindblad. The key design choices:

  **Simulation approach**: For each computational basis input, map to `BLOCKADE_BASIS`, run 3-segment `sesolve` with `build_blockade_hamiltonian`, extract the output amplitude in the computational subspace. The result is a 4×4 unitary (with small leakage to `|rr⟩` at finite U).

  **Mapping between computational and blockade bases**: The computational basis `{|00⟩, |01⟩, |10⟩, |11⟩}` maps to blockade-basis states as:
  - `|00⟩`: dark, no dynamics, amplitude = 1.
  - `|01⟩`: only atom 2 is in `|g⟩`. During the 2π pulse on atom 2, atom 1 is in `|0⟩` (dark), so atom 1 is never excited. This means no `|rr⟩` coupling exists — atom 2 evolves as a single atom. Use `single_atom_hamiltonian` or the 4-state Hamiltonian with atom 1 never coupling. Result: same as ideal.
  - `|10⟩`: only atom 1 is in `|g⟩`. Atom 2 is dark throughout. Atom 1 does π–π round trip as a single atom. No `|rr⟩` coupling. Result: same as ideal.
  - `|11⟩`: **both atoms in `|g⟩`**. This is the only input where `|rr⟩` matters. After the π pulse on atom 1 puts it into `|r⟩`, the 2π pulse on atom 2 can partially excite through `|rr⟩` because the blockade shift U is finite. **This is where the error lives.**

  **Critical insight**: Only the `|11⟩` input is affected by finite blockade. The other three inputs behave identically to the ideal gate. This means:
  1. The 4×4 output unitary is diagonal for `|00⟩`, `|01⟩`, `|10⟩` (same as ideal: `1, -1, -1`).
  2. The `|11⟩` output may have a small amplitude in `|rr⟩` (leakage) and a phase error on the `|gg⟩` component.
  3. The blockade error is entirely in the `|11⟩` column of the gate unitary.

  For fidelity, project the final state back to the computational subspace. Since this is a coherent error, the projected operator is still nearly unitary (up to the small `|rr⟩` leakage). Use `pedersen_fidelity` with `correct_local_z=True`.

  **Handling `|rr⟩` leakage**: At strong blockade, the `|rr⟩` population at the end of the gate is negligible. At weak blockade (small U/Ω), it can be significant. The projected 4×4 unitary will not be exactly unitary when there's leakage — its column norm for `|11⟩` will be < 1. This is physically correct: population has leaked outside the computational subspace.

  The Pedersen formula `F = [Tr(M†M) + |Tr(M)|²] / [d(d+1)]` still works on a sub-unitary matrix, and it naturally penalizes leakage because `Tr(M†M) < d`.

- **Rejected alternatives and why**:
  - *Use the 4-state basis for all inputs*: Unnecessary — `|00⟩`, `|01⟩`, `|10⟩` have no dynamics in the `|rr⟩` subspace. Simulating them in 4 states wastes solver time for the same answer.
  - *Add Lindblad terms*: Wrong — the issue explicitly says "coherent error, no Lindblad needed." Decay is a separate channel (issue #3).
  - *Optimal local-Z correction instead of `LOCAL_Z_PRODUCT`*: Potentially better, but `optimal_local_z_correction` exists and can be used. At strong blockade, the fixed `LOCAL_Z_PRODUCT` and the optimal correction give nearly identical results. At weak blockade, the optimal correction is more meaningful because phases shift. **Recommendation: use `optimal_local_z_correction` (or equivalently `correct_local_z=True`) for the fidelity, which is what `pedersen_fidelity` already does when that flag is set.**

- **Assumptions**:
  - The baseline operating point uses `get_rydberg_params()` values: `U/Ω ≈ 296` at `R = 3.0 µm`.
  - The AGENTS.md mentions `R ≈ 6 µm` as the Evered distance, but `src/params.py` defaults to `R = 3.0 µm` with the comment "3 um gives U/Omega ~= 300 for n=70; 6 um would not be strong blockade." The existing code convention is the source of truth for the baseline.
  - Units: rad/µs for frequencies and shifts, µs for time.

- **Maintainer decision needed, if any**: None. The blockade slice is straightforward with no ambiguous physics choices.

## 6. Implementation handoff

### Ordered implementation steps

1. **Create `src/errors/__init__.py`** (if not already created by issue #3). Empty file.

2. **Add `epsilon_blockade()` to `src/analytical.py`** (create file if needed):
   ```python
   def epsilon_blockade(omega: float, blockade_shift: float) -> float:
       """Perturbative blockade error: Ω²/(8U²)."""
       return omega**2 / (8.0 * blockade_shift**2)
   ```
   At baseline: `25.13² / (8 × 7435.74²) ≈ 1.43 × 10⁻⁶`.

3. **Create `src/errors/blockade.py`**:
   - `run_blockade_gate(omega, blockade_shift, n_steps_per_pi=400)` → dict with unitary, populations, phases
   - Only the `|11⟩` input uses the 4-state `BLOCKADE_BASIS` with `build_blockade_hamiltonian`. The other three inputs reuse the ideal single-atom / ideal-pair logic from `protocol.py`, or equivalently, just return known ideal phases.
   - The `|11⟩` simulation: start in `|gg⟩` of the 4-state basis, run 3 segments with `sesolve`, extract the final state, project back to computational subspace.
   - Return a 4×4 (projected) unitary matrix plus diagnostics (rr leakage, population traces if desired).

4. **Add blockade sweep to `src/sweeps.py`** (create file if needed):
   - `sweep_blockade(omega, blockade_shifts, ...)` → list of result dicts
   - Sweep over at least 20 points of U/Ω from ~5 to ~500 (log scale).
   - For each point, run `run_blockade_gate`, compute `pedersen_fidelity` with local-Z correction, also compute `epsilon_blockade`.

5. **Create `tests/test_blockade.py`**:
   - `test_blockade_large_u_is_ideal`: U/Ω = 1000 → infidelity < 10⁻⁸
   - `test_blockade_baseline_matches_analytical`: numerical vs analytical within 50%
   - `test_blockade_zero_blockade_raises_or_degrades`: U = 0 should either error or give low fidelity
   - `test_blockade_sweep_shape`: ≥20 points, correct range
   - `test_epsilon_blockade_formula`: known input → known output

6. **Update `scripts/run_sweeps.py`** (create if needed): add blockade sweep call, save data.

7. **Update `scripts/generate_figures.py`** (create if needed): generate `figures/fidelity_vs_blockade.png`:
   - x-axis: U/Ω (log scale)
   - y-axis: gate fidelity (or infidelity on log scale)
   - Numerical scatter + analytical curve `F = 1 - Ω²/(8U²)` + vertical line at baseline U/Ω ≈ 296
   - Publication quality matching the style in `scripts/run_ideal.py`

8. **Add "Finite Blockade" section to `site/index.html`**: Use the template at lines 982–1017.
   - Physics: the 2π pulse is detuned by U, not blocked perfectly; small leakage through `|rr⟩`
   - Formula: `ε_block = Ω²/(8U²)` in MathJax
   - Collapsible derivation sketch: second-order perturbation theory — the 2π pulse is a 2-level rotation detuned by U; the transition probability is `(Ω/2)²/((Ω/2)² + (U/2)²)` ≈ `Ω²/U²` for U >> Ω; the factor of 1/8 comes from the specific gate protocol
   - Figure embed
   - Takeaway: "This error dominates when interaction is weak relative to Rabi frequency"

9. **Force-commit the figure**: `git add -f figures/fidelity_vs_blockade.png` (per issue #3 decision).

10. **Full validation**: run targeted + full test suites.

### Files likely to change

| File | Action |
|---|---|
| `src/errors/__init__.py` | New (empty) |
| `src/errors/blockade.py` | New |
| `src/analytical.py` | New or extend (add `epsilon_blockade`) |
| `src/sweeps.py` | New or extend (add `sweep_blockade`) |
| `scripts/run_sweeps.py` | New or extend |
| `scripts/generate_figures.py` | New or extend |
| `tests/test_blockade.py` | New |
| `site/index.html` | Modify (add Finite Blockade section) |
| `figures/fidelity_vs_blockade.png` | Generated |

### Scope boundaries / what not to touch

- Do NOT add Lindblad/dissipative dynamics. This is a purely coherent error.
- Do NOT implement decay, Doppler, scattering, amplitude noise, or combined budget.
- Do NOT refactor `src/protocol.py` unless extracting a small shared helper is clearly cleaner than duplication.
- Do NOT change `src/hamiltonian.py` — `build_blockade_hamiltonian` is already correct.
- Do NOT change `src/params.py` defaults.

### Regression test shape

```python
def test_blockade_large_u_is_ideal():
    """U/Omega = 1000 should give nearly perfect CZ."""
    omega = 25.13
    U = omega * 1000
    result = run_blockade_gate(omega, U)
    F = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
    assert 1.0 - F < 1e-8

def test_blockade_baseline_matches_analytical():
    params = get_rydberg_params()
    result = run_blockade_gate(params.omega_rad_per_us, params.blockade_shift_rad_per_us)
    F = pedersen_fidelity(result["unitary"], CZ_TARGET, correct_local_z=True)
    numerical_error = 1.0 - F
    analytical_error = epsilon_blockade(params.omega_rad_per_us, params.blockade_shift_rad_per_us)
    # At U/Omega ~ 296, analytical gives 1.43e-6
    assert abs(numerical_error - analytical_error) / analytical_error < 0.5

def test_epsilon_blockade_formula():
    assert np.isclose(epsilon_blockade(1.0, 10.0), 1.0 / 800.0)
```

### Validation commands

```bash
# Targeted blockade tests
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q tests/test_blockade.py

# Generate sweep data + figure
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/run_sweeps.py
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/generate_figures.py

# Full regression
PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q
```

### Expected pre-fix signal

- `import src.errors.blockade` → `ModuleNotFoundError`
- `pytest tests/test_blockade.py` → file not found
- `figures/fidelity_vs_blockade.png` → absent
- `site/index.html` → template placeholder only

### Expected post-fix signal

- All new and existing tests pass (21 existing + new blockade tests)
- `U/Ω = 1000` infidelity < 10⁻⁸
- Baseline infidelity ≈ 1.43 × 10⁻⁶ (within 50% of analytical)
- Sweep: ≥20 points, U/Ω from ~5 to ~500
- Figure exists with dots, curve, baseline line
- HTML contains "Finite Blockade" section; existing scaffold tests pass

### Negative/edge validation

- `blockade_shift = 0` → should raise `ValueError` (division by zero in analytical formula, physically meaningless)
- `blockade_shift < 0` → raise `ValueError` (non-physical)
- `U/Ω < 1` → fidelity collapses, perturbative formula breaks down, acceptable if test only checks that fidelity is low
- Very small U/Ω → significant `|rr⟩` leakage at gate end; projected unitary is substantially sub-unitary

## 7. Risks, edge cases, rollback

### Risks

| Risk | Mitigation |
|---|---|
| **Fidelity metric on sub-unitary matrix**: At low U/Ω, `|rr⟩` leakage makes the projected 4×4 matrix non-unitary. | Pedersen formula handles this correctly — `Tr(M†M) < d` naturally penalizes leakage. Verify by checking that infidelity > 0 at low U/Ω. |
| **Local-Z correction on a non-diagonal gate**: With finite blockade, the gate has small off-diagonal elements (population in `|rr⟩` that leaks back). | Use `correct_local_z=True` which calls `optimal_local_z_correction`. If the diagonal entries are near-zero (extreme leakage), it will raise `ValueError` — add a guard for very small U/Ω. |
| **Phase convention mismatch with ideal gate**: The ideal gate gives `diag(1, -1, -1, -1)`. The blockade gate should approach this at large U. | Test explicitly that blockade at U/Ω = 1000 gives phases matching the ideal gate within 10⁻⁸. |
| **Shared file conflicts with issue #3**: Both issues create `src/errors/__init__.py`, `src/analytical.py`, `src/sweeps.py`. | If #3 lands first, extend those files. If not, create them. The files are simple enough that merge conflicts are trivial. |

### Edge cases

- **U/Ω → ∞**: Must reproduce the ideal gate exactly (this is the physical sanity check).
- **U/Ω ~ 1**: Gate completely fails. Perturbative formula is invalid. Numerical simulation should still run and give a very low fidelity.
- **U/Ω ~ 5**: Near the boundary of the sweep range. Perturbative formula starts to deviate from numerics — this is physically interesting and the figure should show it.

### Rollback/revert notes

This slice is purely additive. Revert = remove new files + revert the `site/index.html` modification. No impact on ideal gate or other error channels.

## 8. Open questions

### 1. Should the sweep x-axis be U/Ω or Ω/U?

- **What was found**: The issue says "x-axis: U/Ω (log scale)" and "Vary ratio U/Ω". The analytical formula `ε = Ω²/(8U²) = 1/(8(U/Ω)²)` decreases with increasing U/Ω. AGENTS.md says "Sweep parameter: Vary the ratio U/Ω."
- **Why it matters**: Convention choice for plotting. With U/Ω on x-axis, the fidelity curve goes *up* to the right (better blockade → better fidelity). This is intuitive.
- **Recommended default**: Use U/Ω on x-axis as stated in the issue. Fidelity increases left-to-right, which is the natural reading direction for "more blockade = better."

### 2. Should the blockade simulation reuse `protocol.py` helpers or be self-contained?

- **What was found**: `protocol.py` has private helpers (`_solve_segment`, `_stitch_segments`, `_run_input_11`, etc.) that implement the 3-segment pulse pattern. The blockade module needs the same pattern but with 4-state Hamiltonians for the `|11⟩` input.
- **Why it matters**: Code reuse vs. coupling. Making `protocol.py` helpers public creates a contract. Duplicating them in `blockade.py` means parallel maintenance.
- **Options**:
  1. Import and reuse `_solve_segment` from `protocol.py` (rename to public).
  2. Duplicate the ~10 lines of solver logic in `blockade.py`.
  3. Extract shared helpers to a new `src/solver_utils.py`.
- **Recommended default**: Option 2. The solver logic is ~10 lines. Duplicating it keeps the blockade module self-contained and avoids changing `protocol.py`'s interface. Refactor when a third module needs the same pattern.
