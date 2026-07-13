# Issue #8: Combined error budget + summary table + final polish

## 1. Executive summary

- **Actionability**: actionable, but blocked by #7
- **Symptom**: No combined error simulation, no summary table, no error-budget bar chart. The HTML explainer has 4 of 5 individual error sections but lacks the amplitude noise channel, the combined budget, the summary table, and final polish.
- **User impact / severity**: High — this is the capstone issue for the milestone. Without it, the deliverable is an incomplete analysis. With it, a reader can compare all error sources at the Evered operating point and see whether perturbation theory holds.
- **Recommended next step**: Wait for issue #7 (amplitude noise) to close, then implement the combined simulation + HTML sections. The combined module is the most technically complex piece in the project (9-state Hilbert space, Lindblad + MC hybrid, process-matrix extraction).
- **Confidence**: medium (physics plan is clear; implementation complexity is real)
- **Main blocker**: **Issue #7 (amplitude noise) is still OPEN.** The combined simulation and summary table require all 5 individual error channels. Issues #3, #4, #5, #6 are closed. #7 is the only remaining blocker.

## 2. Affected-path context

### Affected subsystem

The integration/capstone slice: `src/errors/combined.py`, summary scripts, bar chart figure, three new HTML sections (amplitude noise from #7, combined budget, summary table), and final HTML polish.

### Repo vocabulary / components needed to understand this issue

| Term | Meaning |
|---|---|
| DECAY_BASIS | 8-state basis `(00, 0g, 0r, g0, gg, gr, r0, rg)` — projects out `\|rr⟩`. Used by decay + scattering modules. |
| BLOCKADE_BASIS | 4-state basis `(gg, gr, rg, rr)` — includes `\|rr⟩`. Used by blockade module. |
| COMBINED_BASIS | **New**: 9-state basis `(00, 0g, 0r, g0, gg, gr, r0, rg, rr)` — needed for this issue. |
| `mesolve` | QuTiP master-equation solver. Used when Lindblad channels (decay, scattering) are active. |
| MC outer loop | Monte Carlo over Doppler velocities + amplitude noise ε. Each sample runs a full `mesolve` gate. |
| Choi matrix | 16×16 matrix encoding the full CPTP map. Built from 16 `mesolve` propagations of basis operators. |
| Pedersen fidelity (Kraus) | `pedersen_fidelity_kraus(kraus_ops, CZ_TARGET)` — the open-system fidelity metric. |
| Additive approximation | `ε_total ≈ Σ ε_i` — valid when errors are small and independent. |

### Current relevant control/data flow

The combined simulation sits on top of ALL existing modules:

```
src/params.py                    → baseline parameters
src/errors/decay.py              → Lindblad channel (γ decay)
src/errors/scattering.py         → Lindblad channel (dephasing)
src/errors/blockade.py           → finite U Hamiltonian
src/errors/doppler.py            → MC sampling of velocities
src/errors/amplitude.py [#7]     → MC sampling of Ω noise
─────────────────────────────────────────────────────
src/errors/combined.py [NEW]     → unified 9-state mesolve + MC
src/analytical.py                → all ε_i formulas for additive sum
src/sweeps.py                    → point evaluation at Evered baseline
scripts/run_sweeps.py            → combined run
scripts/generate_figures.py      → bar chart
site/index.html                  → combined + summary + polish
```

### Why the key files matter

- **`src/errors/decay.py`**: The process-matrix approach (`_propagate_process_outputs`, `_choi_from_process_outputs`, `_kraus_from_choi`) is the template for the combined module. The combined module reuses these algorithms in a larger Hilbert space.
- **`src/errors/scattering.py`**: Imports and reuses decay.py's internal infrastructure (`_hamiltonian`, `_collapse_operators`, `_embed_computational_operator`, etc.). Shows how to build Lindblad channels on top of the decay basis.
- **`src/errors/doppler.py`**: Shows the MC sampling pattern: `rng.normal(...)`, loop over samples, average `pedersen_fidelity`.
- **`src/hamiltonian.py`**: Has `build_blockade_hamiltonian(omega, blockade_shift, atom_index, detuning)` returning 4×4 in BLOCKADE_BASIS with Doppler detuning support.
- **`site/index.html`**: 1573 lines. Sections for intro, ideal gate, decay, blockade, doppler, scattering. Template for remaining sections at line 1226. Footer at line 1258. Navigation bar at line 680.
- **`src/fidelity.py`**: `pedersen_fidelity_kraus` for open-system channels. `LOCAL_Z_PRODUCT` for phase correction.

### First file a new implementer should open

`src/errors/decay.py` → the process-matrix / Choi / Kraus pipeline. Then `src/errors/doppler.py` → the MC sampling pattern. Then `src/hamiltonian.py` → `build_blockade_hamiltonian` for the finite-U Hamiltonian.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Integration slice. Run all 5 errors simultaneously, compare against additive approximation, produce bar chart and summary table, polish HTML.
- **Reproduction steps**: N/A — feature request.
- **Expected behavior**:
  - Combined simulation at Evered 2023 point with all 5 errors
  - |ε_sum − ε_numerical| / ε_numerical < 20%
  - Bar chart: `figures/error_budget_combined.png`
  - HTML: combined budget section, summary table, all 8 sections present
  - Full `pytest tests/` < 60 seconds
  - Mobile responsive, no raw LaTeX visible
- **Actual behavior**: All required files absent. Issue #7 still open.
- **Blocked by**: #3 (closed), #4 (closed), #5 (closed), #6 (closed), **#7 (OPEN)**
- **Environment**: `qutip==5.0.4`, `ARC==3.9.0`
- **Relevant comments**: None.

### Repository evidence

| File / fact | What it proves |
|---|---|
| `src/errors/decay.py` (250 lines) | Process-matrix pipeline exists: `_propagate_process_outputs` (16 mesolve runs), `_choi_from_process_outputs`, `_kraus_from_choi`. Operates in DECAY_BASIS (8 states, no `\|rr⟩`). |
| `src/errors/scattering.py` (185 lines) | Imports and extends the decay module's internal functions. Shows that Lindblad channels in the 8-state basis are composable (both decay and scattering use `_hamiltonian`, `_collapse_operators` etc.) |
| `src/errors/blockade.py` (200 lines) | Operates in BLOCKADE_BASIS (4 states, includes `\|rr⟩`). Uses `sesolve` (coherent only). Cannot be directly composed with decay/scattering modules due to basis mismatch. |
| `src/errors/doppler.py` (200 lines) | MC loop: samples velocities, calls `projected_gate_for_detunings(omega, delta1, delta2)`, averages `pedersen_fidelity`. Operates in 2- and 3-state ideal basis (no `\|rr⟩`, no Lindblad). |
| `src/errors/amplitude.py` | **ABSENT** — issue #7 not yet implemented |
| `src/hamiltonian.py` | `build_blockade_hamiltonian(omega, blockade_shift, atom_index, detuning)` supports both finite-U and Doppler detuning in 4-state basis. Key building block for the combined Hamiltonian, but needs extension to 9-state. |
| `src/analytical.py` | Has `epsilon_decay_from_gamma`, `epsilon_blockade`, `epsilon_doppler`, `epsilon_scattering`. Missing `epsilon_amplitude` (issue #7). |
| `src/sweeps.py` | Has 4 sweep functions + 4 CSV read/write pairs. No combined or amplitude support. |
| `site/index.html` (1573 lines) | Sections: intro, ideal, decay, blockade, doppler, scattering. Missing: amplitude, combined, summary. Navigation bar has 6 links. Footer says "remaining error sources, combined budget, and summary table are intentionally left to later issue slices." |
| `tests/` (575 lines total) | 8 test files. No `test_combined.py`. `test_html_scaffold.py` asserts specific section IDs exist. |
| `.gitignore` | `figures/*.png` and `figures/*.csv` gitignored; force-commit for publication figures. |

**Baseline values at the Evered operating point** (from existing code and analytical formulas):

| Error source | Analytical formula | Baseline ε (approx) | Status |
|---|---|---|---|
| Rydberg decay | γ × 7π/(4Ω) | 5.9 × 10⁻⁴ | ✅ implemented |
| Finite blockade | Ω²/(8U²) | 1.4 × 10⁻⁶ | ✅ implemented |
| Doppler dephasing | (π²/4) k² kT/(mΩ²) | 9.5 × 10⁻⁵ | ✅ implemented |
| Intermediate scattering | (7π/8)(Γ_e/Δ_p) | 1.7 × 10⁻² | ✅ implemented |
| Amplitude noise | (π²/2)σ² | ~2.0 × 10⁻³ | ❌ #7 open |
| **Additive sum** | | **~1.9 × 10⁻²** | |

Note: intermediate-state scattering at the default Δ_p/2π = 1 GHz dominates the error budget by an order of magnitude.

## 4. Root-cause candidates

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
|---|---|---|---|---|
| Unimplemented capstone slice — all prerequisite channels exist except amplitude (#7). | Required files absent. 4 of 5 blockers resolved. | Issue #7 is still open, making the combined simulation impossible to implement correctly now. | Wait for #7, then implement. | high |
| The 9-state combined basis may introduce new numerical issues that the isolated 8-state and 4-state modules don't have. | No evidence yet — hasn't been built. The 9×9 density matrix in `mesolve` is small enough for any solver. | QuTiP handles 9-state systems easily. The concern is the 9² = 81-element density matrix × 16 process inputs × N MC samples = moderate compute. | Benchmark: single combined-gate evaluation at 9-state should take < 30 seconds with N_MC=100. | low (unlikely issue) |
| The additive approximation may not hold at 20% because scattering (~1.7%) is not "small." | Perturbation theory assumes ε ≪ 1. At ε_scatter ~ 2%, cross-terms with other channels could be O(10⁻⁴), which is comparable to ε_decay. | The cross-term ε_scatter × ε_decay ~ 10⁻⁵ is much smaller than either individual error. The dominant non-additivity would be ε_scatter × ε_amp ~ 3 × 10⁻⁵, still small relative to ε_scatter. | Run the combined simulation and check. If |ε_sum − ε_num|/ε_num < 20%, the additive approximation is validated. | medium (likely holds) |

## 5. Decision

- **Chosen likely cause**: Unimplemented capstone slice. The implementation is blocked by issue #7 (amplitude noise).

- **Recommended fix**: Two-phase implementation. Phase 1: wait for #7 to close. Phase 2: implement the combined module + HTML sections + final polish.

### Phase 2 design: the combined simulation

**The central technical challenge**: Each existing error channel uses a different Hilbert space and solver. The combined simulation must unify them.

| Channel | Hilbert space | Solver | Mechanism |
|---|---|---|---|
| Decay | 8-state DECAY_BASIS | `mesolve` | Lindblad √γ \|g⟩⟨r\| |
| Scattering | 8-state DECAY_BASIS | `mesolve` | Lindblad √γ_eff \|r⟩⟨r\| |
| Blockade | 4-state BLOCKADE_BASIS | `sesolve` | Hamiltonian: U \|rr⟩⟨rr\| |
| Doppler | 2/3-state ideal basis | `sesolve` | Hamiltonian: −δ \|r⟩⟨r\| |
| Amplitude | 2/3-state ideal basis | `sesolve` | Hamiltonian: Ω(1+ε) coupling |

**Solution**: A new **9-state COMBINED_BASIS**:

```
("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg", "rr")
```

This is the DECAY_BASIS + the |rr⟩ state. It accommodates:
- Dark qubit |0⟩ (present in DECAY_BASIS)
- Rydberg decay and scattering (Lindblad channels from DECAY_BASIS)
- Finite blockade (U|rr⟩⟨rr| from BLOCKADE_BASIS)
- Doppler detuning (−δ|r⟩⟨r| per atom)
- Amplitude noise (modified Ω in coupling terms)

**Hamiltonian in COMBINED_BASIS**: For a pulse on atom i with Rabi frequency Ω'=Ω(1+ε) and Doppler detuning δ_i:

```
H = (Ω'/2)(|g⟩⟨r|_i + h.c.) − δ_i |r⟩⟨r|_i + U |rr⟩⟨rr|
```

This is `build_blockade_hamiltonian` extended to 9 states with the dark qubit |0⟩ and detuning support.

**Collapse operators in COMBINED_BASIS**:
- Decay: √γ |g⟩⟨r| for each atom (4 terms: 0r→0g, gr→gg, r0→g0, rg→gg, and the rr→gr/rg terms)
- Scattering dephasing: √γ_eff |r⟩⟨r| on the driven atom only (applied per-segment)

**MC outer loop**:
```python
for sample in range(N):
    eps = rng.normal(0, sigma_omega)     # amplitude noise
    v1, v2 = rng.normal(0, v_rms, 2)    # Doppler velocities
    omega_eff = omega * (1 + eps)
    delta1 = k_eff * v1
    delta2 = k_eff * v2
    
    # Build 9-state Hamiltonian + collapse operators
    # Run 3-segment mesolve
    # Propagate 16 process basis operators → Choi → Kraus → fidelity
    fidelities[sample] = pedersen_fidelity_kraus(kraus_ops, CZ_TARGET)
```

**Performance estimate**: 16 mesolve runs per sample × 3 segments = 48 mesolve calls. With 9×9 density matrices, each takes ~5–10 ms. Per sample: ~0.3 s. At N=100: ~30 s. At N=500 for final figure: ~2.5 min. Acceptable for scripts, must be capped at N≤50 for tests.

### The bar chart and summary table

**Bar chart** (`figures/error_budget_combined.png`):
- 7 bars: ε_decay, ε_blockade, ε_doppler, ε_scatter, ε_amp, ε_sum (additive), ε_numerical
- The first 5 bars use individual-channel numerical values at the Evered point
- ε_sum = Σ ε_i
- ε_numerical = combined simulation result
- Log scale on y-axis (errors span 10⁻⁶ to 10⁻²)
- Color-coded, with ε_sum and ε_numerical in distinct style (hatched/outlined)

**Summary table** (HTML):
| Error source | Analytical ε | Numerical ε | Dominant scaling |
|---|---|---|---|
| Rydberg decay | γ·7π/(4Ω) | (from sweep) | 1/(Ωτ) |
| Finite blockade | Ω²/(8U²) | (from sweep) | (Ω/U)² |
| Doppler | (π²/4)k²kT/(mΩ²) | (from MC) | T/Ω² |
| Scattering | (7π/8)(Γ_e/Δ_p) | (from Lindblad) | 1/Δ_p |
| Amplitude noise | (π²/2)σ² | (from MC) | σ² |
| **Total (additive)** | **Σ ε_i** | | |
| **Total (numerical)** | | **(from combined)** | |

### HTML final polish checklist

1. Add amplitude noise section (prerequisite: #7)
2. Add combined error budget section with bar chart
3. Add summary table section
4. Update navigation bar: add "Amplitude", "Combined", "Summary" links
5. Update footer text
6. Verify all figures load (onerror fallbacks)
7. Verify MathJax renders everywhere (no raw `\(`, `$$`)
8. Verify collapsible `<details>` sections work
9. Verify mobile layout (no horizontal scroll at 375px width)
10. Update `PAGE_TITLE` in JS

- **Rejected alternatives and why**:
  - *Run each error channel independently and just add the numbers*: Misses the entire point. The issue asks whether the additive approximation holds — that requires the combined simulation.
  - *Use the DECAY_BASIS (8-state) for the combined simulation, projecting out |rr⟩*: Wrong — this would omit the finite-blockade error entirely.
  - *Use the BLOCKADE_BASIS (4-state) for the combined simulation*: Missing the dark qubit |0⟩, which is needed to properly handle |00⟩, |01⟩, |10⟩ inputs where one qubit is in the non-coupled state.
  - *Average per-sample unitaries instead of per-sample fidelities*: Wrong for the combined channel because each sample includes dissipative (Lindblad) dynamics, so the per-sample evolution is non-unitary.

- **Assumptions**:
  - The Evered operating point uses `get_rydberg_params()` defaults + σ_Ω = 0.02 + Δ_p/2π = 1 GHz + T = 10 μK.
  - Issue #7 will provide `epsilon_amplitude` in `src/analytical.py` and `run_amplitude_noise_gate` in `src/errors/amplitude.py`.
  - The combined module only evaluates at one operating point (the Evered baseline), not a sweep. A sweep would be impractical given the MC + Lindblad computational cost.

- **Maintainer decision needed**: See Open Questions §8.

## 6. Implementation handoff

### Ordered implementation steps

**Prerequisites**: Issue #7 must be closed first.

1. **Create `src/errors/combined.py`**:
   - Define `COMBINED_BASIS = ("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg", "rr")` — the 9-state space.
   - Build a 9×9 Hamiltonian for each pulse segment: Rabi coupling + Doppler detuning + blockade shift.
   - Build 9×9 collapse operators: decay (√γ |g⟩⟨r| for both atoms, including rr→rg/gr channels) + scattering dephasing (√γ_eff |r⟩⟨r| on driven atom).
   - MC outer loop: sample (ε, v₁, v₂), build Hamiltonian/collapse ops, run 3-segment `mesolve`.
   - Process-matrix extraction: propagate 16 basis operators (following decay.py's `_propagate_process_outputs` pattern), build Choi matrix, extract Kraus operators, compute Pedersen fidelity.
   - Result dataclass: `CombinedGateResult(average_fidelity, std_fidelity, individual_errors: dict, additive_sum, n_samples, ...)`.
   - Main entry: `run_combined_gate(omega, gamma, blockade_shift, gamma_e, delta_p, sigma_omega, k_eff, temperature, mass, *, n_samples=100, seed=None)`.

2. **Collect individual errors at the operating point**:
   - Call each existing `run_*_gate(...)` function with Evered baseline parameters.
   - Store the individual numerical ε_i values.
   - Compute ε_sum = Σ ε_i.
   - These are the "per-channel" bars in the bar chart.

3. **Add `epsilon_total_additive` to `src/analytical.py`**:
   ```python
   def epsilon_total_additive(
       omega, gamma, blockade_shift, gamma_e, delta_p, 
       sigma_omega, k_eff, temperature, mass
   ) -> tuple[float, dict[str, float]]:
       """Return (sum, {channel: value}) for the additive budget."""
   ```

4. **Create `tests/test_combined.py`**:
   - `test_combined_additive_approximation`: |ε_sum − ε_num|/ε_num < 20%
   - `test_combined_all_zero_is_ideal`: zero all error params → F = 1
   - `test_combined_single_error_matches_isolated`: turn on only decay → matches `run_decay_gate` within 10%
   - Use N ≤ 50 for MC, fixed seed, keep test < 30 seconds.

5. **Update `scripts/run_sweeps.py`**: Add combined evaluation at the Evered point. Print the additive vs numerical comparison.

6. **Update `scripts/generate_figures.py`**: Add `plot_error_budget_bar_chart()`:
   - 7 bars on log scale
   - Colors matching existing figure palette
   - Clear labels and legend
   - Publication quality

7. **Add three HTML sections to `site/index.html`**:
   
   a. **"Amplitude Noise" section** (if not already added by #7): `id="sec-error-amplitude"`
   
   b. **"Combined Error Budget" section**: `id="sec-combined-budget"`
     - Explain the additive approximation and why it might fail
     - Embedded bar chart figure
     - Discussion: does ε_sum ≈ ε_numerical? What are the cross-terms?
   
   c. **"Summary Table" section**: `id="sec-summary-table"`
     - HTML `<table>` with all 5 error sources + totals
     - Columns: source, analytical ε, numerical ε, dominant scaling
     - Highlight which error dominates

8. **Final HTML polish**:
   - Update navigation bar (add Amplitude, Combined, Summary links)
   - Update footer text
   - Verify MathJax rendering end-to-end
   - Verify mobile layout
   - Update `PAGE_TITLE` in JS
   - Check all `onerror` fallbacks work for missing figures

9. **Validate test timing**: `pytest tests/` must complete in < 60 seconds. If MC tests are slow, reduce N.

10. **Force-commit**: `git add -f figures/error_budget_combined.png`

### Files likely to change

| File | Action |
|---|---|
| `src/errors/combined.py` | New (~250–350 lines) |
| `src/analytical.py` | Extend (add `epsilon_total_additive`) |
| `src/sweeps.py` | Extend (add combined point evaluation, optional) |
| `scripts/run_sweeps.py` | Extend (add combined run) |
| `scripts/generate_figures.py` | Extend (add bar chart) |
| `tests/test_combined.py` | New |
| `site/index.html` | Major modification (3 new sections, nav bar, footer, polish) |
| `figures/error_budget_combined.png` | Generated |

### Scope boundaries / what not to touch

- Do NOT modify existing error channel modules (#3–#7 slices). The combined module calls them, it doesn't alter them.
- Do NOT add new error physics. The combined module only turns on what already exists.
- Do NOT modify `src/params.py` unless the amplitude noise default σ_Ω is not yet there from #7.
- Do NOT change existing tests. The combined test is additive; existing tests must still pass.
- Do NOT implement parameter sweeps for the combined simulation (too expensive). A single-point evaluation at the Evered baseline is sufficient.

### Regression test shape

```python
def test_combined_all_zero_errors_is_ideal():
    """All error parameters zero → F = 1."""
    result = run_combined_gate(
        omega=DEFAULT_OMEGA, gamma=0.0, blockade_shift=1e12,
        gamma_e=0.0, delta_p=1e6, sigma_omega=0.0,
        k_eff=0.0, temperature=0.0, mass=RB87_MASS_KG,
        n_samples=2, seed=42,
    )
    assert 1.0 - result.average_fidelity < 1e-8

def test_combined_additive_approximation():
    """Additive sum matches full numerical within 20%."""
    result = run_combined_gate(...)  # Evered baseline
    additive = result.additive_sum
    numerical = 1.0 - result.average_fidelity
    assert abs(additive - numerical) / numerical < 0.20
```

### Validation commands

```bash
# Run combined simulation at Evered point
PYTHONPATH=. python scripts/run_sweeps.py

# Generate all figures including bar chart
PYTHONPATH=. python scripts/generate_figures.py

# Combined tests only
PYTHONPATH=. python -m pytest -xvs tests/test_combined.py

# Full regression with timing
PYTHONPATH=. python -m pytest -q --tb=short -q 2>&1 | tail -5
# Must see: "X passed in <60.00s"
```

Note: exact Python executable depends on environment. Tests use `PYTHONPATH=.` convention.

### Expected pre-fix signal

- `import src.errors.combined` → `ModuleNotFoundError`
- `import src.errors.amplitude` → `ModuleNotFoundError` (until #7 lands)
- `figures/error_budget_combined.png` → absent
- `site/index.html` → no `sec-combined-budget` or `sec-summary-table` sections
- Footer says "remaining error sources, combined budget, and summary table are intentionally left"

### Expected post-fix signal

- All tests pass in < 60 seconds
- Combined at Evered point: |ε_sum − ε_num|/ε_num < 20%
- All-zero-errors → F = 1
- Bar chart exists with 7 bars
- HTML has all 8 content sections + summary table
- Navigation bar has all links
- MathJax renders (no raw `\(` visible)
- Mobile responsive (no horizontal scroll at 375px)

### Negative/edge validation

- All individual errors zero → combined = ideal gate (F = 1)
- Only one error active → combined ≈ isolated channel (within 10%)
- All error sources at 10× baseline → fidelity drops significantly; additive approximation may visibly break
- `n_samples < 2` → should work (single sample, meaningless statistics)

## 7. Risks, edge cases, rollback

### Risks

| Risk | Mitigation |
|---|---|
| **Issue #7 not resolved**: Cannot implement combined simulation or complete HTML without amplitude noise channel. | This is a hard blocker. The triage explicitly marks it. Do not begin implementation until #7 is closed. |
| **9-state Hilbert space introduces new Lindblad terms**: The |rr⟩ state needs decay channels (rr→gr and rr→rg). These don't exist in the current 8-state decay module. | Build them fresh in `combined.py`. The physics is: if either atom in |rr⟩ decays, the state goes to |gr⟩ or |rg⟩ with rate γ per atom. |
| **Combined mesolve is slow**: 16 process inputs × 3 segments × N_MC samples = O(50N) mesolve calls. At N=500 for figure-quality data: ~5 minutes. | Cap tests at N≤50 (< 30s). Script data at N=200–500. Cache results in CSV. |
| **Additive approximation test flakiness**: MC variance with N=50 could push the ratio beyond 20%. | Use fixed seed. At N=50, the SE on ε_total is ~σ/√50 ≈ 14% of ε. Combined with the dominant ε_scatter being deterministic (Lindblad, not MC), the MC variance affects only the Doppler + amplitude contributions, which are subdominant. Total MC variance on ε_total is < 5%. |
| **HTML conflicts with concurrent issues**: If issues #5–#7 modified site/index.html, merge conflicts are likely. | This issue is designed to run last. By the time it runs, all error-section HTML should be in place. The combined/summary sections are inserted fresh. |
| **Test suite timing > 60s**: Adding combined tests (with MC) may push the suite over. | Use N ≤ 50 in combined tests. If existing MC tests are slow, reduce their N too. Profile with `pytest --durations=10`. |

### Edge cases

- **Scattering dominates**: At the default Δ_p/2π = 1 GHz, ε_scatter ≈ 1.7% is 10× larger than any other error. The bar chart will be dominated by one bar. This is physically correct but may look odd visually — use log scale.
- **Cross-channel interference**: The main candidate is scattering (dephasing) + decay (population transfer). Both operate on the same Rydberg coherence. The cross-term is O(γ × γ_eff) ≈ O(10⁻⁵), negligible at the Evered point.
- **Finite blockade + decay**: Population in |rr⟩ can decay to |gr⟩ or |rg⟩. This channel doesn't exist in the isolated modules (blockade ignores decay; decay projects out |rr⟩). The combined module naturally captures it, but the effect is tiny (|rr⟩ population ~ 10⁻⁶ × decay during gate ~ 10⁻⁸).

### Rollback/revert notes

This slice is additive. Revert = remove `src/errors/combined.py`, `tests/test_combined.py`, revert script/figure/HTML modifications. No impact on isolated error channels.

## 8. Open questions

### 1. Should the combined simulation use the full 16-input process matrix or a cheaper approximation?

- **What was found**: The decay/scattering modules use the full 16-input approach (`_propagate_process_outputs`). This is the gold standard but costs 16 × 3 = 48 mesolve calls per MC sample.
- **Why it matters**: At N=500, that's 24,000 mesolve calls. At ~5ms each (9×9 matrices), ~2 minutes. Acceptable for scripts but slow for interactive development.
- **Options**:
  1. Full 16-input process matrix per MC sample, averaged over samples. (Correct, ~2 min at N=500)
  2. Per-sample: propagate only the 4 diagonal inputs (|i⟩⟨i|), compute a "diagonal channel fidelity." Average over samples. (Faster: 4×3=12 calls/sample. Approximate.)
  3. Compute the full 16-input Choi for a single "representative" set of (ε, v₁, v₂) = (0, 0, 0), then add MC corrections perturbatively. (Fastest. Most approximate.)
- **Recommended default**: Option 1 for correctness. The combined simulation runs once at the Evered point, not as a sweep. Two minutes is fine for a capstone evaluation.

### 2. How should the `rr` decay channels be handled?

- **What was found**: The existing 8-state DECAY_BASIS doesn't include |rr⟩. The combined 9-state basis adds it. When |rr⟩ is populated (finite blockade), either atom can decay: |rr⟩ → |gr⟩ or |rr⟩ → |rg⟩.
- **Why it matters**: The collapse operators √γ |g⟩⟨r|₁ and √γ |g⟩⟨r|₂ in the 9-state basis naturally include the rr→gr and rr→rg terms. This is correct and requires no special handling — just construct the collapse operators in the 9-state basis.
- **Recommended default**: Build the 9-state collapse operators from scratch in `combined.py`. The existing 8-state operators in decay.py cannot be simply extended.

### 3. Should the combined test assert a specific numerical value for ε_total or just the 20% agreement?

- **What was found**: The individual channels have known baseline values. The additive sum is ~1.9 × 10⁻². The combined simulation should give a similar value.
- **Why it matters**: Asserting a specific value makes the test more informative but also more fragile (depends on ARC version, numerical tolerances, MC seed).
- **Options**:
  1. Assert only relative agreement: |ε_sum − ε_num|/ε_num < 20%.
  2. Additionally assert ε_num is in a plausible range (e.g., 0.5% to 5%).
  3. Assert both.
- **Recommended default**: Option 3. The relative test is the physics check. The absolute range test catches gross bugs (e.g., F = 0 or F = 1 when it shouldn't be).
