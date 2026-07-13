# Issue #7: Amplitude noise: simulation, analytical formula, sweep, figure, HTML section

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: The repository has no amplitude noise error channel. No `src/errors/amplitude.py`, no `epsilon_amplitude()` analytical formula, no amplitude sweep in `src/sweeps.py`, no `figures/fidelity_vs_amplitude_noise.png`, and no "Amplitude Noise" HTML section. This is a feature-implementation issue.
- **User impact / severity**: High for the milestone. This is the fifth (final) individual error channel before the combined budget (issue #8). It completes the single-channel analysis.
- **Recommended next step**: Implement the amplitude noise vertical slice end-to-end. The existing `run_ideal_gate(omega)` in `src/protocol.py` already runs the coherent π–2π–π gate at any Rabi frequency, so the core MC loop is: sample ε, call the gate at Ω(1+ε), compute fidelity, average. The analytical formula is `ε_amp = (π²/2) σ²` — derived below.
- **Confidence**: high
- **Main blocker, if any**: Issues #1 and #2 (blockers listed in issue) are both closed and merged. The existing infrastructure is sufficient. Issues #5 (Doppler) and #6 (scattering) are open but independent — these error channels are designed for parallel development.

## 2. Affected-path context

### Affected subsystem

The fifth error-channel vertical slice: `src/errors/amplitude.py`, `src/analytical.py`, `src/sweeps.py`, sweep/figure scripts, amplitude tests, and the "Amplitude Noise" HTML section.

### Repo vocabulary / components needed to understand this issue

| Term | Meaning |
|---|---|
| `\|0⟩` | Dark qubit state. Does not couple to the Rydberg laser. |
| `\|g⟩` | Qubit `\|1⟩`. Ground state that couples to `\|r⟩`. |
| `\|r⟩` | Auxiliary Rydberg state (`\|70S₁/₂⟩`). |
| `Ω` | Two-photon Rabi frequency. Solver units: rad/μs (`omega_rad_per_us`). |
| `σ_Ω` | Fractional Rabi frequency noise. Ω → Ω(1+ε), ε ~ N(0, σ²). Evered baseline: σ = 0.02. |
| `sesolve` | QuTiP closed-system (Schrödinger) solver. Used because noise is shot-to-shot, not intra-shot. |
| `Monte Carlo` | Average over many coherent simulations, each with a different ε sample. |
| `S_Ω` | Protocol sensitivity factor. For π–2π–π: S_Ω = π²/2 ≈ 4.93. |
| `LOCAL_Z_PRODUCT` | `diag(1, -1, -1, 1)`. Local Z phase correction from raw to canonical CZ. |
| `pedersen_fidelity` | Average gate fidelity metric (Pedersen 2007). Works on sub-unitary projected matrices. |

### Current relevant control/data flow

```
src/params.py
    → baseline Ω (from ARC + experimental defaults)
    → src/protocol.py::run_ideal_gate(omega=Ω(1+ε))   [exists, reusable]
    → src/errors/amplitude.py::run_amplitude_noise_gate [new]
    → src/analytical.py::epsilon_amplitude              [new]
    → src/sweeps.py::sweep_amplitude                    [new]
    → scripts/run_sweeps.py                             [extend]
    → scripts/generate_figures.py                       [extend]
    → site/index.html                                   [extend]
```

### Why the key files matter

- **`src/protocol.py`**: `run_ideal_gate(omega)` already runs the full π–2π–π at arbitrary Ω and returns the 4×4 unitary diagonal. The amplitude noise module wraps this in a Monte Carlo loop — it does not need to rebuild the pulse sequence.
- **`src/fidelity.py`**: `pedersen_fidelity(actual, CZ_TARGET, correct_local_z=True)` computes the fidelity including optimal local Z correction. Works on sub-unitary matrices (correctly penalizes leakage to `|r⟩`).
- **`src/errors/blockade.py` and `src/errors/decay.py`**: Establish the patterns for error module structure, return types, input validation, and solver options. The amplitude module should follow these conventions.
- **`src/sweeps.py`**: Has `sweep_decay`, `sweep_blockade` with dataclass rows, CSV I/O, and grid generation. The amplitude sweep adds to this pattern.
- **`site/index.html`**: Lines 1107–1133 contain the `ERROR-SOURCE SECTION TEMPLATE` with the exact structure for new error sections.

### First file a new implementer should open

`src/protocol.py` → `run_ideal_gate()`. Then `src/errors/blockade.py` for the error module pattern. Then `src/sweeps.py` for the sweep pattern.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Implement end-to-end vertical slice for amplitude noise: Monte Carlo over Gaussian Ω fluctuations, analytical formula ε_amp ∝ σ² × S_Ω, parameter sweep over σ_Ω, figure, HTML section, tests.
- **Reproduction steps**: N/A — feature request.
- **Expected behavior**:
  - σ_Ω = 0: F = 1 (1 − F < 10⁻¹⁰)
  - σ_Ω = 2%: numerical matches analytical within factor 2
  - Quadratic scaling: ε(2σ) ≈ 4 × ε(σ)
  - MC uses N ≥ 500 for figure data, N ≤ 100 for tests (< 60s)
  - Sweep: ≥ 20 points, σ = 0.1% to 20%
  - Figure: `figures/fidelity_vs_amplitude_noise.png` with MC error bars, analytical overlay, Evered line
  - HTML: "Amplitude Noise" section
  - `pytest tests/test_amplitude.py` passes
- **Actual behavior**: All specified files are absent.
- **Environment**: `qutip==5.0.4`, `ARC==3.9.0` in `requirements.txt`. Python imports use `PYTHONPATH=.` convention (no `pyproject.toml`).
- **Relevant comments**: No issue comments.

### Repository evidence

| File | What it proves |
|---|---|
| `src/protocol.py` | `run_ideal_gate(omega)` runs the ideal π–2π–π at any `omega` and returns `unitary_diagonal` (4 complex amplitudes) and `unitary` (4×4 matrix). Uses `sesolve` with `nsteps=10000, atol=1e-12, rtol=1e-12`. Already handles all 4 computational basis states. |
| `src/fidelity.py` | `pedersen_fidelity(actual, CZ_TARGET, correct_local_z=True)` works on sub-unitary matrices. `as_operator_matrix` accepts both 1D diagonal vectors and 2D matrices. `LOCAL_Z_PRODUCT = diag(1,-1,-1,1)` matches the ideal raw gate `diag(1,-1,-1,-1)` to CZ. |
| `src/errors/blockade.py` | Established error module pattern: input validation helpers, `_check_positive_finite`, solver options dict, result dict return, self-contained 3-segment pulse sequence. ~200 lines. |
| `src/errors/decay.py` | Established alternative pattern: `@dataclass(frozen=True)` for results, Kraus-based fidelity for dissipative channels, `_zero_decay_result` boundary shortcut. ~250 lines. |
| `src/sweeps.py` | `DecaySweepRow` and `BlockadeSweepRow` dataclasses, CSV I/O with `_write_dataclass_csv`, log-spaced grid generators. Sweep functions take optional grid parameters with sensible defaults. |
| `src/analytical.py` | Has `epsilon_decay_from_gamma`, `epsilon_blockade`, `rydberg_decay_exposure_time`. Each is a one-function formula with input validation and docstring. No `epsilon_amplitude` yet. |
| `scripts/generate_figures.py` | `plot_decay_sweep` and `plot_blockade_sweep` follow a consistent style: `figsize=(7.0, 4.x)`, `constrained_layout=True`, `dpi=220`, analytical curve + numerical scatter + baseline reference line. Uses `_load_or_create_*` pattern with CSV cache. |
| `scripts/run_sweeps.py` | Generates sweep CSVs and prints summary statistics. Uses `PYTHONPATH=.` convention. |
| `tests/test_blockade.py` | Boundary test (`test_blockade_large_u_is_ideal`), baseline agreement test, analytical formula test, sweep shape test. 5 tests, consistent pattern. |
| `tests/test_decay.py` | Same pattern: zero-γ test, baseline agreement (<50%), large-γ degradation, sweep shape, HTML section check. |
| `site/index.html` | Has completed sections for introduction, ideal gate, Rydberg decay, finite blockade. Template at lines 1107–1133 defines the structure for remaining error sections. Footer notes remaining channels. |
| `.gitignore` | `figures/*.png` and `figures/*.csv` are gitignored. Per issue #4 convention, generated figures are force-committed with `git add -f`. |
| `src/params.py` | No amplitude noise default parameter (σ_Ω). The Evered baseline σ_Ω = 0.02 will need to be either added to `RydbergParams` or defined as a module constant in the amplitude error module. |

**Files absent**: `src/errors/amplitude.py`, `tests/test_amplitude.py`, `figures/fidelity_vs_amplitude_noise.png`.

**Existing tests**: HTML scaffold tests, blockade tests, decay tests, ideal gate tests, fidelity tests, hamiltonian tests, params tests. Import requires `PYTHONPATH=.` (no installed package).

## 4. Root-cause candidates

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
|---|---|---|---|---|
| Unimplemented feature slice — the amplitude noise channel simply hasn't been built yet. | All required files are absent. Issue explicitly requests implementation. Blocking issues #1, #2 are closed. | None — this is the direct state of the repo. | Create the files, make tests pass, generate figure and HTML section. | high |
| `run_ideal_gate` may not work correctly at non-default Ω values, causing MC to produce wrong fidelities. | The function uses `omega` as a parameter and passes it to `sesolve`. No obvious validation issue. `test_ideal_gate.py` only tests at the default Ω. | The solver is general-purpose; there's no reason it would fail at Ω(1+ε) for |ε| < 0.2. The Hamiltonians scale linearly with Ω. | Run `run_ideal_gate(omega * 1.02)` and verify fidelity is slightly below 1. | low (bug unlikely) |
| The analytical sensitivity factor S_Ω might not match the numerical result because the formula neglects correlations between pulse segments. | The π–2π–π protocol has three segments that share the same ε. The analytical formula treats them as independent perturbations, which is only approximate. | For small ε, the cross-terms are O(ε⁴) and negligible. The dominant error is from the |01⟩ and |10⟩ channels individually. | Compare analytical S_Ω = π²/2 against numerical sensitivity (finite differences at small ε). Agreement to ~5% would confirm. | medium (likely agrees) |

## 5. Decision

- **Chosen likely cause**: Unimplemented feature slice. The infrastructure is in place; this is the implementation work.

- **Recommended fix**: Implement the amplitude noise vertical slice end-to-end.

### Key design decisions

**1. Simulation approach — Monte Carlo over `run_ideal_gate`**:

For each MC trial:
1. Sample ε ~ N(0, σ²)
2. Set Ω' = Ω(1+ε)
3. Run `run_ideal_gate(omega=Ω')` or a stripped-down variant that returns only the unitary diagonal
4. Compute F = `pedersen_fidelity(unitary, CZ_TARGET, correct_local_z=True)`
5. Average over N trials

This is correct because amplitude noise is a **random-unitary channel** (each shot is coherent). For such channels, the average of per-shot Pedersen fidelities equals the channel fidelity — no Kraus operators needed.

**Optimization**: `run_ideal_gate` collects population traces at 400 points per π-pulse. For MC with N=500, this is wasteful. A stripped-down helper that evolves only the final state (not storing intermediate states) would be ~10× faster. The implementer should write a `_fast_ideal_unitary(omega)` that returns just the 4-element diagonal, using `sesolve` with only `[0, t_final]` as the time grid. Alternatively, compute |01⟩ and |10⟩ amplitudes analytically (they're just `cos(π(1+ε)`) and simulate only |11⟩.

**2. Analytical formula**: `ε_amp = (π²/2) × σ²`

Derivation sketch (the implementer should verify numerically):

The raw ideal gate gives diagonal amplitudes:
- |00⟩: 1 (dark, no dynamics)
- |01⟩: cos(π(1+ε)) = −cos(πε) (2π rotation of atom 2)
- |10⟩: cos(π(1+ε)) = −cos(πε) (π + π rotation of atom 1)
- |11⟩: −(1 − 2sin⁴(πε/2)) ≈ −1 + O(ε⁴) (blocked; quartic sensitivity)

After optimal local-Z correction (which equals `LOCAL_Z_PRODUCT` for this channel):

U_corrected ≈ diag(1, cos(πε), cos(πε), −1 + O(ε⁴))

M = CZ† × U_corrected ≈ diag(1, cos(πε), cos(πε), 1 − O(ε⁴))

Pedersen fidelity:
```
Tr(M†M) = 1 + 2cos²(πε) + 1 − O(ε⁴) ≈ 4 − 2π²ε²
|Tr(M)|² = (1 + 2cos(πε) + 1)² ≈ (4 − π²ε²)² ≈ 16 − 8π²ε²
F = (20 − 10π²ε²) / 20 = 1 − π²ε²/2
```

Averaging over ε ~ N(0, σ²): **ε_amp = ⟨1−F⟩ = (π²/2)σ²**

So **S_Ω = π²/2 ≈ 4.935** for the π–2π–π protocol.

**Key insight**: Only the |01⟩ and |10⟩ channels contribute at O(ε²). The |11⟩ channel (where blockade acts) is insensitive to amplitude noise through O(ε²) — its error is O(ε⁴). This makes physical sense: the blocked 2π-pulse does nothing whether or not Ω is slightly wrong.

**3. Expected baseline value**:

At σ_Ω = 0.02 (Evered 2023):
```
ε_amp = (π²/2)(0.02)² ≈ 4.93 × 4 × 10⁻⁴ ≈ 1.97 × 10⁻³
```

This is ~0.2% infidelity — **larger than Rydberg decay** (ε_decay ≈ 5.85 × 10⁻⁴) at the project baseline parameters. The issue describes this as "typically negligible at 2%," which may apply to optimized protocols with composite pulses but not to the simple π–2π–π. The implementer should verify this numerically and adjust the HTML takeaway accordingly.

For comparison:
| Error source | Baseline infidelity |
|---|---|
| Rydberg decay | ~5.9 × 10⁻⁴ |
| Finite blockade | ~1.4 × 10⁻⁶ |
| Amplitude noise (predicted) | ~2.0 × 10⁻³ |

**4. Handling σ_Ω as a parameter**:

`src/params.py` has no amplitude noise parameter. Two options:
1. Add `sigma_omega: float = 0.02` to `RydbergParams` (requires modifying `get_rydberg_params`).
2. Define `DEFAULT_SIGMA_OMEGA = 0.02` as a module constant in `src/errors/amplitude.py`.

**Recommendation**: Option 2 (module constant). σ_Ω is an experimental choice, not an ARC-derivable atomic physics parameter. It doesn't belong in `get_rydberg_params()`. Follow the precedent of `n_steps_per_pi` being a module-level default rather than a params field.

- **Rejected alternatives and why**:
  - *Use `mesolve` with a stochastic Hamiltonian*: Wrong for this error. Amplitude noise is **shot-to-shot** (each experimental run has a different Ω), not **within-shot** (Ω doesn't fluctuate during a single gate). `sesolve` per shot + averaging is the correct model.
  - *Build a new pulse sequence simulator from scratch*: Unnecessary — `run_ideal_gate` already does the job. Wrap, don't rewrite.
  - *Use a process-matrix approach (Choi/Kraus)*: Overkill for a random-unitary channel. The per-shot fidelity average is mathematically equivalent and much simpler. The decay module uses Kraus because Lindblad dynamics are genuinely non-unitary within each shot.

- **Assumptions**:
  - The same ε applies to all three pulses in each shot (global intensity fluctuation, not per-pulse). This matches the issue specification: "scale Ω → Ω(1+ε) for the entire pulse sequence."
  - The infinite-blockade model is used (same as `run_ideal_gate`). Finite-blockade cross-terms with amplitude noise are a combined-error effect (issue #8).
  - Units: σ_Ω is dimensionless (fractional noise). Ω is in rad/μs.

- **Maintainer decision needed**: None. The implementation path is clear from existing patterns.

## 6. Implementation handoff

### Ordered implementation steps

1. **Add `epsilon_amplitude` to `src/analytical.py`**:
   ```python
   def epsilon_amplitude(sigma_omega: float, sensitivity: float | None = None) -> float:
       """Leading amplitude-noise infidelity: (π²/2) × σ².
       
       ``sensitivity`` defaults to π²/2 (the π–2π–π protocol value).
       """
       if sensitivity is None:
           sensitivity = math.pi**2 / 2.0
       return sensitivity * sigma_omega**2
   ```
   Verify: `epsilon_amplitude(0.02) ≈ 1.97e-3`.

2. **Create `src/errors/amplitude.py`**:
   - Module constant: `DEFAULT_SIGMA_OMEGA = 0.02`, `DEFAULT_N_SAMPLES = 500`.
   - A fast unitary helper: `_fast_ideal_diagonal(omega)` that returns only the 4-element complex diagonal without population traces. Uses `sesolve` with `tlist=[0, t_final]` for each segment.
   - Main function: `run_amplitude_noise_gate(omega, sigma_omega, n_samples, *, seed=None)` → dataclass result with `average_fidelity`, `fidelity_std`, `fidelity_stderr`, `sample_fidelities` (for error bars), `n_samples`.
   - The seed parameter is important for reproducibility in tests.
   - MC loop: for k in range(n_samples): sample ε, compute diagonal at Ω(1+ε), compute Pedersen fidelity with `correct_local_z=True`, collect results.
   - Return mean, std, standard error of the mean (std/√N).
   - Boundary shortcut: if σ = 0, return F = 1 exactly (no MC needed).

3. **Add amplitude sweep to `src/sweeps.py`**:
   - `AmplitudeSweepRow` dataclass: `sigma_omega`, `numerical_fidelity`, `numerical_error`, `numerical_stderr`, `analytical_error`, `analytical_fidelity`.
   - `amplitude_sigma_grid(num_points=25, min_sigma=0.001, max_sigma=0.20)` → log-spaced σ array.
   - `sweep_amplitude(omega, sigmas, n_samples=500, seed=None)` → list of `AmplitudeSweepRow`.
   - CSV read/write helpers following the `_write_dataclass_csv` pattern.

4. **Create `tests/test_amplitude.py`**:
   ```python
   def test_amplitude_zero_noise_is_ideal():
       """σ = 0 → F = 1 exactly."""
       result = run_amplitude_noise_gate(omega=DEFAULT_OMEGA, sigma_omega=0.0)
       assert 1.0 - result.average_fidelity < 1e-10

   def test_amplitude_baseline_matches_analytical():
       """σ = 2%: numerical ≈ analytical within factor 2."""
       result = run_amplitude_noise_gate(omega=DEFAULT_OMEGA, sigma_omega=0.02, 
                                          n_samples=100, seed=42)
       analytical = epsilon_amplitude(0.02)
       numerical = 1.0 - result.average_fidelity
       assert abs(numerical - analytical) / analytical < 1.0  # factor 2

   def test_amplitude_quadratic_scaling():
       """ε(2σ) ≈ 4 × ε(σ)."""
       r1 = run_amplitude_noise_gate(sigma_omega=0.05, n_samples=100, seed=42)
       r2 = run_amplitude_noise_gate(sigma_omega=0.10, n_samples=100, seed=42)
       ratio = (1 - r2.average_fidelity) / (1 - r1.average_fidelity)
       assert 3.0 < ratio < 5.0  # should be ~4

   def test_epsilon_amplitude_formula():
       assert np.isclose(epsilon_amplitude(0.1), math.pi**2 * 0.01 / 2)
   ```
   - Use N ≤ 100 and fixed seed for tests to keep runtime < 60s.
   - Tolerances on MC results should be loose (factor 2) to avoid flaky tests.

5. **Update `scripts/run_sweeps.py`**: Add amplitude sweep call, save `figures/amplitude_sweep.csv`.

6. **Update `scripts/generate_figures.py`**: Add `plot_amplitude_sweep()`:
   - x-axis: σ_Ω (%) — multiply by 100 for display
   - y-axis: gate infidelity (1 − F), log scale
   - Numerical scatter with error bars (stderr from MC)
   - Analytical curve: ε = (π²/2)σ²
   - Vertical dashed line at σ = 2% (Evered baseline)
   - Style matching existing plots (7.0 × 4.x inches, 220 dpi, constrained_layout)

7. **Add "Amplitude Noise" section to `site/index.html`**:
   - Insert before the footer, after the blockade section (or after whichever error sections exist at merge time)
   - Use the template at lines 1107–1133
   - Section id: `sec-error-amplitude`
   - Physics: laser intensity drifts → Ω varies shot-to-shot → pulse areas wrong → population doesn't fully return to ground state
   - Formula: ε_amp = (π²/2)σ² in MathJax
   - Collapsible derivation: key insight is that only |01⟩ and |10⟩ contribute at O(σ²), because the blocked |11⟩ channel is insensitive to Ω changes
   - Figure embed
   - Takeaway: "Amplitude noise scales quadratically — small at σ < 1%, but can become the dominant error at σ ~ 5%."

8. **Force-commit the figure**: `git add -f figures/fidelity_vs_amplitude_noise.png`

### Files likely to change

| File | Action |
|---|---|
| `src/errors/amplitude.py` | New |
| `src/analytical.py` | Extend (add `epsilon_amplitude`) |
| `src/sweeps.py` | Extend (add `AmplitudeSweepRow`, `sweep_amplitude`, CSV I/O) |
| `scripts/run_sweeps.py` | Extend (add amplitude sweep call) |
| `scripts/generate_figures.py` | Extend (add `plot_amplitude_sweep`) |
| `tests/test_amplitude.py` | New |
| `site/index.html` | Modify (add Amplitude Noise section, update footer) |
| `figures/fidelity_vs_amplitude_noise.png` | Generated |
| `figures/amplitude_sweep.csv` | Generated |

### Scope boundaries / what not to touch

- Do NOT add Lindblad/dissipative dynamics. Amplitude noise is shot-to-shot coherent noise.
- Do NOT implement Doppler, scattering, or combined budget (issues #5, #6, #8).
- Do NOT modify `src/protocol.py` — reuse `run_ideal_gate` or write a standalone fast helper in the amplitude module.
- Do NOT modify `src/params.py` — define σ_Ω defaults in the amplitude module.
- Do NOT modify `src/hamiltonian.py`.

### Regression test shape

See step 4 above. Key tests:
1. **Boundary**: σ = 0 → F = 1 (1−F < 10⁻¹⁰)
2. **Baseline agreement**: σ = 2% → numerical within factor 2 of analytical
3. **Quadratic scaling**: ε(2σ)/ε(σ) ≈ 4
4. **Analytical formula unit test**: known σ → known ε

### Validation commands

```bash
# Targeted amplitude tests (should pass in < 60s with N ≤ 100)
PYTHONPATH=. python -m pytest -xvs tests/test_amplitude.py

# Generate sweep data + figure
PYTHONPATH=. python scripts/run_sweeps.py
PYTHONPATH=. python scripts/generate_figures.py

# Full regression (all existing + new tests)
PYTHONPATH=. python -m pytest -q

# Quick numerical check of the analytical sensitivity
PYTHONPATH=. python -c "
from src.protocol import run_ideal_gate
from src.fidelity import pedersen_fidelity, CZ_TARGET
import numpy as np
omega = 2*np.pi*4.0
eps = 0.001
result = run_ideal_gate(omega*(1+eps))
F = pedersen_fidelity(result['unitary'], CZ_TARGET, correct_local_z=True)
print(f'At eps={eps}: 1-F = {1-F:.6e}, predicted = {np.pi**2*eps**2/2:.6e}')
"
```

Note: the exact Python executable depends on the environment. The existing tests use `PYTHONPATH=.` convention from the repo root.

### Expected pre-fix signal

- `import src.errors.amplitude` → `ModuleNotFoundError`
- `pytest tests/test_amplitude.py` → file not found
- `figures/fidelity_vs_amplitude_noise.png` → absent
- `site/index.html` → no `sec-error-amplitude` section

### Expected post-fix signal

- All new and existing tests pass
- σ_Ω = 0 → infidelity < 10⁻¹⁰
- σ_Ω = 2% → numerical infidelity ≈ 2 × 10⁻³ (analytical: 1.97 × 10⁻³)
- Quadratic scaling verified
- Sweep: ≥ 20 points, σ from 0.1% to 20%
- Figure exists with scatter + error bars + analytical curve + Evered line
- HTML section present, scaffold tests pass

### Negative/edge validation

- σ < 0 → should raise `ValueError`
- σ = 0 → F = 1 exactly (boundary, no MC needed)
- Very large σ (e.g., σ = 1.0 = 100%) → fidelity should drop significantly; perturbative formula breaks down (expected)
- N_samples = 1 → should work but give meaningless statistics; useful for debugging
- Seed reproducibility: same seed → same result

## 7. Risks, edge cases, rollback

### Risks

| Risk | Mitigation |
|---|---|
| **MC statistical noise in tests**: With N=100, the MC estimate has substantial variance (~10% relative error on fidelity). Tests may be flaky. | Use fixed seed for reproducibility. Use loose tolerances (factor 2) for agreement checks. The quadratic-scaling test is more robust (ratio of two MC estimates, systematic errors cancel). |
| **Analytical formula may disagree with numerics at large σ**: The perturbative formula ε = (π²/2)σ² is valid only for σ ≪ 1. At σ = 20%, higher-order terms matter. | This is expected and physically interesting. The sweep figure should show the analytical curve diverging from numerical points at large σ. Restrict the "match within factor 2" acceptance criterion to σ ≤ 5%. |
| **Population leakage to `\|r⟩` makes the projected unitary sub-unitary**: At each shot with ε ≠ 0, the atom doesn't fully return to `\|g⟩`. The 4×4 projected matrix has column norms < 1. | `pedersen_fidelity` handles this correctly — `Tr(M†M) < d` naturally penalizes leakage. No special handling needed. |
| **HTML section ordering**: If issues #5 and #6 are implemented concurrently, the HTML section insertion point may conflict. | Insert using the section ID convention (`sec-error-amplitude`). The template structure makes conflicts resolvable by inspection. |
| **`run_ideal_gate` overhead**: The function collects population traces (400 points × 4 segments) which are discarded in MC. | Write a `_fast_ideal_diagonal(omega)` that returns only the final-state amplitudes. This is a ~5× speedup per shot, making N=500 practical in seconds. |

### Edge cases

- **σ → 0**: Must reproduce F = 1 exactly. The MC loop should short-circuit at σ = 0 (return identity gate).
- **σ → ∞**: Gate completely fails. Perturbative formula is invalid. Numerical simulation should still run and give low fidelity.
- **σ ≈ 0.01 (1%)**: Good perturbative regime. Analytical and numerical should agree closely.
- **σ ≈ 0.2 (20%)**: Boundary of the sweep range. Perturbative formula starts breaking down — visible in the figure.
- **Negative ε samples**: The Gaussian distribution includes negative ε (Ω slightly below nominal). This is correct and handled naturally.
- **Very negative ε (ε < −1, so Ω < 0)**: At large σ, a rare sample could give Ω < 0. The `sesolve` Hamiltonian with negative Ω is still valid (it reverses the rotation direction). However, `run_ideal_gate` validates `omega > 0`. Either clamp samples to ε > −1, or handle this edge case. At σ = 0.20, P(ε < −1) ≈ 3 × 10⁻⁷ — negligible.

### Rollback/revert notes

This slice is purely additive. Revert = remove new files + revert `src/analytical.py`, `src/sweeps.py`, script, and `site/index.html` modifications. No impact on ideal gate or other error channels.

## 8. Open questions

### 1. Should the fast unitary helper be in `amplitude.py` or extracted to a shared module?

- **What was found**: `run_ideal_gate` in `protocol.py` computes the full gate with population traces. The amplitude MC needs only the final unitary diagonal. A fast helper would skip intermediate time steps.
- **Why it matters**: N=500 shots at ~0.5s each = ~4 minutes with full `run_ideal_gate`. With a fast helper, ~30 seconds.
- **Options**:
  1. Write `_fast_ideal_diagonal(omega)` in `amplitude.py` (self-contained).
  2. Add a `fast=True` parameter to `run_ideal_gate` in `protocol.py`.
  3. Extract to `src/solver_utils.py`.
- **Recommended default**: Option 1. Keep the amplitude module self-contained, following the blockade.py precedent of duplicating solver logic rather than modifying shared modules. If issue #5 (Doppler MC) needs the same helper, refactor then.

### 2. Should the analytical formula use a hard-coded S_Ω or compute it numerically?

- **What was found**: S_Ω = π²/2 is the analytically derived sensitivity for π–2π–π. This is a clean closed-form result.
- **Why it matters**: A hard-coded formula is simpler and faster. A numerically computed S_Ω (via finite differences) would be more general but adds complexity.
- **Recommended default**: Hard-code S_Ω = π²/2 in `epsilon_amplitude`. Add a docstring noting the derivation. Include an optional `sensitivity` parameter for future protocol variants. Verify numerically in the test suite.

### 3. Should the HTML takeaway say "negligible" or "non-negligible" at σ = 2%?

- **What was found**: The issue suggests "typically negligible at σ ~ 2%." But the analytical formula gives ε_amp ≈ 2 × 10⁻³ at σ = 2%, which is **larger** than the Rydberg decay error (5.9 × 10⁻⁴) at the project baseline.
- **Why it matters**: The HTML should state correct physics. Calling a dominant error source "negligible" would mislead the reader.
- **Recommended default**: Let the numerical results speak. If ε_amp ≈ 2 × 10⁻³ is confirmed, the takeaway should be: "At σ_Ω = 2%, amplitude noise contributes ~0.2% infidelity — comparable to or larger than Rydberg decay for this protocol. More sophisticated pulse sequences (composite pulses, echo techniques) can reduce this sensitivity."
