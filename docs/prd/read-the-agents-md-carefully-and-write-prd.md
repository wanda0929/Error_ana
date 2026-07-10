# PRD: Rydberg CZ Gate Error Analysis for Rb-Rb

## 1. Overview

### 1.1 Problem statement

Rydberg blockade CZ gates are the workhorse entangling operation for neutral-atom quantum computers. Understanding how each experimental imperfection degrades gate fidelity — and validating analytical scaling laws against exact numerics — is essential for identifying the dominant error sources and guiding hardware improvements.

No self-contained, pedagogical resource exists that both derives the analytical error formulas *and* validates them numerically in one place, with clear visualizations accessible to physicists outside the Rydberg community.

### 1.2 Solution

Build a numerical simulation of the π–2π–π Rydberg blockade CZ gate for two ⁸⁷Rb atoms, study five error sources individually and combined, and present the results in an interactive HTML explainer.

Each error source gets:
- A closed-form analytical scaling law (from literature)
- A full QuTiP numerical simulation
- A plot overlaying both, with the Evered 2023 experimental operating point marked

### 1.3 Deliverables

| # | Artifact | Format | Purpose |
|---|---|---|---|
| 1 | Simulation library | Python (`src/`) | Reproducible numerical results |
| 2 | Parameter sweep data | NumPy arrays / CSV | Raw data for all plots |
| 3 | Publication-quality figures | PNG (`figures/`) | Embedded in HTML |
| 4 | Interactive HTML explainer | Single HTML file (`site/index.html`) | **Primary deliverable** — reader-facing |
| 5 | Test suite | pytest (`tests/`) | Correctness guarantees |

### 1.4 Success criteria

| Criterion | Metric | Threshold |
|---|---|---|
| Ideal gate correctness | Infidelity of zero-error simulation | < 10⁻¹⁰ |
| Analytical-numerical agreement | Max relative error at Evered 2023 point | < factor of 2 for all 5 error sources |
| Combined budget validity | \|ε_sum − ε_numerical\| / ε_numerical | < 20% (additive approximation holds) |
| HTML completeness | All 8 content sections present with figures | 100% |
| Test coverage | Every public function in `src/` has ≥1 test | 100% |

---

## 2. System and physics specification

### 2.1 Physical system

- **Atoms**: Two ⁸⁷Rb atoms in optical tweezers
- **Levels per atom**: |g⟩ (ground, qubit |1⟩) and |r⟩ (Rydberg |70S₁/₂⟩)
- **Dark state**: Qubit |0⟩ does not couple to the laser
- **Excitation**: Two-photon (780 nm + 480 nm, counter-propagating)
- **Hilbert space**: 4-dimensional {|gg⟩, |gr⟩, |rg⟩, |rr⟩}

### 2.2 Baseline parameters

| Parameter | Symbol | Value | Source |
|---|---|---|---|
| Rydberg state | \|n,l,j⟩ | \|70, 0, 0.5⟩ | Evered 2023 |
| Rabi frequency | Ω/2π | ~4 MHz | Evered 2023 |
| Rydberg lifetime (300K) | τ | ~230 μs | ARC at runtime |
| C₆ coefficient | C₆ | ~860 GHz·μm⁶ | ARC at runtime |
| Interatomic distance | R | ~6 μm | U >> Ω regime |
| Blockade strength | U = C₆/R⁶ | >> Ω | Blockade regime |
| Atom temperature | T | 10 μK | Evered 2023 |
| Intermediate detuning | Δ_p/2π | ~1 GHz | Evered 2023 |
| Amplitude noise | σ_Ω | 2% | Evered 2023 |
| Effective wavevector | k_eff | 0.80 μm⁻¹ | Counter-propagating |

All physical constants (τ, C₆, Γ_e) from ARC (`arc.Rubidium87`) — never hard-coded.

### 2.3 Gate protocol

**π–2π–π sequential** (Jaksch 2000, Saffman 2010, Zhang 2012):

1. π-pulse on atom 1: |g⟩₁ → |r⟩₁
2. 2π-pulse on atom 2: full rotation (unblocked) or blocked (if atom 1 in |r⟩)
3. π-pulse on atom 1: |r⟩₁ → |g⟩₁

**Target unitary**: CZ = diag(1, 1, 1, −1)

**Total Rydberg occupation time**: T_R = 7π/(4Ω)

### 2.4 Fidelity metric

Pedersen average gate fidelity:

```
F_avg = [Σ_k Tr(M_k†M_k) + Σ_k |Tr(M_k)|²] / [d(d+1)]
```

where M_k = U_target† G_k, {G_k} are Kraus operators, d = 4.

Properties: global-phase-invariant, sensitive to both population and phase errors, standard in the field.

---

## 3. Error sources

### 3.1 Error channel summary

| # | Error source | Type | Model | Sweep parameter |
|---|---|---|---|---|
| 1 | Rydberg decay | Incoherent | Lindblad: L = √γ \|g⟩⟨r\| | Decay rate γ (or lifetime τ) |
| 2 | Finite blockade | Coherent | Include \|rr⟩ at energy U | Ratio U/Ω |
| 3 | Doppler dephasing | Incoherent (ensemble) | Monte Carlo over thermal velocities | Temperature T |
| 4 | Intermediate-state scattering | Incoherent | Lindblad effective dephasing | Detuning Δ_p |
| 5 | Amplitude noise | Incoherent (ensemble) | Monte Carlo over Ω fluctuations | Noise σ_Ω |

### 3.2 Analytical scaling laws

| Error source | Formula | Scaling |
|---|---|---|
| Rydberg decay | ε_decay = (1/τ) × 7π/(4Ω) | ∝ 1/(Ω·τ) |
| Finite blockade | ε_block = Ω²/(8U²) | ∝ (Ω/U)² |
| Doppler | ε_Doppler ∝ k_eff² k_B T t² / m | ∝ T × t_gate² |
| Intermediate scattering | ε_scatter ∝ Γ_e/Δ_p × (pulse time) | ∝ 1/Δ_p |
| Amplitude noise | ε_amp ∝ σ_Ω² × S_Ω | ∝ σ² |

### 3.3 Combined analysis

All five errors turned on simultaneously. Compare:
- ε_numerical (full simulation)
- ε_sum = Σᵢ εᵢ (additive approximation)

Discrepancy reveals cross-channel interference.

---

## 4. User stories

### 4.1 Simulation core

| ID | Story | Acceptance criteria |
|---|---|---|
| S-01 | As a developer, I can query ARC for Rb87 Rydberg state parameters (τ, C₆) given (n, l, j) | `params.py` returns τ ≈ 230 μs and C₆ ≈ 860 GHz·μm⁶ for \|70,0,0.5⟩ at 300K |
| S-02 | As a developer, I can construct the 2-atom Hamiltonian with per-atom addressing and blockade interaction | `hamiltonian.py` returns a 4×4 QuTiP Qobj; matrix elements match hand calculation |
| S-03 | As a developer, I can run the π–2π–π protocol on the ideal system and get F = 1 | `protocol.py` ideal run yields 1 − F < 10⁻¹⁰ |
| S-04 | As a developer, I can compute Pedersen average gate fidelity from simulation output | `fidelity.py` returns F = 1.0 for CZ vs CZ, F = 0.5 for CZ vs identity |
| S-05 | As a developer, I can extract the 4×4 process matrix or Kraus operators from density matrix outputs | Process matrix is trace-preserving; Kraus ops satisfy Σ G_k†G_k = I for ideal case |

### 4.2 Error modules

| ID | Story | Acceptance criteria |
|---|---|---|
| E-01 | As a developer, I can simulate Rydberg decay by adding a Lindblad collapse operator | At γ=0: F=1. At Evered 2023 γ: numerical ε_decay matches analytical to within 50% |
| E-02 | As a developer, I can simulate finite blockade by including \|rr⟩ with finite U | At U→∞ (U/Ω=1000): F→1. At Evered 2023 U/Ω: numerical ε_block matches Ω²/(8U²) to within 50% |
| E-03 | As a developer, I can simulate Doppler dephasing via Monte Carlo velocity sampling | At T=0: F=1. At T=10μK: numerical ε_Doppler matches analytical to within factor 2 |
| E-04 | As a developer, I can simulate intermediate-state scattering as effective Lindblad dephasing | At Δ_p→∞: F=1. At Δ_p=1GHz: numerical ε_scatter matches analytical to within factor 2 |
| E-05 | As a developer, I can simulate amplitude noise via Monte Carlo Ω sampling | At σ_Ω=0: F=1. At σ_Ω=2%: numerical ε_amp matches analytical to within factor 2 |
| E-06 | As a developer, I can run all 5 errors simultaneously and compare against additive sum | Combined infidelity within 20% of Σ εᵢ at Evered 2023 point |

### 4.3 Analytical formulas

| ID | Story | Acceptance criteria |
|---|---|---|
| A-01 | As a developer, I can compute ε_decay analytically given (Ω, τ) | Returns value matching Zhang 2012 formula |
| A-02 | As a developer, I can compute ε_block analytically given (Ω, U) | Returns value matching Saffman 2010 formula |
| A-03 | As a developer, I can compute ε_Doppler analytically given (k_eff, T, m, t_gate) | Returns value consistent with Evered 2023 reported Doppler error |
| A-04 | As a developer, I can compute ε_scatter analytically given (Γ_e, Δ_p, Ω) | Returns value matching de Léséleuc 2018 formula |
| A-05 | As a developer, I can compute ε_amp analytically given (σ_Ω, protocol) | Returns value consistent with quadratic scaling |

### 4.4 Parameter sweeps

| ID | Story | Acceptance criteria |
|---|---|---|
| P-01 | As a developer, I can sweep decay rate and get fidelity vs γ data | Output: arrays of (γ_values, F_numerical, F_analytical) with ≥20 points |
| P-02 | As a developer, I can sweep U/Ω ratio and get fidelity vs blockade data | Output: arrays of (U_Ω_ratio, F_numerical, F_analytical) with ≥20 points |
| P-03 | As a developer, I can sweep temperature and get fidelity vs T data | Output: arrays of (T_values, F_numerical, F_analytical) with ≥20 points |
| P-04 | As a developer, I can sweep Δ_p and get fidelity vs detuning data | Output: arrays of (Delta_p_values, F_numerical, F_analytical) with ≥20 points |
| P-05 | As a developer, I can sweep σ_Ω and get fidelity vs noise data | Output: arrays of (sigma_values, F_numerical, F_analytical) with ≥20 points |

### 4.5 Figures

| ID | Story | Acceptance criteria |
|---|---|---|
| F-01 | As a reader, I can see population dynamics during the ideal gate | Plot: Rydberg population vs time for each basis state {00,01,10,11}. 4 curves, correct pulse timing visible |
| F-02 | As a reader, I can see fidelity vs decay rate with analytical overlay | Scatter (numerical) + solid line (analytical) + dashed vertical (Evered 2023). Axes labeled with units |
| F-03 | As a reader, I can see fidelity vs U/Ω with analytical overlay | Same format as F-02 |
| F-04 | As a reader, I can see fidelity vs temperature with analytical overlay | Same format as F-02 |
| F-05 | As a reader, I can see fidelity vs Δ_p with analytical overlay | Same format as F-02 |
| F-06 | As a reader, I can see fidelity vs σ_Ω with analytical overlay | Same format as F-02 |
| F-07 | As a reader, I can see the combined error budget as a bar chart | Bars: individual ε_i values + total ε_sum + numerical ε_combined. Clear legend |

### 4.6 HTML explainer

| ID | Story | Acceptance criteria |
|---|---|---|
| H-01 | As a reader, I can understand what a Rydberg CZ gate is from the introduction alone | Section defines: qubit encoding, Rydberg state, blockade, CZ truth table. No undefined jargon |
| H-02 | As a reader, I can see the ideal gate working before errors are introduced | Section shows truth table + population dynamics plot + F=1 confirmation |
| H-03 | As a reader, I can understand each error source independently | 5 sections, each with: physics explanation, formula, derivation sketch (collapsible), plot, takeaway |
| H-04 | As a reader, I can see the combined error budget and whether additive approximation holds | Bar chart + discussion of agreement/disagreement |
| H-05 | As a reader, I can see a summary table of all errors at the Evered 2023 point | Table with columns: error source, analytical ε, numerical ε, dominant scaling |
| H-06 | As a reader, I can expand collapsible derivation sections for more detail | Click-to-expand works; derivations show key steps, not just citations |
| H-07 | As a reader, the HTML renders correctly with MathJax equations | All LaTeX renders; no raw `\frac{}{}` visible on page |
| H-08 | As a reader, the HTML is readable on mobile | Responsive layout; figures scale; no horizontal scroll |

### 4.7 Testing

| ID | Story | Acceptance criteria |
|---|---|---|
| T-01 | As a developer, the ideal gate test catches any Hamiltonian or protocol bug | `test_ideal_gate.py` asserts 1 − F < 10⁻¹⁰; fails loudly if broken |
| T-02 | As a developer, each error module has a zero-error → F=1 boundary test | 5 tests, one per error module, all pass with F > 1 − 10⁻¹⁰ |
| T-03 | As a developer, the fidelity function has known-input tests | CZ vs CZ → 1.0; CZ vs I → known value; random unitary → bounded value |
| T-04 | As a developer, analytical formulas reproduce published values | Each formula tested against at least one literature-reported number |
| T-05 | As a developer, I can run the full test suite in < 60 seconds | `pytest tests/` completes within timeout (no Monte Carlo with N>100 in tests) |

---

## 5. Architecture

### 5.1 Technology stack

| Layer | Choice | Rationale |
|---|---|---|
| Simulation engine | QuTiP (mesolve, sesolve) | Industry standard for open quantum systems; no wrapper libraries |
| Atomic parameters | ARC (arc.Rubidium87) | Authoritative Rydberg data; no hard-coded constants |
| Fidelity | Custom implementation | Pedersen formula is 10 lines; no need for external dependency |
| Plotting | Matplotlib | Publication-quality figures; full control |
| HTML | Hand-written + MathJax CDN | Single-file, no build step, no framework |
| Testing | pytest | Standard Python testing |
| Language | Python ≥ 3.10 | Type hints, match statements if needed |

### 5.2 Module dependency graph

```
params.py (ARC queries)
    ↓
hamiltonian.py (builds QuTiP operators)
    ↓
protocol.py (runs π–2π–π time evolution)
    ↓                    ↓
errors/*.py          fidelity.py
(modify protocol       (scores output
 with noise)           against CZ)
    ↓                    ↓
sweeps.py (orchestrates: error module × sweep range → fidelity array)
    ↓
analytical.py (overlays closed-form curves)
    ↓
generate_figures.py → figures/*.png → site/index.html
```

### 5.3 Data flow for one error sweep

```
Input: error_type="decay", sweep_range=np.logspace(-6, -3, 30)
  │
  ├─ For each γ in sweep_range:
  │     1. params.py → get base parameters
  │     2. hamiltonian.py → build H
  │     3. protocol.py + errors/decay.py → run π–2π–π with Lindblad L=√γ|g⟩⟨r|
  │     4. fidelity.py → compute F_avg from output ρ(T)
  │     5. Store (γ, F_avg)
  │
  ├─ analytical.py → compute ε_decay(γ) = γ × 7π/(4Ω) for same range
  │
  └─ Output: (γ_array, F_numerical_array, F_analytical_array)
       → passed to generate_figures.py
```

### 5.4 Testing seams

Every module boundary is a **pure-function seam**: call function with known inputs, assert on outputs. No mocking, no infrastructure, no network.

| Seam | Input | Output | Oracle |
|---|---|---|---|
| `params.py` | (species, n, l, j) | (τ, C₆, Γ_e) | ARC documentation / literature |
| `hamiltonian.py` | (Ω, Δ, U, atom_index) | 4×4 Qobj | Hand-computed matrix elements |
| `protocol.py` | (ideal params) | final state/unitary | F = 1 (exact) |
| `fidelity.py` | (U_actual, U_target) | F_avg float | Known matrix pair → known F |
| `errors/*.py` | (params + error_strength) | density matrix | At error=0: F=1; at finite error: matches analytical |
| `analytical.py` | (physical params) | ε float | Literature-reported values |

---

## 6. Milestones and execution order

### Phase 1: Foundation (simulation core)

| Step | Delivers | Depends on | Validates |
|---|---|---|---|
| 1.1 | `params.py` | ARC installed | Known τ, C₆ for |70S₁/₂⟩ |
| 1.2 | `hamiltonian.py` | 1.1 | Matrix elements match hand calc |
| 1.3 | `protocol.py` (ideal) | 1.2 | F = 1 to machine precision |
| 1.4 | `fidelity.py` | — | Known-input tests pass |
| 1.5 | `test_ideal_gate.py` | 1.3, 1.4 | Green CI |

**Exit criterion**: `pytest tests/test_ideal_gate.py` passes. Ideal CZ confirmed.

### Phase 2: Error modules (one at a time)

| Step | Delivers | Depends on | Validates |
|---|---|---|---|
| 2.1 | `errors/decay.py` + `analytical.py` (decay) | Phase 1 | Zero→F=1; Evered point matches formula |
| 2.2 | `errors/blockade.py` + `analytical.py` (blockade) | Phase 1 | Zero→F=1; Evered point matches formula |
| 2.3 | `errors/doppler.py` + `analytical.py` (Doppler) | Phase 1 | Zero→F=1; Evered point matches formula |
| 2.4 | `errors/scattering.py` + `analytical.py` (scatter) | Phase 1 | Zero→F=1; Evered point matches formula |
| 2.5 | `errors/amplitude.py` + `analytical.py` (amplitude) | Phase 1 | Zero→F=1; Evered point matches formula |
| 2.6 | Combined simulation | 2.1–2.5 | Additive approx within 20% |

**Exit criterion**: All error modules pass boundary tests. Analytical-numerical agreement confirmed at Evered 2023 point.

### Phase 3: Sweeps and figures

| Step | Delivers | Depends on |
|---|---|---|
| 3.1 | `sweeps.py` | Phase 2 |
| 3.2 | `scripts/run_sweeps.py` → sweep data | 3.1 |
| 3.3 | `scripts/generate_figures.py` → `figures/*.png` | 3.2 |

**Exit criterion**: 7 publication-quality PNGs generated. Each has analytical + numerical + Evered reference.

### Phase 4: HTML explainer

| Step | Delivers | Depends on |
|---|---|---|
| 4.1 | HTML structure + intro + ideal gate section | Phase 1 figures |
| 4.2 | 5 error source sections with embedded figures | Phase 3 |
| 4.3 | Combined budget section + summary table | Phase 3 |
| 4.4 | Collapsible derivations + MathJax + responsive CSS | 4.1–4.3 |

**Exit criterion**: `site/index.html` opens in a browser, all figures visible, MathJax renders, mobile-responsive.

---

## 7. Constraints and non-goals

### Constraints

- **No external gate libraries**: Raw QuTiP only. No Pulser, no Cirq, no custom wrappers.
- **No hard-coded physics**: All atomic constants from ARC at runtime.
- **Single HTML file**: No build system, no npm, no framework. MathJax via CDN is the only external dependency.
- **Sequential protocol only**: No GRAPE, no optimized pulses. The π–2π–π protocol is the only gate studied.
- **Python ≥ 3.10**: Use type hints throughout.

### Non-goals

- Optimizing the pulse (this is an analysis project, not an optimization project)
- Multi-qubit gates beyond two qubits
- Other atomic species (Cs, Yb, Sr)
- Phase noise / laser linewidth (deferred — subdominant at 3×10⁻⁴)
- m_F sublevel structure (would require expanding Hilbert space beyond 4 levels)
- Atom loss as a separate channel from decay
- Experimental comparison beyond Evered 2023 operating point
- Deployment / hosting of the HTML (just a local file)

---

## 8. Risks

| Risk | Impact | Mitigation |
|---|---|---|
| ARC returns unexpected values for |70S₁/₂⟩ | Wrong baseline parameters | Cross-check ARC output against Evered 2023 supplemental |
| Analytical Doppler formula has poor literature consensus | A-03 hard to validate | Use Evered 2023 reported value as ground truth; flag if >2× off |
| Monte Carlo errors (Doppler, amplitude) converge slowly | Noisy numerical dots on plots | Use N≥500 samples for final figures; N≤100 for tests (speed) |
| QuTiP mesolve numerics for very small γ may not resolve decay | ε_decay underestimated at small γ | Increase n_timesteps; verify convergence by doubling |
| Intermediate-state scattering model (effective dephasing) is approximate | ε_scatter off by constant factor | Acceptable — the effective model is standard (de Léséleuc 2018). Flag as approximation in HTML |
| Additive error approximation breaks down | Combined ε ≠ Σεᵢ | This is a result, not a failure. Report and discuss in HTML |

---

## 9. Dependencies

### Python packages

```
qutip >= 5.0
arc-alkali-rydberg-calculator >= 3.0
numpy
scipy
matplotlib
pytest
```

### External

- MathJax CDN (loaded by HTML at view time; no build dependency)
- ARC database (bundled with the arc package; no network needed after install)

---

## 10. Glossary

| Term | Definition |
|---|---|
| CZ gate | Controlled-Z: two-qubit gate that applies a π phase to \|11⟩ and leaves other basis states unchanged |
| Rydberg state | Highly excited atomic state (large principal quantum number n) with strong long-range interactions |
| Blockade | When one atom in \|r⟩ shifts the Rydberg level of a nearby atom so far off-resonance that it cannot be excited |
| π-pulse | A laser pulse of area π that drives complete population transfer between two levels |
| Lindblad equation | Master equation for open quantum systems: ρ̇ = −i[H,ρ] + Σ(LρL† − ½{L†L,ρ}) |
| Pedersen fidelity | Average gate fidelity formula that is invariant to global phase |
| ARC | Alkali Rydberg Calculator — Python library for atomic physics data |
| QuTiP | Quantum Toolbox in Python — library for quantum dynamics simulation |
| Evered 2023 | State-of-the-art Rb-Rb CZ gate experiment reaching 99.5% fidelity |
| k_eff | Effective wavevector seen by the atom; determines sensitivity to atomic motion |
