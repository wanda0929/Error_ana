# AGENTS.md — Rydberg CZ Gate Error Analysis for Rb-Rb

## What this project does

Simulate a **two-qubit CZ gate** between two ⁸⁷Rb atoms using the **Rydberg blockade effect**, then systematically study how each experimental error source degrades gate fidelity.

The project produces:
1. **Python simulation code** (raw QuTiP, no external gate libraries)
2. **A self-contained interactive HTML explainer** with embedded figures

The HTML is the primary deliverable. A physicist with no prior Rydberg knowledge should be able to read it and understand: what the gate does, how each error source works, and how badly each one hurts.

---

## Physics summary

### System

Two ⁸⁷Rb atoms in optical tweezers, separated by distance R.

Each atom has two relevant levels:
- **|g⟩**: ground state (qubit state |1⟩; state |0⟩ is dark and never couples to the laser)
- **|r⟩**: Rydberg state |70S₁/₂⟩

Excitation scheme: **two-photon** (780 nm + 480 nm), matching the Evered et al. 2023 experimental platform.

### Baseline parameters (from Evered 2023)

| Parameter | Symbol | Value | Source |
|---|---|---|---|
| Rydberg state | \|n,l,j⟩ | \|70, 0, 0.5⟩ | Evered 2023 |
| Rabi frequency | Ω/2π | ~4 MHz | Evered 2023 |
| Rydberg lifetime (300K) | τ | ~230 μs | ARC (compute at runtime) |
| C₆ coefficient | C₆ | ~860 GHz·μm⁶ | ARC (compute at runtime) |
| Interatomic distance | R | ~6 μm | chosen so U >> Ω |
| Blockade strength | U = C₆/R⁶ | >> Ω | blockade regime |
| Atom temperature | T | 10 μK | Evered 2023 |
| Intermediate detuning | Δ_p/2π | ~1 GHz | Evered 2023 |
| Amplitude noise | σ_Ω | 2% | Evered 2023 |
| Effective wavevector | k_eff | \|1/λ₇₈₀ − 1/λ₄₈₀\| ≈ 0.80 μm⁻¹ | counter-propagating two-photon |

**Key regime**: U/Ω >> 1 (blockade regime). The doubly-excited state |rr⟩ is energetically forbidden during the gate, which is what makes the CZ gate work.

Exact values for τ and C₆ should be computed from ARC (`arc.Rubidium87`) at runtime, not hard-coded. The table above gives approximate values for orientation.

### Gate protocol: π–2π–π sequential

The textbook Rydberg blockade CZ gate (Jaksch 2000, Saffman 2010):

```
Step 1:  π-pulse on atom 1     |g⟩₁ → |r⟩₁
Step 2:  2π-pulse on atom 2    |g⟩₂ → |r⟩₂ → |g⟩₂  (if atom 1 not in |r⟩)
                                blocked              (if atom 1 in |r⟩)
Step 3:  π-pulse on atom 1     |r⟩₁ → |g⟩₁
```

Effect on computational basis states:
- |00⟩ → |00⟩ (both dark, no dynamics, phase 0)
- |01⟩ → |01⟩ (atom 2 does 2π rotation, picks up phase −1)
- |10⟩ → |10⟩ (atom 1 does π–π round-trip, atom 2 dark, picks up phase −1)
- |11⟩ → |11⟩ (atom 1 does π–π, atom 2 blocked → no 2π phase, net phase +1 relative)

Result: CZ = diag(1, 1, 1, −1) up to single-qubit phases.

Total Rydberg occupation time: T_R = 7π/(4Ω) (Zhang 2012).

### Target unitary

```
CZ = diag(1, 1, 1, -1)
```

Single-qubit Z phases may appear; the fidelity metric (Pedersen average gate fidelity) is invariant to global phase. Local phase corrections are reported as diagnostics.

---

## Approach: Analytical + Numerical (Approach C)

For each error source, the project provides **both**:

1. **Analytical formula**: closed-form scaling law from perturbation theory (Zhang 2012, Saffman 2010, de Léséleuc 2018). Gives physical insight — you see *why* the error scales the way it does.

2. **Numerical simulation**: full QuTiP master-equation or Monte Carlo simulation with that error source turned on (all others ideal). Gives the exact answer within the model.

Each error source gets a plot with the **analytical curve overlaid on numerical data points**. Agreement validates the analytical formula; disagreement reveals where perturbation theory breaks down.

A vertical dashed line marks the **Evered 2023 operating point** on every plot, so the reader sees where real experiments sit.

---

## Error sources

Five error channels, studied individually then combined.

### 1. Rydberg decay

**Physics**: The Rydberg state |r⟩ has a finite lifetime τ. While the atom is in |r⟩ during the gate, it can spontaneously decay back to |g⟩ (or to other states), destroying coherence.

**Analytical formula**:
```
ε_decay = Γ_eff × T_R = (1/τ) × 7π/(4Ω)
```
Scales as 1/(Ω·τ): faster gates or longer-lived states → less decay error.

**Numerical model**: Lindblad master equation with collapse operator L = √γ |g⟩⟨r| where γ = 1/τ.

**Sweep parameter**: Vary decay rate γ (equivalently, lifetime τ).

### 2. Finite blockade

**Physics**: The blockade is not perfect. When atom 1 is in |r⟩, the 2π-pulse on atom 2 is detuned by U (the interaction strength), not infinitely shifted. There is a small probability of doubly exciting to |rr⟩.

**Analytical formula**:
```
ε_block = Ω² / (8U²)
```
Scales as (Ω/U)²: stronger blockade → less error.

**Numerical model**: Include the |rr⟩ state in the Hamiltonian with energy shift U. No Lindblad needed — this is a coherent error.

**Sweep parameter**: Vary the ratio U/Ω.

### 3. Doppler dephasing

**Physics**: Atoms in the trap have thermal motion. An atom moving with velocity v sees the laser frequency shifted by k_eff × v. This effective detuning causes phase errors. For counter-propagating two-photon Rb excitation, k_eff = |1/λ₇₈₀ − 1/λ₄₈₀| ≈ 0.80 μm⁻¹ (small — this is an advantage of two-photon schemes).

**Analytical formula**:
```
ε_Doppler ∝ k_eff² × k_B × T_atom × t_gate² / m
```
Scales with temperature and gate time squared.

**Numerical model**: Monte Carlo — sample atomic velocities from the thermal (Maxwell-Boltzmann) distribution, add detuning δ = k_eff × v to each shot, average fidelity over many samples.

**Sweep parameter**: Vary atom temperature T.

### 4. Intermediate-state scattering

**Physics**: The two-photon Rydberg excitation passes through an intermediate state |e⟩ (5P₃/₂ for the 780 nm leg). If the intermediate-state detuning Δ_p is not large enough, the atom can scatter a photon from |e⟩, decohering the gate.

**Analytical formula**:
```
ε_scatter ≈ (Γ_e / Δ_p) × (Ω_eff / Ω₁) × (π / Ω_eff) × (number of π-pulses)
```
where Γ_e is the linewidth of |e⟩, Ω₁ is the single-photon Rabi frequency on the lower leg, and Ω_eff is the effective two-photon Rabi frequency. Scales as 1/Δ_p: larger intermediate detuning → less scattering. See de Léséleuc 2018 for the precise form.

**Numerical model**: Add an effective dephasing Lindblad channel with rate γ_scatter = Γ_e × (Ω₁/(2Δ_p))² during the pulse, representing photon scattering from the intermediate state.

**Sweep parameter**: Vary intermediate-state detuning Δ_p.

### 5. Amplitude noise

**Physics**: Shot-to-shot fluctuations in laser intensity cause the Rabi frequency to vary: Ω → Ω(1 + ε) where ε ~ N(0, σ²). A π-pulse becomes a (π + δ)-pulse.

**Analytical formula**:
```
ε_amp ∝ σ_Ω² × S_Ω
```
where S_Ω is a pulse-shape-dependent sensitivity factor. For the π–2π–π protocol, this can be computed analytically. Scales as σ² — quadratic in noise strength.

**Numerical model**: Monte Carlo — sample ε from a Gaussian, scale Ω → Ω(1+ε) for the entire pulse, run coherent simulation, average fidelity over many samples.

**Sweep parameter**: Vary σ_Ω (fractional Rabi frequency noise).

### 6. Combined (all errors on)

Turn on all five error sources simultaneously at the Evered 2023 operating point. Compare the numerical fidelity against the additive approximation:

```
ε_total ≈ ε_decay + ε_block + ε_Doppler + ε_scatter + ε_amp
```

If they match, the perturbative error budget is valid. If they disagree, the discrepancy reveals cross-talk between error channels — a genuine physics result worth discussing.

---

## Fidelity metric

**Pedersen average gate fidelity** (Pedersen et al. 2007):

```
F_avg = [Tr(M†M) + |Tr(M)|²] / [d(d+1)]
```

where M = U_target† × U_actual and d = 4 (two-qubit Hilbert space dimension).

For simulations with Lindblad dynamics (decay, scattering), extract Kraus operators {G_k} from the final density matrices and use:

```
F_avg = [Σ_k Tr(M_k†M_k) + Σ_k |Tr(M_k)|²] / [d(d+1)],   M_k = U_target† G_k
```

This metric is:
- Invariant to global phase (good — CZ is defined up to global phase)
- Sensitive to both population and phase errors
- Standard in the field (used by Evered 2023, Levine 2019, etc.)

---

## Simulation framework

**Raw QuTiP** — no Pulser, no verifier dependency, no external gate libraries.

### Hilbert space

- **Ideal / blockade-regime**: Each atom has 2 levels {|g⟩, |r⟩}. Two atoms → 4-dimensional space: {|gg⟩, |gr⟩, |rg⟩, |rr⟩}.
- **Finite blockade study**: Same 4-level space, but |rr⟩ is *not* projected out — it sits at energy U above |gr⟩ and |rg⟩.

### Hamiltonian

For the two-atom system during a pulse on atom i with Rabi frequency Ω_i and detuning Δ_i:

```
H = (Ω_i/2)(|g⟩⟨r|_i + h.c.) − Δ_i |r⟩⟨r|_i + U |rr⟩⟨rr|
```

The π–2π–π protocol applies pulses sequentially to atom 1, then atom 2, then atom 1 again. Each step is a separate time evolution.

### Code requirements

- Use `qutip.mesolve` for Lindblad dynamics (decay, scattering errors).
- Use `qutip.sesolve` for coherent-only simulations (ideal, finite blockade, amplitude noise Monte Carlo).
- All physical constants from ARC (`arc.Rubidium87`) — no hard-coded lifetimes or C₆ values.
- Fidelity computation: implement Pedersen formula directly from the 4×4 process matrix. Do not use external fidelity libraries.

---

## Repo structure

```
Error_ana/
├── AGENTS.md              # This file
├── TASK.md                # Current task for the agent
│
├── src/                   # Python source code
│   ├── params.py          # Physical parameters from ARC + Evered 2023 defaults
│   ├── hamiltonian.py     # Build H for 1-atom and 2-atom systems
│   ├── protocol.py        # π–2π–π sequential gate: ideal + with errors
│   ├── errors/            # One module per error source
│   │   ├── __init__.py
│   │   ├── decay.py       # Rydberg decay (Lindblad)
│   │   ├── blockade.py    # Finite blockade (coherent)
│   │   ├── doppler.py     # Doppler dephasing (Monte Carlo)
│   │   ├── scattering.py  # Intermediate-state scattering (Lindblad)
│   │   └── amplitude.py   # Amplitude noise (Monte Carlo)
│   ├── fidelity.py        # Pedersen average gate fidelity
│   ├── analytical.py      # Closed-form error scaling laws
│   └── sweeps.py          # Parameter sweep runners (vary error param → fidelity)
│
├── scripts/
│   ├── run_ideal.py       # Perfect-condition simulation, verify CZ
│   ├── run_sweeps.py      # Run all parameter sweeps, save data
│   └── generate_figures.py # Produce all plots as PNGs
│
├── figures/               # Generated PNGs (gitignored, regenerated by scripts)
│
├── site/
│   └── index.html         # Self-contained interactive HTML explainer
│
└── tests/
    ├── test_ideal_gate.py     # Ideal CZ fidelity == 1 (to numerical precision)
    ├── test_analytical.py     # Analytical formulas match known literature values
    └── test_fidelity.py       # Pedersen fidelity: known inputs → known outputs
```

---

## HTML explainer requirements

The HTML file (`site/index.html`) is the primary deliverable. Requirements:

### Content structure

1. **Introduction**: What is a Rydberg CZ gate? Two Rb atoms, ground ↔ Rydberg coupling, blockade effect. Define every term before using it.

2. **Ideal gate**: The π–2π–π protocol explained step by step. Show the truth table. Show numerical simulation confirming F = 1 (or very close). Include a population dynamics plot (Rydberg population vs. time for each basis state).

3. **Error sources** (one section each):
   - Plain-language explanation of the physics
   - The analytical formula with a derivation sketch (not just "from Zhang 2012" — show the key step)
   - The plot: fidelity vs. error parameter, with analytical curve + numerical dots + Evered 2023 reference line
   - One-sentence takeaway: "This error dominates when ___"

4. **Combined error budget**: All errors on. Bar chart of individual contributions vs. the additive sum vs. the full numerical result. Discuss whether the additive approximation holds.

5. **Summary table**: Error budget at the Evered 2023 operating point. One row per error source, columns for analytical estimate, numerical result, and dominant scaling.

### Style

- Self-contained single HTML file. All CSS/JS inline. No external dependencies except MathJax CDN for equations.
- Clean, readable typography. No flashy frameworks.
- Figures embedded as base64 PNGs or generated inline.
- MathJax for all equations.
- Collapsible derivation sections (click to expand) so the main flow stays clean.
- Mobile-friendly (responsive layout).

### Audience

A physicist who knows quantum mechanics but has never worked with Rydberg atoms. They can follow a Hamiltonian and a master equation. They cannot read Python. Do not reference code — explain the physics.

---

## Standing orders

1. **Read `TASK.md` first** — it contains your current assignment. Do that task. Don't freelance.
2. **Ideal gate first** — before adding any error source, the ideal (zero-error) simulation must produce F = 1 to numerical precision (< 10⁻¹⁰ infidelity). This is the sanity check that everything downstream depends on.
3. **One error at a time** — implement and validate each error source in isolation before combining. Never debug two new error channels simultaneously.
4. **Analytical before numerical** — for each error source, implement the analytical formula first, compute the expected value at the Evered 2023 operating point, then check that the numerical simulation agrees. If they disagree by more than a factor of 2, something is wrong — investigate before moving on.
5. **ARC for physics, not hard-coded numbers** — lifetimes, C₆, transition frequencies come from `arc.Rubidium87` at runtime. The only hard-coded values are the quantum numbers (n=70, l=0, j=0.5) and experimental parameters (T=10μK, σ_Ω=0.02, etc.).
6. **Test at boundaries** — every error module must pass: (a) zero error parameter → F = 1, (b) large error parameter → F drops significantly, (c) intermediate value matches analytical formula.
7. **Figures are publication-quality** — axis labels with units, legends, consistent color scheme, no default matplotlib styling. The HTML reader will judge the project by the figures.
8. **When writing the HTML**, explain physics, not code. The reader never sees `src/`. They see equations, figures, and plain English.

---

## Key references

| Short name | Full reference | Used for |
|---|---|---|
| Evered 2023 | Evered et al., Nature 622, 268 (2023) | Baseline experimental parameters, error budget template |
| Saffman 2010 | Saffman, Walker & Mølmer, Rev. Mod. Phys. 82, 2313 (2010) | Blockade gate theory, error scaling |
| Zhang 2012 | Zhang, Robicheaux & Saffman, PRA 84, 043408 (2011) | π–2π–π protocol, Rydberg occupation time, decay error formula |
| de Léséleuc 2018 | de Léséleuc et al., PRA 97, 053803 (2018) | Intermediate-state scattering formula |
| Levine 2019 | Levine et al., PRL 123, 170503 (2019) | Experimental CZ gate, phase calibration |
| Pedersen 2007 | Pedersen et al., Phys. Lett. A 367, 47 (2007) | Average gate fidelity formula |
| Jaksch 2000 | Jaksch et al., PRL 85, 2208 (2000) | Original Rydberg blockade gate proposal |
