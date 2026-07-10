# Triage Explanation: Issue #1 — Foundation: Ideal CZ Gate Simulation + Test Harness

## TL;DR

Issue #1 is the greenfield "build everything from scratch" ticket for a quantum-computing simulation project. No source code exists yet. The goal is to create the core simulation infrastructure—parameter lookup, Hamiltonian construction, gate protocol execution, fidelity computation, tests, and a figure—for an **ideal controlled-Z (CZ) Rydberg gate**. This is a **blocker** for all seven downstream issues (#2–#8), which add error models, sweeps, and documentation on top of this foundation.

The triage plan prescribes implementing seven new files across `src/`, `tests/`, and `scripts/`, using the existing conda environment (QuTiP 5.0.4 + ARC 3.9.0). The single most important design decision is to use a **3-level Hilbert space** (projecting out the doubly-excited Rydberg state |rr⟩) for the ideal gate simulation, because the alternative 4-level model with large interaction strength either fails the fidelity threshold or crashes the ODE solver. Two parameter discrepancies between the project spec (`AGENTS.md`) and the physics library (ARC) were discovered during triage and flagged for maintainer review.

---

## Decision map

```
Issue #1: "Build the simulation foundation"
│
├─ Q: Does any code exist?
│   └─ No → Greenfield implementation. All files are new.
│
├─ Q: What Hilbert space for the ideal gate?
│   ├─ Option A: 4-level {|gg⟩,|gr⟩,|rg⟩,|rr⟩} with very large U
│   │   └─ REJECTED: U/Ω=10⁴ gives infidelity 3.7×10⁻⁹ (fails 10⁻¹⁰ threshold);
│   │      U/Ω=10⁶ crashes the ODE solver (stiff system)
│   │
│   └─ Option B: 3-level {|gg⟩,|gr⟩,|rg⟩} with |rr⟩ projected out ← CHOSEN
│       └─ Gives infidelity = 0 at machine precision. Passes threshold trivially.
│
├─ Q: AGENTS.md says lifetime τ ≈ 230 μs. ARC says 146.5 μs. Which to trust?
│   └─ Trust ARC (the physics library), document the discrepancy. ← CHOSEN
│       (Neither ARC's 0K value of 374 μs nor its 300K value of 146.5 μs matches 230 μs.)
│
├─ Q: AGENTS.md says R ≈ 6 μm, but U/Ω = 4.6 there (weak blockade). What default?
│   └─ Default to R = 4 μm (U/Ω ≈ 53, strong blockade). ← CHOSEN
│       (Ideal gate uses U = ∞ and doesn't need R, but params.py should set a
│        sensible physical default for downstream error-analysis issues.)
│
└─ Q: Environment — local venv or existing conda env?
    └─ Create a `requirements.txt` for reproducibility. ← CHOSEN
```

---

## State machine trace

This section describes how the triage process moved through its stages. A **triage branch** is a Git branch created to hold the analysis document; it does not contain implementation code. A **triage document** is the Markdown artifact (`docs/triage/issue-1.md`) recording the analysis.

| Stage | What happened |
|---|---|
| **Issue intake** | Issue #1 identified as a greenfield implementation ticket—no bug to reproduce, no existing code to inspect. Classified as *actionable* with *high confidence*. |
| **Context gathering** | Repository found to be empty of source code. `AGENTS.md` reviewed for the project specification. Conda environment located with QuTiP 5.0.4 and ARC 3.9.0. |
| **Parameter verification** | ARC was queried directly for Rb87 |70S₁/₂⟩ lifetime and C₆ coefficient. Two discrepancies found with `AGENTS.md` (lifetime and blockade regime at R=6 μm). |
| **Prototype simulation** | Multiple ideal-gate simulations were run during triage to determine whether the 4-level or 3-level approach meets the 10⁻¹⁰ infidelity threshold. The 3-level model was validated as the correct choice. |
| **Design decision** | 3-level model chosen for ideal gate. 4-level model deferred to issue #4 (blockade error). Implementation order established. |
| **Handoff** | Ordered implementation steps, file list, scope boundaries, test shapes, and validation commands documented. Open questions flagged for maintainer. |

---

## Evidence

The triage document is unusually thorough because the triager ran actual simulations. Here is the key evidence:

### 1. The 4-level ideal model is infeasible

A table of simulation results was produced by sweeping the ratio U/Ω (interaction strength to Rabi frequency) in the 4-level Hilbert space:

| U/Ω | Infidelity | Outcome |
|---|---|---|
| 1,000 | 3.7 × 10⁻⁷ | Fails threshold |
| 10,000 | 3.7 × 10⁻⁹ | Fails threshold (barely) |
| 100,000 | 3.7 × 10⁻¹¹ | Passes, but numerical stability degrades |
| 1,000,000 | — | **ODE solver crash** (`IntegratorException`) |
| **3-level (U=∞)** | **0** | **Machine precision, no solver issues** |

This is the strongest evidence in the document. It justifies the central design decision and was obtained empirically, not assumed.

### 2. ARC parameter discrepancies

The triager queried ARC for the Rydberg state |70, 0, 0.5⟩ of Rubidium-87:

- **Lifetime**: ARC returns 146.5 μs at 300K. `AGENTS.md` claims ~230 μs. Neither ARC's 0K radiative lifetime (373.9 μs) nor its 300K value matches. **Unresolved.**
- **C₆ coefficient**: ARC returns −862.7 GHz·μm⁶. `AGENTS.md` claims ~860 GHz·μm⁶. **Match** (sign aside).
- **Blockade at R=6 μm**: U/Ω = 4.6. This is *not* the strong-blockade regime ("U >> Ω") that `AGENTS.md` claims. R ≈ 3–4 μm is needed for strong blockade.

### 3. Phase structure of the π–2π–π protocol

The protocol (three sequential laser pulses: π on atom 1, 2π on atom 2, π on atom 1) produces the diagonal unitary `diag(1, −1, −1, −1)`. This equals `CZ × (Z₁ ⊗ Z₂)`—the target CZ gate up to known single-qubit Z rotations. This is standard in the Rydberg gate literature and means the fidelity module must include a local-Z-correction step.

### 4. QuTiP 5 API confirmed

QuTiP 5.0.4 in the conda env was verified to support all needed functions (`sesolve`, `mesolve`, `basis`, `tensor`, `Qobj`). One API break from v4 was noted: solver options must be passed as a plain `dict`, not an `Options()` object.

---

## Implementation plan

The plan creates the simulation from nothing. The terms below are project-specific:

- **`params.py`**: The module that queries the ARC (Alkali Rydberg Calculator) library for physical constants of the Rydberg atom.
- **`hamiltonian.py`**: Constructs the quantum Hamiltonian matrices (as QuTiP `Qobj` objects) that describe the laser-atom interaction.
- **`protocol.py`**: Executes the π–2π–π pulse sequence by calling QuTiP's Schrödinger equation solver (`sesolve`) three times in succession.
- **`fidelity.py`**: Computes the Pedersen average gate fidelity, a standard metric that is insensitive to global phase.
- **Pedersen fidelity**: F = [Tr(M†M) + |Tr(M)|²] / [d(d+1)], where M = U_target† · U_actual and d is the Hilbert space dimension.

### Ordered steps

| Step | File(s) | What it does | Why this order |
|---|---|---|---|
| 0 | `requirements.txt` | Pin QuTiP ≥ 5.0, ARC, numpy, matplotlib, pytest | Everything else needs imports to work |
| 1 | `src/params.py` | Query ARC for τ, C₆, Γ. Compute Ω, U, U/Ω. Return a frozen dataclass. | Every other module imports physical parameters from here |
| 2 | `src/hamiltonian.py` | Build 3-level (ideal, U=∞) and 4-level (finite blockade) Hamiltonians as QuTiP Qobj matrices | Protocol needs Hamiltonians to run |
| 3 | `src/fidelity.py` | Implement Pedersen formula for unitary and Kraus-operator inputs. Include local-Z-correction optimizer. | Needed to score the protocol output |
| 4 | `src/protocol.py` | Run the π–2π–π gate: three `sesolve` calls per input state. Return unitary diagonal, phases, and time-resolved populations. | Depends on steps 1–3 |
| 5 | `tests/test_fidelity.py` | Unit tests: CZ vs CZ → 1.0; CZ vs I → 0.4; random unitary → F ∈ [0,1] | Validates the metric before it's used to score the gate |
| 6 | `tests/test_ideal_gate.py` | Integration test: run ideal gate, assert 1−F < 10⁻¹⁰, assert entangling phase = π | The acceptance criterion |
| 7 | `scripts/run_ideal.py` + `figures/` | Run simulation, plot population dynamics for all 4 basis states, save PNG | Visual validation and deliverable figure |

### Critical path

```
params.py → hamiltonian.py → protocol.py → fidelity.py → tests → figure
```

Each module depends on the one before it. There is no parallelism in this chain.

### Scope boundaries (what NOT to build)

- No error/noise models (those are issues #3–#7: spontaneous decay, finite blockade, Doppler shifts, Raman scattering, amplitude noise).
- No HTML explainer (issue #2).
- No sweep or analytical comparison modules.
- The 4-level Hamiltonian builder can be *included* in `hamiltonian.py` for convenience, but it is not *tested or used* until issue #4.

---

## Validation

### Primary acceptance test

```python
def test_ideal_gate_fidelity():
    result = run_ideal_gate(omega=2*np.pi*4.0)  # 4 MHz Rabi frequency
    F = pedersen_fidelity(result["unitary_diagonal"], CZ_TARGET)
    assert 1 - F < 1e-10
```

This directly maps to acceptance criterion #3 from the issue: "protocol.py ideal run yields 1 − F < 10⁻¹⁰."

### Secondary tests

| Test | What it checks |
|---|---|
| `test_ideal_entangling_phase` | The accumulated phase difference φ₁₁ − φ₀₁ − φ₁₀ + φ₀₀ = π (to within 10⁻⁸) |
| `test_cz_vs_cz` | Pedersen fidelity of CZ against itself = 1.0 exactly |
| `test_cz_vs_identity` | Pedersen fidelity of I against CZ = 0.4 (hand-computed) |
| `test_random_unitary_bounded` | F ∈ [0, 1] for randomly generated unitaries |
| `test_params_invalid_quantum_numbers` | `params.py` raises a clear error for unphysical inputs (e.g., l ≥ n) |

### Validation commands

```bash
conda activate qutip
pytest tests/ -v
python scripts/run_ideal.py
# → figures/population_dynamics.png should exist and show 4 curves with 3 pulse regions
```

### Expected signals

- **Before implementation**: `pytest` fails with `ModuleNotFoundError` (no `src/` directory).
- **After implementation**: All tests pass. Figure file exists. Infidelity reported as < 10⁻¹⁰.

---

## Caveats and open questions

These are preserved verbatim from the triage document. None are resolved; all require human judgment.

### Q1: Lifetime discrepancy (HIGH priority for downstream issues)

`AGENTS.md` specifies τ ≈ 230 μs for the Rydberg state. ARC returns **146.5 μs at 300K** and **373.9 μs at 0K**. Neither matches. The triage recommends trusting ARC (per `AGENTS.md` standing order #5: "use ARC as the authoritative source"), but this affects every error estimate in issues #3–#8. **The project owner should confirm this is acceptable.**

### Q2: Default interatomic distance

`AGENTS.md` specifies R ≈ 6 μm, but at that distance the blockade ratio U/Ω = 4.6—far from the "U >> Ω" regime the spec claims. The ideal gate (this issue) uses U = ∞ and doesn't need R, but `params.py` must pick a sensible default for downstream use. Triage recommends R = 4 μm (U/Ω ≈ 53). **This is a physics judgment call, not a code decision.**

### Q3: Environment management

A conda environment with the right packages exists, but the project has no `requirements.txt` or setup instructions. Triage recommends creating a `requirements.txt` for reproducibility. This is low-risk but needs a decision on whether to also document the conda path.

### Caveat: 3-level vs. 4-level model confusion risk

The 3-level model (U = ∞, |rr⟩ projected out) is correct for the ideal gate but is a *different* Hilbert space than the 4-level model used by error-analysis issues. Implementers of issue #4 (blockade error) will need to switch to the 4-level model. The triage plan mitigates this by including both builders in `hamiltonian.py` with clear naming and documentation, but the conceptual discontinuity remains.

### Caveat: Phase convention

The π–2π–π protocol produces `diag(1, −1, −1, −1)` = CZ × (Z₁⊗Z₂), not bare CZ. The fidelity module must apply local Z corrections before scoring. If this correction logic is wrong, the fidelity will be reported as < 1 even for a perfect gate. The triage document flags this but does not provide a reference implementation of the correction optimizer.

---

## Revision targets

These are the specific decisions in the triage plan that a human annotator might reasonably want to change, override, or request more information about:

| # | Decision | What to review | Risk if wrong |
|---|---|---|---|
| **R1** | Use 3-level Hilbert space for ideal gate | Is projecting out \|rr⟩ consistent with the project's pedagogical goals? Issue #1 may want to show all 4 levels even if fidelity is slightly imperfect. | Low — only affects ideal test threshold |
| **R2** | Trust ARC lifetime (146.5 μs) over AGENTS.md (230 μs) | The 230 μs may come from a specific paper the PI wants to match. Changing it affects all decay-related error estimates in #3–#8. | **High** — propagates through every error module |
| **R3** | Default R = 4 μm instead of AGENTS.md's 6 μm | The PI may have chosen 6 μm deliberately (e.g., to match a specific experimental setup). | Medium — affects blockade error (#4) quantitatively |
| **R4** | Include 4-level Hamiltonian in `hamiltonian.py` now (but don't test it) | Could introduce dead code that drifts. Alternative: defer entirely to issue #4. | Low |
| **R5** | Create `requirements.txt` vs. document conda env | If other collaborators use a different OS or package manager, a requirements.txt alone may not suffice. | Low |
| **R6** | Pedersen fidelity with local Z correction as the metric | Some references use bare process fidelity or entanglement fidelity instead. Confirm the metric matches the project's target publication/comparison. | Medium — wrong metric = wrong conclusions |
| **R7** | QuTiP solver options (`nsteps`, tolerances) | The triage doesn't prescribe specific solver tolerances. For the 3-level ideal model this is irrelevant (exact result), but downstream issues with stiff Lindblad equations may need careful tuning established here. | Low for this issue, medium for later issues |