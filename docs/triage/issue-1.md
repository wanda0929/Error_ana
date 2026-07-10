# Issue #1: Foundation — Ideal CZ gate simulation + test harness

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: Greenfield implementation — no code exists yet. This is the first issue in the milestone.
- **User impact / severity**: Blocker for all 7 downstream issues (#2–#8). Nothing else can proceed until this lands.
- **Recommended next step**: Implement in the order below. The critical path is params → hamiltonian → protocol → fidelity → tests → figure.
- **Confidence**: high
- **Main blocker**: QuTiP is not installed in the system Python. A conda environment `qutip` exists with QuTiP 5.0.4 + ARC 3.9.0 at `/opt/homebrew/Caskroom/miniconda/base/envs/qutip/`. The project needs either a local venv or to document the conda env.

## 2. Affected-path context

### Affected subsystem

The entire `src/` directory (does not exist yet). This issue creates the simulation foundation that every error-analysis slice (#3–#8) imports.

### Repo vocabulary / components

| Term | Meaning |
|---|---|
| `params.py` | ARC parameter lookup — the only module that talks to the `arc` library |
| `hamiltonian.py` | Builds QuTiP `Qobj` Hamiltonians for 1-atom and 2-atom systems |
| `protocol.py` | Runs the π–2π–π sequential gate as three `sesolve` calls |
| `fidelity.py` | Computes Pedersen average gate fidelity from unitary or Kraus operators |
| Ideal gate | The CZ gate with infinite blockade (U → ∞), zero decay, zero noise |
| Pedersen fidelity | Global-phase-invariant average gate fidelity: F = [Tr(M†M) + |Tr(M)|²] / [d(d+1)] |

### Current relevant control/data flow

None — empty repo. After this issue, the flow is:

```
params.py → hamiltonian.py → protocol.py → fidelity.py
                                  ↓
                          state trajectories → population dynamics figure
```

### Why the key files matter

- `params.py`: every downstream module needs τ, C₆, Γ_e. If ARC returns wrong values, everything is wrong.
- `protocol.py`: the π–2π–π simulation is the core object that error modules modify. Its correctness (F=1 for ideal case) is the foundation test.
- `fidelity.py`: the metric used by every sweep, every test, every figure. Must match the Pedersen formula exactly.

### First file a new implementer should open

`AGENTS.md` (repo-level spec), then start implementation with `src/params.py`.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Build simulation foundation with 7 files covering params, Hamiltonian, protocol, fidelity, tests, and a figure.
- **Reproduction steps**: N/A — greenfield.
- **Expected behavior**: Ideal gate yields F = 1 (infidelity < 10⁻¹⁰). Tests pass. Figure generated.
- **Acceptance criteria** (from issue):
  1. `params.py` returns τ ≈ 374 μs (at T=0K) and C₆ ≈ 863 GHz·μm⁶ for |70,0,0.5⟩
  2. `hamiltonian.py` matrix elements match hand calculation
  3. `protocol.py` ideal run yields 1 − F < 10⁻¹⁰
  4. `fidelity.py` passes known-input tests
  5. `pytest tests/` passes
  6. Population dynamics figure saved
  7. No hard-coded physics

### Repository evidence

**ARC parameter verification** (run during triage on this machine):

```
Rb87 |70, 0, 0.5⟩:
  Radiative lifetime (T=0K):  373.9 μs   ← use this (no BBR)
  Total lifetime (T=300K):    146.5 μs   (with BBR at room temp)
  Decay rate (0K):            2674 Hz
  C₆:                        −862.7 GHz·μm⁶  ✓ (matches AGENTS.md ~860)
  5P₃/₂ linewidth:           2π × 6.07 MHz
```

**Use T=0K (radiative-only lifetime) for `params.py`.** The simulation models Rydberg decay as a separate error channel — BBR-induced transitions are not part of the Lindblad collapse operator √γ|g⟩⟨r|, so the radiative lifetime is the correct input.

**AGENTS.md says τ ≈ 230 μs.** This value does not match n=70 at any temperature. It exactly matches **n=60 at T=0K** (ARC returns 230.0 μs) — likely a copy error from a different parameter set. The correct value for n=70 at T=0K is **τ ≈ 374 μs**. C₆ ≈ 863 GHz·μm⁶ is confirmed for n=70. This does not block the ideal gate (which has no decay).

**Blockade strength at baseline parameters**:

```
R = 6 μm:  U = 18.5 MHz,  U/Ω = 4.6   ← NOT "U >> Ω"
R = 4 μm:  U = 211 MHz,   U/Ω = 53
R = 3 μm:  U = 1183 MHz,  U/Ω = 296
```

The AGENTS.md says R ≈ 6 μm and "U >> Ω", but at R = 6 μm, U/Ω ≈ 4.6 — this is the finite-blockade regime, not the blockade limit. For the ideal simulation (U → ∞) this doesn't matter, but `params.py` should compute and report U/Ω so downstream issues (#4 blockade error) use a physically sensible R. A default of R ≈ 3–4 μm better matches "U >> Ω".

**QuTiP environment verified**:

```
Conda env: /opt/homebrew/Caskroom/miniconda/base/envs/qutip/
QuTiP:     5.0.4
ARC:       3.9.0
Python:    3.9.x
```

QuTiP 5 API confirmed: `sesolve`, `mesolve`, `basis`, `tensor`, `Qobj` all available. Options passed as `dict`, not `Options()` class (v4 → v5 breaking change).

**Ideal gate simulation verified during triage**:

| Approach | Hilbert space | Infidelity | Notes |
|---|---|---|---|
| 4-level, U/Ω = 1000 | {gg,gr,rg,rr} | 3.7 × 10⁻⁷ | Finite blockade leakage |
| 4-level, U/Ω = 10000 | {gg,gr,rg,rr} | 3.7 × 10⁻⁹ | Still above 10⁻¹⁰ |
| 4-level, U/Ω = 100000 | {gg,gr,rg,rr} | 3.7 × 10⁻¹¹ | Passes but numerical issues appearing |
| 4-level, U/Ω = 1000000 | {gg,gr,rg,rr} | **ODE fails** | Stiff system — `IntegratorException` |
| **3-level (U = ∞)** | **{gg,gr,rg}** | **0** (machine precision) | **Recommended for ideal test** |

**Critical finding**: The 3-level model (project out |rr⟩) is the correct approach for the ideal gate. The 4-level model with large U hits stiffness problems before reaching 10⁻¹⁰ reliably.

**Phase structure verified**:

```
π–2π–π protocol produces: diag(1, −1, −1, −1)
This equals CZ × (Z₁ ⊗ Z₂), i.e., CZ up to known local Z corrections.
Entangling phase = π (exact in 3-level model).
Standard in literature — local Z corrections absorbed into gate definition.
```

**Existing tests/fixtures**: None (empty repo).

**Relevant scripts**: Verifier code at `Auto_opt/verifier/` uses similar QuTiP patterns but is NOT a dependency. Reference only.

## 4. Root-cause candidates

N/A — this is a greenfield implementation issue, not a bug. The "root cause" is "no code exists yet." The relevant analysis is about implementation risks.

| Risk | Evidence for | Evidence against / uncertainty | How to mitigate | Severity |
|---|---|---|---|---|
| 4-level ideal model fails the 10⁻¹⁰ threshold | Verified: U/Ω=10000 gives 3.7×10⁻⁹ (fails). U/Ω=1e6 causes ODE crash. | 3-level model gives exactly 0. | Use 3-level model for ideal gate; 4-level for blockade error (issue #4). | high |
| AGENTS.md lifetime (230 μs) is wrong | ARC returns 146.5 μs at 300K. | Could be a different temperature or source. | Trust ARC. Document discrepancy. | medium |
| AGENTS.md interatomic distance (6 μm) gives weak blockade | U/Ω = 4.6 at R=6μm. | Issue #1 doesn't need a physical R (ideal = U=∞). | Note in params.py that R should be chosen to give U/Ω >> 1 for physical simulations. | low (for this issue) |
| QuTiP 5 API differences from v4 | Options class deprecated. Dict required. | sesolve/mesolve core API unchanged. | Use `options={"nsteps": N}` dict syntax. | low |

## 5. Decision

- **Chosen approach**: Implement all 7 files as specified, with one key design decision: **use a 3-level Hilbert space {|gg⟩, |gr⟩, |rg⟩} for the ideal gate** (|rr⟩ projected out = infinite blockade). The 4-level model goes into `hamiltonian.py` as well, ready for issue #4.

- **Recommended implementation**:

  1. `params.py`: Query ARC. Return a dataclass with all physical params. Compute and report U/Ω at the chosen R so the user sees whether they're in the blockade regime.

  2. `hamiltonian.py`: Provide both `build_hamiltonian_ideal(Omega, atom_index)` (3-level, U=∞) and `build_hamiltonian_blockade(Omega, U, atom_index)` (4-level, finite U). The ideal version is used in this issue; the blockade version waits for issue #4.

  3. `protocol.py`: Three sequential `sesolve` calls. For each computational input {00, 01, 10, 11}, return the final state and time-resolved populations. The |00⟩ case is trivial (both dark = identity). The |01⟩ and |10⟩ cases are single-atom evolutions. Only |11⟩ needs the full 2-atom (or 3-level) Hamiltonian.

  4. `fidelity.py`: Implement Pedersen formula. Accept either a unitary matrix (for ideal case) or Kraus operators (for later Lindblad use). Include local-Z-correction logic since the π–2π–π protocol produces CZ × (Z₁⊗Z₂), not bare CZ.

  5. Tests: `test_ideal_gate.py` asserts 1−F < 10⁻¹⁰. `test_fidelity.py` checks known matrix pairs.

  6. Figure: Population dynamics for all 4 basis states.

- **Rejected alternatives**:
  - 4-level model with very large U for ideal test: rejected because it either fails the 10⁻¹⁰ threshold or crashes the ODE solver.
  - Importing from verifier: rejected per AGENTS.md ("raw QuTiP, no verifier dependency").

- **Assumptions**:
  - The conda env `qutip` is the target environment. A `requirements.txt` or env setup instruction should be added.
  - Gate fidelity is computed with optimal local Z corrections (standard practice).

- **Maintainer decision needed**: The AGENTS.md lifetime value (230 μs) contradicts ARC (146.5 μs at 300K). The implementer should use ARC and document the discrepancy, but the project owner should confirm this is acceptable.

## 6. Implementation handoff

### Ordered implementation steps

1. **Environment setup**: Create `requirements.txt` (or document conda env). Verify `import qutip` and `import arc` work.

2. **`src/params.py`** (S-01):
   - Query `arc.Rubidium87().getStateLifetime(70, 0, 0.5, temperature=0, includeLevelsUpTo=100)` (T=0K: radiative-only, no BBR).
   - Query `PairStateInteractions` for C₆ with `(n=70, l=0, j=0.5)` symmetric pair, `m1=m2=0.5`.
   - Return a frozen dataclass: `RydbergParams(tau, gamma, c6, omega, R, U, U_over_Omega, ...)`.
   - Print/log the ARC values on construction so discrepancies are visible.

3. **`src/hamiltonian.py`** (S-02):
   - `build_ideal_hamiltonian(omega, atom_index)`: 3×3 Qobj in {|gg⟩, |gr⟩, |rg⟩} space.
     - Atom 1 pulse: couples |gg⟩ ↔ |rg⟩ only. H = (Ω/2)(|gg⟩⟨rg| + h.c.)
     - Atom 2 pulse: couples |gg⟩ ↔ |gr⟩ only. H = (Ω/2)(|gg⟩⟨gr| + h.c.)
     - The key: |rg⟩ ↔ |rr⟩ and |gr⟩ ↔ |rr⟩ couplings are absent because |rr⟩ is projected out.
   - `build_blockade_hamiltonian(omega, U, atom_index)`: 4×4 Qobj in {|gg⟩, |gr⟩, |rg⟩, |rr⟩}. For issue #4.

4. **`src/fidelity.py`** (S-04, S-05):
   - `pedersen_fidelity(U_actual, U_target)` for unitary matrices.
   - `pedersen_fidelity_kraus(kraus_ops, U_target)` for open-system case.
   - Internal: `optimal_local_z_correction(U_actual)` that finds the Z₁⊗Z₂ correction minimizing distance to CZ.
   - Test vectors: CZ vs CZ → 1.0; CZ vs I → 0.2 (d=4: F = (d + |Tr(I)|²)/(d(d+1)) with Tr(I)=4 gives 1.0 — wait, need to check this).

5. **`src/protocol.py`** (S-03):
   - `run_ideal_gate(omega, n_steps=2000)` → dict with keys `{unitary_diagonal, phases, populations, times}`.
   - For each input |ij⟩:
     - |00⟩: trivial, phase = 0.
     - |01⟩: single-atom 2π evolution → phase.
     - |10⟩: single-atom π + π evolution → phase.
     - |11⟩: 3-level, 3-step sesolve → phase + populations.
   - Record Rydberg population vs time for the population dynamics figure.

6. **`tests/test_fidelity.py`** (T-03):
   - `CZ vs CZ → 1.0`
   - `I vs CZ → known value` (compute by hand: F = (Tr(CZ†) Tr(CZ) + |Tr(CZ†·I)|²)/(d(d+1)) where CZ† = CZ. Tr(CZ) = 2, so |Tr(CZ·I)|² = 4. Tr(CZ†·CZ) = Tr(I) = 4. F = (4+4)/20 = 0.4)
   - Random unitary → bounded in [0, 1]

7. **`tests/test_ideal_gate.py`** (T-01):
   - Run `run_ideal_gate()`, compute Pedersen fidelity against CZ (with Z corrections).
   - Assert `1 - F < 1e-10`.
   - Assert entangling phase = π to within 1e-8.

8. **`scripts/run_ideal.py`** + **`figures/population_dynamics.png`** (F-01):
   - Run the simulation, save population arrays, generate matplotlib figure.
   - 4 curves (one per basis state), 3 pulse regions marked.

### Files likely to change

All new files:
- `src/__init__.py`
- `src/params.py`
- `src/hamiltonian.py`
- `src/protocol.py`
- `src/fidelity.py`
- `tests/test_ideal_gate.py`
- `tests/test_fidelity.py`
- `scripts/run_ideal.py`
- `figures/` directory (gitignored or tracked PNGs)
- `requirements.txt` or environment documentation

### Scope boundaries

- Do **not** implement any error modules (decay, blockade, Doppler, scattering, amplitude). Those are issues #3–#7.
- Do **not** create the HTML explainer. That is issue #2.
- Do **not** implement `sweeps.py` or `analytical.py`. Those come with the error slices.
- The 4-level Hamiltonian builder can be included in `hamiltonian.py` for completeness, but it is not tested or used until issue #4.

### Regression test shape

```python
def test_ideal_gate_fidelity():
    result = run_ideal_gate(omega=2*np.pi*4.0)
    F = pedersen_fidelity(result["unitary_diagonal"], CZ_TARGET)
    assert 1 - F < 1e-10

def test_ideal_entangling_phase():
    result = run_ideal_gate(omega=2*np.pi*4.0)
    ent_phase = result["phases"][3] - result["phases"][1] - result["phases"][2] + result["phases"][0]
    assert abs(abs(ent_phase) - np.pi) < 1e-8
```

### Validation command(s)

```bash
# Activate the conda environment first
conda activate qutip
# Or: source .venv/bin/activate (if a local venv is created)

# Run tests
pytest tests/ -v

# Run ideal simulation and generate figure
python scripts/run_ideal.py
```

These commands are inferred from the project structure; confirm `pytest` is installed in the environment.

### Expected pre-fix signal

No code exists → `pytest` would fail with `ModuleNotFoundError`.

### Expected post-fix signal

```
tests/test_ideal_gate.py::test_ideal_gate_fidelity PASSED
tests/test_ideal_gate.py::test_ideal_entangling_phase PASSED
tests/test_fidelity.py::test_cz_vs_cz PASSED
tests/test_fidelity.py::test_cz_vs_identity PASSED
figures/population_dynamics.png exists and shows 4 curves
```

### Negative/edge validation

- Test that Pedersen fidelity returns 0 ≤ F ≤ 1 for random unitaries.
- Test that `params.py` raises a clear error for invalid quantum numbers (e.g., l ≥ n).
- Test that the ideal gate with Ω → 0 doesn't crash (degenerate case — infinite pulse time).

## 7. Risks, edge cases, rollback

### Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| ARC version mismatch between system (3.10.2) and conda env (3.9.0) gives different τ, C₆ | Low | Medium | Pin ARC version in requirements; document which values to expect. Use T=0K. |
| QuTiP 5 `sesolve` option syntax confuses implementer | Medium | Low | Use `options={"nsteps": 5000}` dict, not `Options()` class |
| 3-level model for ideal gate confuses later blockade work | Low | Low | Document clearly in `hamiltonian.py` that 3-level = ideal, 4-level = finite blockade |
| Population dynamics figure for |11⟩ needs all 3 steps stitched together | Medium | Low | Concatenate time arrays and population arrays from the 3 sesolve calls |

### Edge cases

- **|00⟩ input**: Both qubits dark. No simulation needed — identity. Population dynamics = flat line at 0 Rydberg population.
- **Very small Ω**: Pulse time → ∞. Not a practical concern (Ω is a user parameter with physical defaults), but `sesolve` should not be called with tlist spanning thousands of μs.
- **Pedersen fidelity with non-unitary input**: The Kraus-operator version should handle trace-decreasing maps gracefully (clamp F to [0,1]).

### Rollback/revert notes

This is the first code in the repo. Rollback = `git revert` the implementation commit. No downstream code depends on it yet (issues #2–#8 are unstarted).

## 8. Open questions

### Q1: AGENTS.md lifetime discrepancy (RESOLVED)

- **Question**: AGENTS.md says τ ≈ 230 μs. ARC returns 373.9 μs at T=0K and 146.5 μs at T=300K for Rb87 |70S₁/₂⟩. Neither matches.
- **What was found**: The 230 μs value exactly matches **n=60 at T=0K** (ARC: 230.0 μs). This is a copy error in AGENTS.md — the value was likely carried over from a different parameter set.
- **Resolution**: Use T=0K (radiative-only lifetime). The correct value for n=70 at T=0K is **τ ≈ 374 μs**. Document in `params.py`. Use `temperature=0` in the ARC call.

### Q2: Default interatomic distance

- **Question**: AGENTS.md says R ≈ 6 μm, but U/Ω = 4.6 at that distance — not "U >> Ω".
- **What was found**: R = 3 μm gives U/Ω ≈ 296 (strong blockade). R = 4 μm gives U/Ω ≈ 53 (reasonable).
- **Why it is ambiguous**: The ideal gate (this issue) doesn't use R at all (U = ∞). But `params.py` should pick a sensible default for downstream issues.
- **Options**: (A) Default R = 4 μm (U/Ω ≈ 53). (B) Default R = 3 μm (U/Ω ≈ 296). (C) Default R = 6 μm per AGENTS.md, note that it's weak blockade.
- **Recommended default**: (A) — R = 4 μm is a compromise. Strong enough blockade for the π–2π–π protocol to work, physically realistic tweezer spacing.

### Q3: Environment management

- **Question**: Should the project create a local `.venv` or document the existing conda env?
- **What was found**: Conda env `qutip` at `/opt/homebrew/Caskroom/miniconda/base/envs/qutip/` has QuTiP 5.0.4 + ARC 3.9.0. System Python has ARC 3.10.2 but no QuTiP.
- **Options**: (A) Create `requirements.txt` + local venv. (B) Document conda env. (C) Both.
- **Recommended default**: (A) — a `requirements.txt` makes the project self-contained and reproducible.
