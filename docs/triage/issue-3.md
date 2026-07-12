# Issue #3: Rydberg decay: simulation, analytical formula, sweep, figure, HTML section

## 1. Executive summary

- Actionability: actionable
- Symptom: the repository has the ideal CZ foundation and HTML scaffold, but no Rydberg-decay error slice yet: no decay module, no analytical formula module, no sweep runner, no decay figure, and no HTML decay section.
- User impact / severity: high for the milestone. Issue #3 is the first error-channel vertical slice and establishes conventions that later error slices will copy.
- Recommended next step: implement the decay slice narrowly, with analytical formula first, a Lindblad `mesolve` channel simulation second, then the sweep/figure/HTML/test layer. Do not infer average gate fidelity from only four basis-state populations; build or extract a channel/Kraus representation so the Pedersen open-system fidelity is meaningful.
- Confidence: medium
- Main blocker, if any: no external blocker remains in `origin/main`; issues #1 and #2 appear merged. The main open decision is which lifetime to label as the “Evered operating point”: the current code uses ARC radiative-only `|70S_1/2>` lifetime, about 373.9 µs, while the issue text asks for about 230 µs.

## 2. Affected-path context

- Affected subsystem: the first non-ideal error-channel vertical slice: `src/errors/decay.py`, `src/analytical.py`, `src/sweeps.py`, plotting scripts, decay tests, and the error-source section in `site/index.html`.
- Repo vocabulary / components needed to understand this issue:
  - `|0⟩`: dark qubit state. It does not couple to the Rydberg laser.
  - `|g⟩`: qubit `|1⟩`, the ground state that couples to `|r⟩`.
  - `|r⟩`: auxiliary Rydberg state used during the gate and returned to `|g⟩` by the end of an ideal pulse sequence.
  - `Ω`: angular Rabi frequency. Current repo units are rad/µs for `Ω`, decay rates, and Hamiltonians.
  - `γ = 1/τ`: Rydberg decay rate. In this issue the requested collapse operator is `L = sqrt(γ) |g⟩⟨r|` for each atom.
  - `π–2π–π`: the gate sequence implemented in `src/protocol.py`: π pulse on atom 1, 2π pulse on atom 2, π pulse on atom 1.
  - Pedersen fidelity: the project’s average gate fidelity metric. For an open-system map this requires Kraus operators or an equivalent channel representation, not just final populations.
- Current relevant control/data flow:

  ```text
  src/params.py
      -> baseline Ω, τ, γ
      -> src/analytical.py::epsilon_decay(omega, tau)
      -> src/errors/decay.py::run_decay_gate(...)
      -> src/sweeps.py decay sweep
      -> scripts/run_sweeps.py writes tabular sweep data
      -> scripts/generate_figures.py writes figures/fidelity_vs_decay_rate.png
      -> site/index.html embeds the figure and explanation
  ```

- Why the key files matter:
  - `src/params.py`: establishes units and baseline parameters. It currently computes ARC values at runtime and returns `rydberg_lifetime_us`, `rydberg_decay_rate_per_us`, and `omega_rad_per_us`.
  - `src/protocol.py`: establishes pulse durations and ideal-gate timing. Decay should use the same `t_pi = π/Ω`, total time `4π/Ω`, and local-Z phase convention.
  - `src/fidelity.py`: already contains the Pedersen unitary and Kraus formulas, plus local-Z correction utilities for the raw blockade gate phase.
  - `site/index.html`: already contains an error-source section template for issues #3–#7. The decay section should plug into that shape instead of inventing a new layout.
  - `scripts/run_ideal.py`: establishes plotting style and project-root import conventions for scripts.
- First file a new implementer should open: `src/protocol.py`, then `src/fidelity.py`, then the error-source template near the end of `site/index.html`.

## 3. Observed facts

### Issue evidence

- Reporter claim: implement an end-to-end vertical slice for the Rydberg decay channel: Lindblad simulation, analytical formula, sweep, figure, HTML section, and tests.
- Reproduction steps, if provided: none; this is a feature slice, not a reported runtime bug.
- Expected behavior:
  - `γ = 0` gives ideal-gate fidelity with `1 - F < 1e-10`.
  - At the baseline/Evered decay rate, numerical decay error matches `ε_decay = Γ × 7π/(4Ω)` within 50%.
  - Sweep has at least 20 points over at least 2 orders of magnitude in `γ`.
  - Figure is saved at `figures/fidelity_vs_decay_rate.png` with numerical points, analytical curve, and a baseline vertical line.
  - HTML gains a “Rydberg Decay” section with explanation, formula, derivation sketch, figure, and takeaway.
  - `pytest tests/test_decay.py` passes.
- Actual behavior: as of the inspected branch, the decay slice does not exist.
- Logs/errors/screenshots, if provided: none in the issue payload.
- Environment/version/config details, if provided: none in the issue payload, but the repository pins `qutip==5.0.4` and `ARC-Alkali-Rydberg-Calculator==3.9.0` in `requirements.txt`.
- Reproducibility rate, if known: N/A.
- Relevant comments or corrections: no issue comments were present in the provided payload.

### Repository evidence

For every code claim, include file references when possible.

- File references and what each proves:
  - `src/params.py`: current baseline constants are `DEFAULT_OMEGA_MHZ = 4.0`, `DEFAULT_OMEGA_RAD_PER_US = 2π×4`, and `DEFAULT_DISTANCE_UM = 3.0`. `get_rydberg_params()` queries ARC and returns radiative-only Rydberg lifetime at `temperature=0`. On the verified local conda environment it reports `τ ≈ 373.9 µs`, `γ ≈ 0.0026747 µs⁻¹`, and analytical `ε_decay ≈ 5.85e-4` at `Ω/2π = 4 MHz`.
  - `src/hamiltonian.py`: contains coherent Hamiltonian builders only: `single_atom_hamiltonian`, `build_ideal_hamiltonian`, and `build_blockade_hamiltonian`. There are no Lindblad collapse operators and no dark-state two-atom open-system basis.
  - `src/protocol.py`: implements the ideal π–2π–π sequence with `qutip.sesolve`. It is branch-specific for the four computational inputs, which is fine for ideal population plots but insufficient by itself for an open-system average gate fidelity on arbitrary superpositions.
  - `src/fidelity.py`: implements `pedersen_fidelity()` and `pedersen_fidelity_kraus()`. It also defines `LOCAL_Z_PRODUCT` and `apply_local_z_correction()`; this matters because the raw ideal blockade sequence is `diag(1, -1, -1, -1)`, not canonical `CZ = diag(1, 1, 1, -1)`, until local Z phases are corrected.
  - `scripts/run_ideal.py`: shows the script convention: add repository root to `sys.path`, use Matplotlib Agg backend, write generated PNGs under `figures/`, and print a short numerical report.
  - `site/index.html`: exists and includes the introduction and ideal-gate sections. Near the end it contains `ERROR-SOURCE SECTION TEMPLATE FOR ISSUES #3-#7`, including slots for physics text, analytical formula, collapsible derivation, figure, and takeaway.
  - `.gitignore`: ignores `figures/*.png`. Generated decay figures are expected artifacts but are not automatically tracked by git unless the ignore policy changes.
  - Missing files/directories: `src/errors/`, `src/analytical.py`, `src/sweeps.py`, `scripts/run_sweeps.py`, `scripts/generate_figures.py`, and `tests/test_decay.py` are absent.
- Existing tests/fixtures covering nearby behavior:
  - `tests/test_ideal_gate.py`: verifies ideal fidelity with local-Z correction and validates population traces.
  - `tests/test_fidelity.py`: verifies CZ fidelity, raw local-Z phase correction, Kraus formula for a unitary, and bounded fidelity for random unitaries.
  - `tests/test_hamiltonian.py`: verifies coherent Hamiltonian matrix elements.
  - `tests/test_params.py`: verifies ARC-derived lifetime/C6 ranges for the current parameter convention.
  - `tests/test_html_scaffold.py`: verifies the HTML scaffold, MathJax, ideal figure reference, and downstream error-section template.
- Relevant scripts or validation commands found:
  - Verified command in the project conda environment:

    ```bash
    PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q
    ```

    Result during triage: `21 passed in 1.55s`.

  - System Python does not have QuTiP installed. Running plain `pytest -q` from the default shell failed during collection, including `ModuleNotFoundError: No module named 'qutip'` when importing directly. Use the pinned environment or install `requirements.txt` before validating implementation.
- Important evidence not found:
  - No existing open-system protocol function.
  - No helper to reconstruct Kraus operators from a simulated process matrix.
  - No persisted sweep-data format yet.
  - No project convention yet for whether generated figure PNGs should stay ignored or be force-added for publication artifacts.
  - No maintainer comment resolving the `τ≈230 µs` issue wording versus current `src/params.py` `τ≈373.9 µs` convention.

## 4. Root-cause candidates

Include 2–4 plausible candidates, or 1 if evidence is conclusive. Do not invent evidence.

| Candidate cause | Evidence for | Evidence against / uncertainty | How to falsify | Confidence |
| --- | --- | --- | --- | --- |
| The decay slice is simply unimplemented after the foundation/scaffold work landed. | `src/errors/`, `src/analytical.py`, `src/sweeps.py`, decay scripts, and `tests/test_decay.py` are absent. Issue #3 is explicitly a vertical-slice implementation request. | None. This is the direct state of the repo. | Create a branch with the requested files and show `tests/test_decay.py`, sweep generation, figure generation, and HTML scaffold checks passing. | high |
| The numerical-fidelity part can be implemented incorrectly if it uses only four final density matrices. | Issue asks for density matrices for four basis inputs, but Pedersen average gate fidelity for a Lindblad channel depends on coherences/superpositions. `src/fidelity.py` has a Kraus formula, not a “basis populations only” fidelity. | If the acceptance test only checks basis-state survival, a weaker implementation might appear to pass while giving the wrong process fidelity. | Add a regression test at `γ=0` and for a coherent superposition/process reconstruction. Propagate the 16 computational operator basis elements or otherwise extract Kraus operators, then compare to the corrected ideal CZ. | medium |
| The baseline vertical line and analytical value are ambiguous because of lifetime convention. | Issue text says Evered operating point `τ≈230 µs`; current `src/params.py` deliberately uses ARC radiative-only `|70S_1/2>` lifetime at `T=0`, about `373.9 µs`. At current params, `ε_decay≈5.85e-4`; at `τ=230 µs`, `ε_decay≈9.51e-4`. | The issue may intentionally use `τ≈230 µs` as an approximate experimental effective lifetime, while code uses a physics-model lifetime. No issue comment resolves it. | Maintainer chooses either: use `get_rydberg_params()` as source of truth for the baseline line, or add a clearly named Evered effective lifetime constant for the figure/HTML. | medium |
| Environment/setup can masquerade as implementation failure. | Default system Python lacks QuTiP; full tests pass in `/opt/homebrew/Caskroom/miniconda/base/envs/qutip`. | `requirements.txt` exists and a clean venv should work; this is not a code design blocker. | Run validation inside the conda env or install `requirements.txt` into a local environment. | medium |

No reproduction details are needed to start; this is an implementation issue. The only missing decision is the baseline lifetime/provenance for the figure’s vertical line and HTML wording.

## 5. Decision

- Chosen likely cause: unimplemented decay vertical slice, with two implementation traps: correct open-system fidelity extraction and baseline lifetime labeling.
- Recommended fix or investigation-first plan:
  1. Implement the analytical formula first in `src/analytical.py`:

     ```python
     epsilon_decay(omega, tau) = (1 / tau) * 7*pi/(4*omega)
     ```

     Keep units explicit: if `omega` is rad/µs, `tau` must be µs, or accept `gamma` directly in µs⁻¹.

  2. Implement `src/errors/decay.py` as a Lindblad channel simulation using `qutip.mesolve`, not as a population-only shortcut. The cleanest bounded model is a common two-atom basis with dark states and infinite blockade, for example:

     ```text
     DECAY_BASIS = ("00", "0g", "0r", "g0", "gg", "gr", "r0", "rg")
     ```

     This is the 3×3 two-atom dark/ground/Rydberg basis with `|rr⟩` projected out. Hamiltonian couplings then match the ideal blockade model while permitting arbitrary computational superpositions:
     - atom-1 pulse: `|g0⟩↔|r0⟩`, `|gg⟩↔|rg⟩`; omit `|gr⟩↔|rr⟩`.
     - atom-2 pulse: `|0g⟩↔|0r⟩`, `|gg⟩↔|gr⟩`; omit `|rg⟩↔|rr⟩`.
     - atom-1 collapse: `sqrt(γ)(|g0⟩⟨r0| + |gg⟩⟨rg|)`.
     - atom-2 collapse: `sqrt(γ)(|0g⟩⟨0r| + |gg⟩⟨gr|)`.

  3. For numerical fidelity, propagate the 16 computational operator basis elements `|i⟩⟨j|` through the three `mesolve` segments, project the output map back to the four computational basis states, reconstruct Kraus operators from the Choi matrix, apply the same local-Z correction used by the ideal gate, and call `pedersen_fidelity_kraus()`. Basis-input density matrices can still be returned for diagnostics, but they are not enough for `F_avg`.

  4. Implement a decay sweep over `γ` on a log grid with at least 20 points over at least two orders of magnitude around the chosen baseline. Save data in a simple, inspectable format (`.csv` or `.npz`) before plotting.

  5. Generate `figures/fidelity_vs_decay_rate.png` with numerical dots, analytical curve `F≈1-ε_decay`, and a vertical baseline line. If PNGs remain ignored by `.gitignore`, document that the script regenerates it; if the HTML deliverable requires committed assets, force-add the PNG intentionally rather than changing ignore rules broadly.

  6. Add the Rydberg Decay section to `site/index.html` using the existing template, and update any navigation/table-of-contents entries if the scaffold has them.

- Rejected alternatives and why:
  - Reuse `src/protocol.py` branch-by-branch and estimate fidelity from four final density matrices: rejected because it loses coherences and can pass population checks while failing true average gate fidelity.
  - Model finite blockade in the decay slice: rejected because issue #4 owns finite blockade. For issue #3, keep `|rr⟩` projected out so decay is isolated.
  - Rewrite the whole protocol abstraction before adding decay: rejected. A small decay-specific common-basis helper is enough; generalization can wait until multiple error channels prove the need.
  - Silently change `src/params.py` lifetime to 230 µs: rejected unless the maintainer explicitly chooses that convention. It would contradict the current merged tests and comments.
- Assumptions:
  - The current repository convention is source of truth unless the maintainer overrides it: ARC radiative-only `τ≈373.9 µs` for `|70S_1/2>` at `T=0`.
  - The canonical target remains `CZ = diag(1,1,1,-1)`, with local Z phases corrected as in the ideal tests.
  - Solver units remain rad/µs and µs.
- Maintainer decision needed, if any: choose the baseline lifetime displayed as the “Evered 2023 operating point.” Recommended default if unanswered: use `get_rydberg_params()` for numerical baselines and write the HTML caption as “ARC radiative lifetime for the project `|70S_1/2>` baseline,” not “τ≈230 µs.”

Confidence is not low; the implementation path is bounded. The ambiguity is labeling/provenance, not core code structure.

## 6. Implementation handoff

- Ordered implementation steps:
  1. Create `src/errors/__init__.py`.
  2. Create `src/analytical.py` with `epsilon_decay(omega, tau)` and/or `epsilon_decay_from_gamma(omega, gamma)`. Add unit checks and a test for the current baseline value. With current parameters, expect about `5.85e-4`; with `τ=230 µs`, expect about `9.51e-4`.
  3. Create `src/errors/decay.py`:
     - define the common decay basis and computational-basis projection;
     - build atom-selective Hamiltonians in that basis;
     - build collapse operators `L1`, `L2` for `γ`;
     - run the three pulse segments via `qutip.mesolve`;
     - expose basis-input final density matrices as diagnostics;
     - expose a channel/Kraus or process-output path for fidelity.
  4. Extend `src/fidelity.py` only if needed, narrowly:
     - add `correct_local_z` support for Kraus fidelity, or pre-multiply all actual Kraus operators by `LOCAL_Z_PRODUCT` before calling the existing `pedersen_fidelity_kraus()`;
     - add a small helper to convert a propagated process/Choi matrix into Kraus operators if that logic does not fit cleanly in `src/errors/decay.py`.
  5. Create `src/sweeps.py` with a reusable but simple decay sweep function. Do not over-generalize for future channels beyond what decay needs now.
  6. Add `tests/test_decay.py`:
     - `γ=0` gives `1-F < 1e-10` after local-Z correction;
     - an artificially large `γ` lowers fidelity significantly;
     - baseline numerical error matches `epsilon_decay()` within 50%;
     - sweep returns at least 20 points and covers at least two decades.
  7. Create or update `scripts/run_sweeps.py` to run the decay sweep and save tabular data.
  8. Create or update `scripts/generate_figures.py` to create `figures/fidelity_vs_decay_rate.png`, following the visual style in `scripts/run_ideal.py`.
  9. Insert the Rydberg Decay section into `site/index.html` after the ideal-gate section and before the template comment. Include plain-language physics, the formula, a `<details>` derivation sketch, the figure, and the required takeaway.
  10. Run targeted tests, figure generation, and then the full test suite in the QuTiP environment.
- Files likely to change:
  - `src/errors/__init__.py` (new)
  - `src/errors/decay.py` (new)
  - `src/analytical.py` (new)
  - `src/sweeps.py` (new)
  - `src/fidelity.py` (possibly narrow extension for open-system local-Z correction / Choi-to-Kraus helper)
  - `scripts/run_sweeps.py` (new)
  - `scripts/generate_figures.py` (new)
  - `tests/test_decay.py` (new)
  - `site/index.html` (modify)
  - `figures/fidelity_vs_decay_rate.png` (generated; currently ignored by `.gitignore`)
- Scope boundaries / what not to touch unless new evidence appears:
  - Do not implement finite blockade, Doppler, scattering, amplitude noise, combined budget, or summary table.
  - Do not refactor `src/protocol.py` unless a tiny shared helper avoids duplication without changing ideal behavior.
  - Do not change `src/params.py` lifetime convention just to match the issue text; make the baseline choice explicit instead.
  - Do not replace the HTML scaffold or style system.
- Regression test shape:

  ```python
  def test_epsilon_decay_current_baseline():
      params = get_rydberg_params()
      eps = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)
      assert np.isclose(eps, 5.85e-4, rtol=0.03)

  def test_decay_zero_gamma_is_ideal():
      result = run_decay_gate(omega=DEFAULT_OMEGA_RAD_PER_US, gamma=0.0)
      assert 1.0 - result.average_gate_fidelity < 1e-10

  def test_decay_baseline_matches_analytical_within_50_percent():
      params = get_rydberg_params()
      result = run_decay_gate(params.omega_rad_per_us, params.rydberg_decay_rate_per_us)
      numerical_error = 1.0 - result.average_gate_fidelity
      analytical_error = epsilon_decay(params.omega_rad_per_us, params.rydberg_lifetime_us)
      assert abs(numerical_error - analytical_error) / analytical_error < 0.5

  def test_decay_sweep_shape_and_range():
      rows = sweep_decay(...)
      assert len(rows) >= 20
      assert rows[-1].gamma / rows[0].gamma >= 100
  ```

- Validation command(s):
  - Verified current full-suite command:

    ```bash
    PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q
    ```

  - After implementation, targeted decay validation should be:

    ```bash
    PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python -m pytest -q tests/test_decay.py
    PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/run_sweeps.py
    PYTHONPATH=. /opt/homebrew/Caskroom/miniconda/base/envs/qutip/bin/python scripts/generate_figures.py
    ```

    These script names are specified by the issue but do not exist yet; first implementation should create them.

  - If using another environment, first install project dependencies from `requirements.txt` and confirm `python -c "import qutip, arc"` succeeds.
- Expected pre-fix signal:
  - `import src.errors.decay` fails because `src/errors/` does not exist.
  - `pytest tests/test_decay.py` cannot run because `tests/test_decay.py` does not exist.
  - `figures/fidelity_vs_decay_rate.png` is absent.
  - `site/index.html` has only the template comment, not an actual Rydberg Decay section.
- Expected post-fix signal:
  - `tests/test_decay.py` passes.
  - `γ=0` average gate fidelity passes the ideal threshold with local-Z correction.
  - Baseline numerical decay error is within 50% of `ε_decay`.
  - Sweep output contains at least 20 rows spanning at least 100× in `γ`.
  - `figures/fidelity_vs_decay_rate.png` exists after running figure generation and shows numerical dots, analytical curve, and baseline vertical line.
  - `site/index.html` contains a real `Rydberg Decay` section and existing HTML scaffold tests still pass.
- Negative/edge validation, if applicable:
  - `gamma < 0` should raise `ValueError`.
  - `omega <= 0` or non-finite values should raise `ValueError`.
  - Very large `γ` should not be part of the production figure if it makes the perturbative analytical curve nonsensical; keep the acceptance “large error drops fidelity” test separate from the publication sweep if needed.
  - The projected computational map may be trace-decreasing for large `γ` because residual Rydberg population/leakage can remain at final time. The Kraus/process code should handle this without assuming trace preservation.

## 7. Risks, edge cases, rollback

- Risks:
  - Process-fidelity bug: using only four basis input density matrices loses coherences and can overstate fidelity. Mitigation: propagate all 16 computational operator basis elements or otherwise construct Kraus operators from the channel.
  - Local-Z convention bug: the raw ideal gate is not canonical CZ. Mitigation: apply `LOCAL_Z_PRODUCT` to actual Kraus/operators before comparing with `CZ_TARGET`, or add a tested `correct_local_z` path for Kraus fidelity.
  - Lifetime-labeling bug: figure says “Evered τ≈230 µs” while code uses `τ≈373.9 µs`. Mitigation: make the baseline source explicit and get a maintainer decision if the label must say Evered 2023.
  - Solver/runtime risk: sweep points with many `mesolve` process runs can be slow. Mitigation: start with modest time grids and targeted tests; only increase resolution if numerical/analytical agreement is noisy.
  - Generated figure tracking risk: `.gitignore` ignores PNG figures. Mitigation: either force-add the publication PNG intentionally or document that it is regenerated by scripts; do not broadly unignore all generated files without maintainer intent.
- Edge cases:
  - `γ=0`: collapse list should be empty or zero-safe; must reduce exactly to the ideal local-Z-corrected gate.
  - `γ` so large that perturbation theory fails: numerical points should depart from the analytical line; this is acceptable if the plot and tests distinguish perturbative baseline from stress tests.
  - Residual Rydberg population at final time for high decay/noisy solver settings: account for leakage in the projected channel and tests.
  - HTML without generated PNG: keep or add a graceful fallback/caption similar to the ideal population figure.
- Rollback/revert notes:
  - This slice should be mostly additive. If it causes problems, revert the single implementation commit or remove the new files plus the `site/index.html` decay section. It should not require rolling back the ideal gate or HTML scaffold.

## 8. Open questions

- Question: Should the decay figure’s vertical line use current ARC radiative-only `τ≈373.9 µs` or issue-stated `τ≈230 µs`?
  - What was found: `src/params.py` and `tests/test_params.py` currently encode the ARC radiative-only `|70S_1/2>` convention at `T=0`, about 373.9 µs. The issue text asks for an Evered operating point of about 230 µs.
  - Why it is ambiguous: both values can be physically defensible if one is an ARC model lifetime and the other is an experimental/effective lifetime, but the repo has no separate constant or comment for a 230 µs value.
  - Options:
    1. Use `get_rydberg_params()` for the vertical line and label it as the project ARC baseline.
    2. Add a named `EVERED_EFFECTIVE_RYDBERG_LIFETIME_US = 230.0` for figure/HTML reference only, while leaving ARC params unchanged.
    3. Change `src/params.py` to a different ARC temperature/model, but only if tests and docs are updated intentionally.
  - Recommended default if unanswered: option 1. It is consistent with merged code and tests.

- Question: Should `figures/fidelity_vs_decay_rate.png` be committed even though `figures/*.png` is ignored?
  - What was found: `.gitignore` ignores generated PNGs; the issue acceptance asks for a saved figure and the HTML wants to embed it.
  - Why it is ambiguous: generated artifacts are often ignored, but the HTML deliverable is reader-facing and can break visually without committed figures.
  - Options:
    1. Keep PNG ignored and rely on `scripts/generate_figures.py`.
    2. Force-add this publication figure with `git add -f figures/fidelity_vs_decay_rate.png`.
    3. Embed the PNG as base64 in `site/index.html`.
  - Recommended default if unanswered: option 1 for implementation PRs unless the maintainer requires the HTML to render fully from a fresh checkout without running scripts; then option 2. Avoid option 3 for now.

- Question: Where should Choi/process-to-Kraus reconstruction live?
  - What was found: `src/fidelity.py` already has `pedersen_fidelity_kraus()` but no reconstruction helper. `src/errors/decay.py` will need process reconstruction if it propagates 16 operator basis elements.
  - Why it is ambiguous: the helper may be useful for later Lindblad channels, but premature abstraction is how small projects get moldy.
  - Options:
    1. Keep a private helper in `src/errors/decay.py` now.
    2. Add a small public `kraus_from_process()` helper to `src/fidelity.py` with tests.
  - Recommended default if unanswered: option 2 only if tests cover it immediately; otherwise option 1 is fine for the first slice.

- Question: Should `src/sweeps.py` be generic now?
  - What was found: future issues will add more sweeps, but only decay exists after this slice.
  - Why it is ambiguous: a generic sweep framework could reduce duplication later, but current repo has no convention to preserve.
  - Options:
    1. Implement `sweep_decay()` directly and refactor once a second error channel proves the pattern.
    2. Build a generic `run_sweep(simulator, params, analytical_fn)` now.
  - Recommended default if unanswered: option 1 with a simple row dataclass or dictionaries. Earn the abstraction.
