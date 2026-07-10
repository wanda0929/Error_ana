# Issue #2: HTML scaffold + introduction + ideal gate section

## 1. Executive summary

- **Actionability**: actionable
- **Symptom**: No HTML file exists yet. This issue creates the presentation shell and first two content sections.
- **User impact / severity**: Blocks issues #3–#8 (every error-source slice adds an HTML section to this scaffold).
- **Recommended next step**: Use the `interactive-html-explainer` skill, following the Tailwind + dark-theme conventions established in the owner's other repo (`Finite_blockade_pulse`). Add MathJax v3 for LaTeX equations.
- **Confidence**: high
- **Main blocker**: Issue #1 must be merged first — the population dynamics figure (`figures/population_dynamics.png`) is needed for the ideal gate section. The figure exists locally and has been visually verified as correct during triage.

## 2. Affected-path context

### Affected subsystem

`site/index.html` — a single new file. No existing code is modified.

### Repo vocabulary / components

| Term | Meaning |
|---|---|
| `site/index.html` | The primary deliverable — self-contained interactive HTML explainer |
| Population dynamics figure | `figures/population_dynamics.png` from issue #1 — embedded in the ideal gate section |
| Error-source section template | A repeatable HTML structure that issues #3–#7 each instantiate for their error channel |
| Collapsible derivation | A click-to-expand `<details>/<summary>` block containing the mathematical derivation sketch |

### Current relevant control/data flow

```
figures/population_dynamics.png  →  embedded in site/index.html
(from issue #1)                     (this issue creates)
                                         ↓
                                    issues #3–#7 each append
                                    one error-source section
                                         ↓
                                    issue #8 adds combined
                                    budget + summary table
```

### Why the key files matter

`site/index.html` is the only reader-facing artifact in the entire project. Every other file (simulation code, scripts, tests) exists to produce data that ultimately appears here. The scaffold must be right because 6 downstream issues add content to it.

### First file a new implementer should open

`AGENTS.md` → "HTML explainer requirements" section, then inspect `Finite_blockade_pulse/czz-cole-gate-explainer.html` as the style reference.

## 3. Observed facts

### Issue evidence

- **Reporter claim**: Create HTML shell with infrastructure (CSS, MathJax, collapsibles, responsive) plus introduction and ideal gate sections.
- **Expected behavior**: `site/index.html` opens in browser, equations render, collapsibles work, mobile-responsive, intro defines all terms, ideal gate shows protocol + figure + F=1.
- **Acceptance criteria** (from issue):
  1. Opens in Chrome/Firefox/Safari, all content visible
  2. MathJax equations render (no raw LaTeX)
  3. Collapsible sections work
  4. Responsive at 375px width
  5. Introduction defines: Rydberg state, blockade, qubit encoding, CZ truth table
  6. Ideal gate shows: protocol steps, truth table, population dynamics figure, F=1
  7. QM-literate reader can follow without external references

### Repository evidence

**Population dynamics figure verified** (exists locally at `figures/population_dynamics.png`, untracked):

The figure shows Rydberg population vs. time for all 4 basis states {|00⟩, |01⟩, |10⟩, |11⟩} during the π–2π–π protocol. Three pulse regions marked with dashed lines. Visually correct:
- |00⟩: flat at 0 (both dark)
- |01⟩: 0 → 1 → 0 during 2π₂ pulse only
- |10⟩: 0 → 1 during π₁, stays at 1, returns during final π₁
- |11⟩: same as |10⟩ — atom 2 is blocked, only atom 1's Rydberg population shows

Publication quality: labeled axes, legend, clean styling.

**Established HTML convention** from `Finite_blockade_pulse/` (3 existing explainer files):

| Convention | Details |
|---|---|
| CSS framework | Tailwind CSS via CDN (`cdn.tailwindcss.com`) |
| Theme | Dark mode (`class="dark"`, custom color palette) |
| Color palette | `dark-900: #0f1117`, `accent: blue`, semantic colors for `rydberg`, `ground`, `blockade` |
| Math rendering | Hand-rendered with `.math-block` CSS class (NOT MathJax) |
| Collapsibles | Custom accordion with `max-height` transition |
| Layout | `max-w-4xl mx-auto`, card-based sections |
| Typography | Inter font (body), Times New Roman (math) |
| Interactive features | Text-selection quoting, local commenting, "Copy as LLM prompt", notes panel |
| Navigation | Sticky top nav with section title + utility buttons |

**Tension with AGENTS.md**: The AGENTS.md says "MathJax CDN for equations" but the existing explainer convention uses hand-rendered math. MathJax is the better choice here because:
1. The error analysis has many LaTeX equations (Lindblad operators, fidelity formulas, scaling laws)
2. Hand-rendering `\frac{7\pi}{4\Omega\tau}` is fragile and ugly
3. MathJax v3 is lightweight (async, only loads what's needed)

**Recommendation**: Use Tailwind + dark theme (matching the owner's convention) but **add MathJax v3** for equations (matching the AGENTS.md spec). This combines the best of both.

**Interactive HTML explainer skill available** at `/Users/wanda/pi-excalibur/skills/interactive-html-explainer/SKILL.md`. This skill produces exactly the kind of file needed: single self-contained HTML with Tailwind, dark theme, text-selection quoting, commenting threads, collapsibles, and accordion sections. It matches the existing convention.

**Existing tests/fixtures**: None for HTML.

## 4. Root-cause candidates

N/A — greenfield implementation. The relevant analysis is design decisions.

| Decision | Options | Recommendation | Rationale |
|---|---|---|---|
| Math rendering | (A) Hand-rendered `.math-block` spans (B) MathJax v3 CDN (C) KaTeX CDN | **(B) MathJax v3** | Many complex equations needed downstream; hand-rendering would be fragile. MathJax v3 is async and only ~100KB for tex-chtml. |
| Figure embedding | (A) Base64 inline (B) Relative `<img>` path to `../figures/` | **(B) Relative path** | Simpler, debuggable. Base64 bloats the HTML. The figure file is regeneratable from scripts. |
| Section template | (A) Freeform per error source (B) Rigid template with fixed slots | **(B) Rigid template** | 5 error slices must plug in identically. A documented template prevents drift. |
| Interactive features | (A) Minimal (just collapsibles) (B) Full notebook (selection quoting, notes, LLM prompts) | **(B) Full notebook** | Matches the owner's established pattern in Finite_blockade_pulse. The skill handles this. |

## 5. Decision

- **Chosen approach**: Build with the `interactive-html-explainer` skill using Tailwind dark theme + MathJax v3. Follow the structural conventions from `Finite_blockade_pulse/czz-cole-gate-explainer.html`.

- **Recommended implementation**:

  The HTML file has three layers: infrastructure, content, and the error-source template.

  ### Layer 1: Infrastructure (built once, never touched again)

  - **Head**: Tailwind CDN, MathJax v3 CDN (`tex-chtml.js`), custom dark-mode color config, responsive meta viewport
  - **CSS**: Dark theme palette matching existing explainers (`dark-900`, `accent-400`, semantic `rydberg`/`ground`/`blockade` colors), `.accordion-content` transitions, print styles
  - **JS**: Text-selection toolbar, comment modal, notes panel, localStorage persistence, "Copy as LLM prompt", collapsible toggle, MathJax auto-typeset on accordion open
  - **Nav**: Sticky top bar with page title + "Notes" + "Print" buttons
  - **Footer**: Export notes, print/PDF, credit line

  ### Layer 2: Content sections (this issue fills sections 1–2)

  **Section 1 — Introduction (H-01)**:

  Must define from scratch (audience: QM-literate physicist, no Rydberg background):

  1. **Rydberg states**: atoms excited to high principal quantum number n. Orbital radius scales as n², so electron is far from nucleus → huge electric dipole → strong long-range interactions between nearby Rydberg atoms.

  2. **Qubit encoding**: computational |0⟩ is a dark state (doesn't couple to the Rydberg laser). Computational |1⟩ = |g⟩ (ground state that does couple). The Rydberg state |r⟩ is auxiliary — used during the gate but population returns to {|0⟩, |1⟩} afterward.

  3. **Two-photon excitation**: 780 nm (5S → 5P) + 480 nm (5P → nS) for ⁸⁷Rb. Counter-propagating beams give small effective k_eff, suppressing Doppler sensitivity.

  4. **Blockade effect**: When atom 1 is in |r⟩, it shifts atom 2's Rydberg level by the interaction energy U. If U >> Ω (Rabi frequency), atom 2 cannot be excited → "blocked." This is the mechanism that makes the CZ gate work.

  5. **CZ gate**: Controlled-Z applies a π phase to |11⟩ and leaves |00⟩, |01⟩, |10⟩ unchanged. Target unitary: diag(1, 1, 1, −1). Include the truth table.

  Key equations to typeset in MathJax:
  - `V(R) = -C_6 / R^6` (van der Waals interaction)
  - `U = C_6 / R^6` (blockade strength)
  - `CZ = \text{diag}(1, 1, 1, -1)`

  **Section 2 — Ideal gate (H-02)**:

  1. **π–2π–π protocol**: Step-by-step walkthrough. Three colored boxes (matching the pulse-step styling in existing explainers).
     - Step 1: π-pulse on atom 1 (|g⟩₁ → |r⟩₁)
     - Step 2: 2π-pulse on atom 2 (blocked if atom 1 in |r⟩)
     - Step 3: π-pulse on atom 1 (|r⟩₁ → |g⟩₁)

  2. **Truth table**: What happens to each input state. Show phases.
     - |00⟩ → +1 (both dark)
     - |01⟩ → −1 (2π rotation)
     - |10⟩ → −1 (π + π round trip)
     - |11⟩ → −1 (atom 1 round trip × atom 2 blocked)
     - Net: diag(1, −1, −1, −1) = CZ × (Z₁⊗Z₂) — local Z corrections absorbed.

  3. **Population dynamics figure**: `<img src="../figures/population_dynamics.png">` with descriptive caption.

  4. **F = 1 confirmation**: State that the ideal simulation yields Pedersen average gate fidelity F = 1.0 to machine precision (infidelity < 10⁻¹⁰). Collapsible section explaining Pedersen fidelity formula.

  ### Layer 3: Error-source section template (empty, documented for slices #3–#7)

  Each error-source section follows this exact structure:

  ```html
  <section id="sec-error-{name}" data-section="Error: {Title}" class="mb-16">
    <h2>0X — {Title}</h2>
    <div class="callout"><!-- one-sentence physics summary --></div>

    <!-- Physics explanation: 2-3 paragraphs, plain language -->

    <!-- Analytical formula in MathJax, highlighted box -->

    <details class="derivation">
      <summary>Derivation sketch ▸</summary>
      <!-- Collapsible: key steps of the derivation -->
    </details>

    <!-- Figure: <img src="../figures/fidelity_vs_{param}.png"> -->
    <!-- Caption describing what the reader sees -->

    <div class="takeaway">
      <!-- "This error dominates when ___" -->
    </div>
  </section>
  ```

  Include this as an HTML comment in the file so implementers of #3–#7 know where and how to add their section.

- **Rejected alternatives**:
  - Hand-rendered math (as in Finite_blockade_pulse): rejected because the Error_ana project has ~20 distinct equations. MathJax is necessary.
  - Minimal HTML without interactive features: rejected because the owner's established convention includes text-selection quoting and commenting.
  - Base64 figure embedding: rejected — bloats the file, harder to debug. Relative paths are fine since all figures are regeneratable.

- **Assumptions**:
  - The `interactive-html-explainer` skill is used during implementation.
  - Figures referenced via `../figures/` relative path (HTML lives in `site/`, figures in `figures/`).
  - Issue #1 is merged before this issue's figure reference is finalized.

## 6. Implementation handoff

### Ordered implementation steps

1. **Load the `interactive-html-explainer` skill** — it produces exactly the right kind of file.

2. **Set up infrastructure** in `site/index.html`:
   - Tailwind CDN + dark-mode config (copy color palette from `czz-cole-gate-explainer.html`)
   - MathJax v3: `<script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>`
   - MathJax config: `MathJax = { tex: { inlineMath: [['$','$'], ['\\(','\\)']] } }`
   - Semantic colors: `rydberg: '#f472b6'`, `ground: '#34d399'`, `blockade: '#fbbf24'` (matching existing convention)
   - Selection toolbar, comment modal, notes panel (from the skill)
   - Sticky nav, footer with export/print

3. **Write Introduction section (H-01)**:
   - Define every term before using it
   - Use MathJax for all equations
   - Include a simple energy-level diagram (inline SVG or text-based) showing |g⟩, |r⟩, and the two-photon path
   - End with the CZ truth table

4. **Write Ideal gate section (H-02)**:
   - Three colored step boxes for π–2π–π
   - Truth table with phases
   - `<img src="../figures/population_dynamics.png">` with caption
   - F = 1 statement + collapsible Pedersen fidelity explanation

5. **Add error-source template** as HTML comment (documented for #3–#7)

6. **Test**:
   - Open in Chrome, Firefox, Safari
   - Verify MathJax renders (check a `\frac{}{}` in the Pedersen formula)
   - Verify collapsible sections work
   - Test at 375px viewport width (no horizontal scroll)
   - Verify figure loads from relative path

### Files likely to change

| File | Action |
|---|---|
| `site/index.html` | **Create** — the entire deliverable |

### Scope boundaries

- Do **not** create any error-source sections (those are issues #3–#7)
- Do **not** create the combined error budget or summary table (issue #8)
- Do **not** modify any Python code
- The error-source template is an HTML comment — a documented pattern, not filled content

### Regression test shape

No automated tests — this is a presentation artifact. Manual verification:

```
1. open site/index.html in browser
2. Ctrl+F "\\frac" → 0 results (all MathJax rendered)
3. Click a collapsible derivation → expands smoothly
4. Resize to 375px width → no horizontal scroll
5. Select text → floating toolbar appears → "Ask about this" works
6. Figure loads in ideal gate section
```

### Validation command(s)

```bash
# Open in default browser
open site/index.html

# Check file exists and is non-trivial
wc -l site/index.html  # expect 300+ lines
```

No `pytest` for HTML — validation is visual.

### Expected pre-fix signal

`site/index.html` does not exist → `open site/index.html` fails.

### Expected post-fix signal

`open site/index.html` renders a dark-themed page with:
- Introduction section with MathJax equations
- Ideal gate section with population dynamics figure
- Functional collapsibles, selection toolbar, notes panel

### Negative/edge validation

- Verify MathJax works offline-ish: CDN-dependent, so no network → equations won't render. This is expected and documented in AGENTS.md ("MathJax via CDN is the only external dependency").
- Verify that the figure `../figures/population_dynamics.png` path resolves correctly when opening from `site/index.html`.

## 7. Risks, edge cases, rollback

### Risks

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| MathJax CDN slow or blocked | Low | Medium (equations don't render) | Accept — AGENTS.md permits this as the only external dependency. Add a `<noscript>` fallback note. |
| Figure path breaks when served vs. opened as file | Medium | Low | Use relative path `../figures/`. Test both `file://` and simple HTTP server. |
| Tailwind CDN changes / breaks | Low | Medium | Pin version in CDN URL if available. |
| Error-source template drifts across 5 slices | Medium | Medium | Document the template as a commented HTML block with clear instructions. |

### Edge cases

- **No network**: MathJax and Tailwind won't load. Page will be unstyled with raw LaTeX. Acceptable per AGENTS.md.
- **Very long equations**: MathJax display-mode equations may overflow on mobile. Use `overflow-x: auto` on equation containers.
- **Dark mode printing**: Add `@media print` CSS to switch to white background / black text (existing convention in Finite_blockade_pulse).

### Rollback/revert notes

Single new file. Rollback = `git rm site/index.html`.

## 8. Open questions

No blocking questions. All design decisions are resolvable from existing conventions + AGENTS.md.

The only dependency is issue #1 (population dynamics figure). This is a stated blocker in the issue itself and the figure has been verified to exist locally.
