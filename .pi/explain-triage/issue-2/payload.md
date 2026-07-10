## TL;DR

Issue #2 asks for the creation of `site/index.html` — the **single reader-facing deliverable** of the entire `Error_ana` repository. This repository is a quantum-computing project that analyzes error sources in a Rydberg-atom CZ gate. Every other file in the repo (simulation scripts, figures, tests) exists only to produce data that eventually appears inside this HTML file.

Issue #2 specifically creates:

1. **The scaffold** — all shared infrastructure (CSS, JavaScript, math rendering, interactive features) that every subsequent issue depends on.
2. **The Introduction section** — defines Rydberg states, qubit encoding, the blockade effect, and the CZ gate for a physicist who hasn't worked with Rydberg atoms before.
3. **The Ideal Gate section** — walks through the π–2π–π pulse protocol, shows a truth table, embeds a population-dynamics figure from issue #1, and confirms that fidelity equals 1 in the ideal case.

Issues #3 through #8 each add one section to this scaffold. **Nothing else can proceed until this file exists.** The triage plan recommends building it with the Tailwind CSS framework, MathJax v3 for equation rendering, and a suite of interactive features (text-selection quoting, comment threads, collapsible derivations) that match the repository owner's established conventions in a sibling project.

---

## Decision map

Below is a summary of every significant design choice recorded in the triage document, the alternatives considered, and the rationale for each recommendation.

| Decision | Chosen option | Rejected alternative(s) | Rationale |
|---|---|---|---|
| **Math rendering** | MathJax v3 via CDN | (A) Hand-rendered HTML spans with CSS classes; (B) KaTeX | The project will have ~20 complex equations (Lindblad operators, fidelity formulas, scaling laws). Hand-rendering is fragile. MathJax v3 is async, loads only needed modules (~100 KB), and is specified in `AGENTS.md`. |
| **CSS framework** | Tailwind CSS via CDN, dark theme | Custom CSS or another framework | Matches the owner's existing convention in the sibling repo `Finite_blockade_pulse`, which has three explainer HTML files already using Tailwind + a dark color palette. |
| **Figure embedding** | Relative path (`../figures/population_dynamics.png`) | Base64 inline embedding | Relative paths keep the HTML small and debuggable. Figures are regeneratable from scripts, so there's no portability concern. |
| **Section template for errors** | Rigid, documented HTML template (inserted as an HTML comment) | Freeform per error source | Five error-source issues (#3–#7) must each add an identically structured section. A commented template prevents structural drift. |
| **Interactive features** | Full notebook-style (selection quoting, comment modal, notes panel, "Copy as LLM prompt") | Minimal (collapsibles only) | The owner's sibling project already uses these features. The `interactive-html-explainer` **skill** (a reusable build recipe available in the Pi toolchain) produces exactly this type of file. |

---

## State machine trace

Some terms for context:

- **Triage branch**: A Git branch created to investigate an issue, gather evidence, and write a plan — without making code changes. Here it is `triage/issue-2-2-html-scaffold-introduction-ideal-gate-section`.
- **Resolve branch**: A subsequent branch where the plan is actually implemented. This has not been created yet.
- **Skill**: A reusable recipe in the Pi extension framework (a VS Code / CLI automation layer) that knows how to produce a particular kind of artifact. The `interactive-html-explainer` skill produces self-contained dark-themed HTML files.

The triage proceeded through these states:

1. **Issue intake** — Issue #2 was read. It requests an HTML scaffold plus two content sections.
2. **Dependency check** — Issue #1 (population dynamics figure) was identified as a blocker. The figure file `figures/population_dynamics.png` was found to exist locally and was visually verified as correct, but it is untracked in Git (meaning issue #1 hasn't been merged yet).
3. **Convention survey** — Three existing HTML explainer files in the sibling repo `Finite_blockade_pulse` were inspected. Their conventions (Tailwind CDN, dark theme, color palette, interactive features) were catalogued.
4. **AGENTS.md review** — The project's own instructions specify MathJax for equations. A tension was noted: the sibling repo uses hand-rendered math, but AGENTS.md says MathJax. The triage resolves this in favor of MathJax because of equation volume.
5. **Skill match** — The `interactive-html-explainer` skill was identified as the right build tool, since it produces exactly the kind of file described.
6. **Plan written** — A six-step implementation plan was recorded with explicit scope boundaries (do not create error-source sections, do not modify Python code).
7. **Triage committed** — The plan was saved to `docs/triage/issue-2.md` on the triage branch.

---

## Evidence

The triage document cites the following concrete evidence. Each item is something the triager observed in the repository or its environment, not a guess.

### 1. The file does not exist yet
`site/index.html` is entirely new. The pre-fix signal is that `open site/index.html` fails.

### 2. Population dynamics figure verified
The file `figures/population_dynamics.png` exists locally (untracked). The triager describes its content in detail: Rydberg population vs. time for all four two-qubit basis states during the π–2π–π protocol, with three pulse regions marked by dashed lines. The behavior of each state is consistent with the physics (|00⟩ flat, |01⟩ shows a 2π rotation, |10⟩ and |11⟩ show π round-trips with |11⟩ exhibiting blockade).

### 3. Sibling-repo conventions catalogued
A detailed table of conventions was extracted from `Finite_blockade_pulse/czz-cole-gate-explainer.html`:

- Tailwind via `cdn.tailwindcss.com`
- Dark mode with specific hex colors (`#0f1117` for background)
- Semantic colors for physics concepts: `rydberg` (pink), `ground` (green), `blockade` (amber)
- Custom accordion with `max-height` CSS transition
- Interactive features: text-selection toolbar, comment modal, notes panel, localStorage persistence

### 4. AGENTS.md specification
AGENTS.md says to use "MathJax CDN for equations." This conflicts with the sibling repo's hand-rendered math but is the better choice for this project due to equation complexity.

### 5. Skill availability
The `interactive-html-explainer` skill exists at a known filesystem path and is documented. It produces files matching the needed conventions.

---

## Implementation plan

The triage prescribes six ordered steps. Only one file is created: `site/index.html`.

### Step 1 — Load the `interactive-html-explainer` skill
This skill (a reusable build recipe in the Pi extension framework) handles the boilerplate for single-file interactive HTML explainers. It is the starting point rather than writing from scratch.

### Step 2 — Set up infrastructure in `site/index.html`
This is the "scaffold" layer — built once and never modified by downstream issues.

| Component | Details |
|---|---|
| **Tailwind CSS** | CDN link, dark-mode class on `<html>`, custom color config matching sibling repo |
| **MathJax v3** | Async CDN load of `tex-chtml.js`; inline math delimiters `\( ... \)`, display math `$$ ... $$` |
| **Semantic colors** | `rydberg: #f472b6`, `ground: #34d399`, `blockade: #fbbf24` |
| **Interactive JS** | Text-selection toolbar, comment modal, notes panel with localStorage, "Copy as LLM prompt" button |
| **Navigation** | Sticky top bar with page title, "Notes" button, "Print" button |
| **Footer** | Export notes, print/PDF, credit line |
| **Print styles** | `@media print` CSS switching to white background / black text |

### Step 3 — Write Introduction section (H-01)
Defines all physics terms from scratch for a QM-literate reader without Rydberg background:

1. Rydberg states (high-n atomic excitations, huge dipole, strong interactions)
2. Qubit encoding (|0⟩ = dark state, |1⟩ = |g⟩, |r⟩ is auxiliary)
3. Two-photon excitation (780 nm + 480 nm for ⁸⁷Rb)
4. Blockade effect (interaction energy U >> Rabi frequency Ω prevents double excitation)
5. CZ gate definition (diag(1, 1, 1, −1) with truth table)

Key MathJax equations: van der Waals interaction `V(R) = -C_6/R^6`, blockade strength, CZ matrix.

### Step 4 — Write Ideal Gate section (H-02)

1. Three colored step-boxes for the π–2π–π protocol
2. Truth table showing phases for all four input states
3. Embedded figure: `<img src="../figures/population_dynamics.png">` with descriptive caption
4. Statement that ideal fidelity F = 1.0 (infidelity < 10⁻¹⁰)
5. Collapsible `<details>` block explaining the Pedersen average gate fidelity formula

### Step 5 — Insert error-source template as HTML comment
A commented-out HTML block showing the exact structure that issues #3–#7 must follow:

```
<section id="sec-error-{name}"> → <h2> → callout → physics paragraphs →
  MathJax formula box → collapsible derivation → figure → takeaway
</section>
```

This prevents structural drift across five independent implementation efforts.

### Step 6 — Manual testing
See **Validation** below.

---

## Validation

There are no automated tests for this deliverable — it is a presentation artifact. Validation is manual:

| Check | Method | Expected result |
|---|---|---|
| **Renders in browsers** | Open `site/index.html` in Chrome, Firefox, Safari | All content visible, styled correctly |
| **MathJax works** | Ctrl+F for `\frac` in the rendered page | 0 results (all LaTeX compiled to rendered math) |
| **Collapsibles work** | Click a `<details>` / `<summary>` block | Content expands/collapses smoothly |
| **Mobile responsive** | Resize viewport to 375px width | No horizontal scrollbar |
| **Selection toolbar** | Select any text | Floating toolbar appears with "Ask about this" |
| **Figure loads** | Scroll to ideal gate section | `population_dynamics.png` visible with caption |
| **File size sanity** | `wc -l site/index.html` | 300+ lines (non-trivial) |

**Pre-fix signal**: `site/index.html` does not exist.  
**Post-fix signal**: Dark-themed page with two content sections, rendered equations, working interactivity.

**Negative / edge-case validation**:
- No network → MathJax and Tailwind CDN fail → page is unstyled with raw LaTeX. This is **expected and accepted** per AGENTS.md (CDN is the only external dependency).
- Figure path `../figures/population_dynamics.png` must resolve when opened as a `file://` URL and when served from a local HTTP server. Both should be tested.

---

## Caveats and open questions

### Blocker dependency
Issue #1 must be merged before this issue's figure reference works in the committed repo. The figure exists locally and has been visually verified, but it is **untracked** — meaning someone could implement issue #2 before #1 is merged and the figure would load locally but fail for others who clone the repo.

### AGENTS.md vs. sibling-repo tension
AGENTS.md says "MathJax CDN." The sibling repo uses hand-rendered math. The triage resolves this in favor of MathJax, but this is a judgment call that could be overridden.

### CDN fragility
Both Tailwind and MathJax are loaded from CDNs. If CDNs are unreachable, the page degrades to unstyled HTML with raw LaTeX. The triage accepts this as a known limitation. Pinning CDN versions is recommended but not mandated.

### Template drift risk
Five separate error-source issues (#3–#7) will each add a section to this scaffold. The HTML-comment template is the only mechanism preventing structural inconsistency. There is no automated linting or schema enforcement.

### No automated regression tests
The triage explicitly states "No `pytest` for HTML — validation is visual." If the project later wants CI, an HTML validator or Lighthouse check would need to be added separately.

### Triage document is truncated
The MathJax configuration block in Section 6, Step 2 is cut off mid-sentence (`MathJax = { tex: { inlineMath: [['`). The full configuration is inferrable from the earlier discussion (inline delimiters `\(…\)`, display delimiters `$$…$$`), but the triage document is technically incomplete at that point.

---

## Revision targets

These are the specific decisions a human reviewer might reasonably want to change or annotate:

1. **MathJax v3 vs. KaTeX vs. hand-rendered math** — The triage chose MathJax v3. KaTeX renders faster but supports fewer LaTeX commands. A reviewer might prefer KaTeX if equation complexity is lower than estimated, or hand-rendered math if CDN dependence is unacceptable.

2. **Relative figure paths vs. Base64 embedding** — The triage chose relative paths (`../figures/`). A reviewer concerned about single-file portability (e.g., emailing the HTML) might prefer Base64.

3. **Full interactive notebook features vs. minimal HTML** — The triage includes text-selection quoting, comment threads, notes panel, and "Copy as LLM prompt." A reviewer might consider this over-engineered for what is essentially a physics explainer document.

4. **Scope of the Introduction section** — The triage specifies that two-photon excitation details (780 nm + 480 nm wavelengths, counter-propagating beams, Doppler suppression) should be included. A reviewer might consider this too detailed for an introduction, or not detailed enough.

5. **Rigid error-source template vs. flexible sections** — The triage mandates a fixed HTML structure for all five error channels. A reviewer might argue that different error sources (e.g., spontaneous emission vs. Doppler broadening) have different enough pedagogical needs to warrant flexibility.

6. **CDN version pinning** — The triage recommends but does not require pinning Tailwind and MathJax CDN versions. A reviewer might want to make this mandatory.

7. **Dependency ordering** — The triage says issue #1 must be merged first. A reviewer might decide to proceed in parallel, using a placeholder image and swapping in the real figure later.

8. **Semantic color palette** — The specific hex values for `rydberg`, `ground`, and `blockade` are copied from the sibling repo. A reviewer might want different colors for this project.