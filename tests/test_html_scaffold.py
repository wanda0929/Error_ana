from pathlib import Path
import re


HTML_PATH = Path(__file__).resolve().parents[1] / "site" / "index.html"


def _html() -> str:
    return HTML_PATH.read_text(encoding="utf-8")


def test_html_scaffold_exists_with_required_sections():
    html = _html()
    assert HTML_PATH.exists()
    assert len(html.splitlines()) >= 300
    assert 'id="sec-introduction"' in html
    assert 'id="sec-ideal-gate"' in html
    assert "Rydberg state" in html
    assert "Qubit encoding" in html
    assert "Blockade" in html
    assert "Controlled-Z truth table" in html
    assert "π–2π–π" in html


def test_html_uses_mathjax_and_inline_project_css():
    html = _html()
    assert 'id="MathJax-script"' in html
    assert "tex-chtml.js" in html
    assert "cdn.tailwindcss.com" not in html
    assert "<style>" in html
    assert "--dark-900" in html


def test_html_has_ideal_figure_and_graceful_fallback():
    html = _html()
    assert '../figures/population_dynamics.png' in html
    assert 'id="population-fallback"' in html
    assert "showPopulationFallback" in html
    assert "$F_{\\mathrm{avg}} = 1.0000000000$" in html
    assert "infidelity below $10^{-10}$" in html


def test_html_interactive_and_downstream_template_scaffold():
    html = _html()
    assert 'id="selection-toolbar"' in html
    assert 'id="comment-modal"' in html
    assert 'id="notes-panel"' in html
    assert "Copy as LLM prompt" in html
    assert html.count("<details>") >= 3
    assert "ERROR-SOURCE SECTION TEMPLATE FOR ISSUES #3-#7" in html
    assert 'id="sec-error-{name}"' in html


def test_svg_fallbacks_do_not_contain_raw_mathjax_delimiters():
    html = _html()
    svg_blocks = re.findall(r"<svg[\s\S]*?</svg>", html)
    assert svg_blocks
    assert all("$" not in svg for svg in svg_blocks)
