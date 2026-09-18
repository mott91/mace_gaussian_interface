"""RED tests for Phase 23 HTMLReportGenerator Plotly overhaul."""

import re
from html.parser import HTMLParser

from mace_gaussian.analysis.html_report_generator import HTMLReportGenerator


class _WellFormedChecker(HTMLParser):
    def __init__(self):
        super().__init__()
        self.stack = []
        self.errors = []

    def handle_starttag(self, tag, attrs):
        if tag not in ("br", "img", "meta", "link", "input", "hr"):
            self.stack.append(tag)

    def handle_endtag(self, tag):
        if self.stack and self.stack[-1] == tag:
            self.stack.pop()


def _generate(fake_analysis_results, tmp_path, mode="anharmonic"):
    gen = HTMLReportGenerator(
        molecule_name="water",
        output_dir=tmp_path,
        mode=mode,
        plotly_js="cdn",
        bandwidth_fwhm=10.0,
    )
    fake_analysis_results["mode"] = mode
    fake_analysis_results["output_dir"] = str(tmp_path)
    gen.generate_report(fake_analysis_results)
    return (tmp_path / "report.html").read_text(encoding="utf-8")


def test_per_method_section_has_plotly_div(fake_analysis_results, tmp_path):
    html = _generate(fake_analysis_results, tmp_path)
    assert "plotly-graph-div" in html or "plotly_" in html or 'class="plotly' in html


def test_plotlyjs_emitted_once(fake_analysis_results, tmp_path):
    html = _generate(fake_analysis_results, tmp_path)
    # CDN mode: script tag referencing plotly CDN should appear exactly once
    cdn_count = html.count("cdn.plot.ly/plotly")
    assert cdn_count == 1, f"Expected plotly CDN referenced exactly once, got {cdn_count}"


def test_mode_flag_harmonic_skips_overtones(fake_analysis_results, tmp_path):
    html_h = _generate(fake_analysis_results, tmp_path / "h", mode="harmonic")
    html_a = _generate(fake_analysis_results, tmp_path / "a", mode="anharmonic")
    # anharmonic contains overtones section; harmonic does not
    assert "overtone" in html_a.lower() or "combination" in html_a.lower()
    assert "overtone" not in html_h.lower()


def test_css_shared_across_modes(fake_analysis_results, tmp_path):
    html_h = _generate(fake_analysis_results, tmp_path / "h", mode="harmonic")
    html_a = _generate(fake_analysis_results, tmp_path / "a", mode="anharmonic")
    css_h = re.search(r"<style>(.*?)</style>", html_h, re.DOTALL).group(1)
    css_a = re.search(r"<style>(.*?)</style>", html_a, re.DOTALL).group(1)
    assert css_h == css_a


def test_timing_section_renders(fake_analysis_results, tmp_path):
    html = _generate(fake_analysis_results, tmp_path)
    # Timing values from fixture should appear
    assert "11.5" in html or "58.0" in html or "speedup" in html.lower()


def test_experimental_trace_present_when_experimental_available(fake_analysis_results, tmp_path):
    html = _generate(fake_analysis_results, tmp_path)
    assert "NIST" in html or "xperimental" in html


def test_html_well_formed(fake_analysis_results, tmp_path):
    html = _generate(fake_analysis_results, tmp_path)
    parser = _WellFormedChecker()
    parser.feed(html)
    # Remaining stack should be empty (all tags closed)
    assert len(parser.stack) == 0, f"Unclosed tags: {parser.stack}"


def test_html_escapes_molecule_name(fake_analysis_results, tmp_path):
    fake_analysis_results["molecule"] = "<script>alert(1)</script>"
    gen = HTMLReportGenerator(
        molecule_name="<script>alert(1)</script>",
        output_dir=tmp_path,
        mode="anharmonic",
        plotly_js="cdn",
        bandwidth_fwhm=10.0,
    )
    gen.generate_report(fake_analysis_results)
    html = (tmp_path / "report.html").read_text()
    # Raw script tag must NOT appear unescaped in output
    assert "<script>alert(1)</script>" not in html
    assert "&lt;script&gt;" in html
