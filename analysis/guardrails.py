"""Executable enforcement of CLAUDE.md.

Two roles:

1. **Import-time API** used by every analysis script that makes a group
   comparison (`require_mast_subject`, `comparison_order`, `check_model_terms`).
2. **Repo audit** — `python3 analysis/guardrails.py` checks the contract against
   the repository and exits non-zero on violation (CLAUDE.md §9).

The checks are deliberately blunt. They are meant to fail loudly rather than to
be clever.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"

# --------------------------------------------------------------------------
# §2  The subject is mast cells
# --------------------------------------------------------------------------
CONTROL_POPULATIONS = {
    "T/NK", "T cells", "T", "Macrophages", "Macrophage", "DC",
    "Keratinocytes", "Fibroblasts", "VEC", "LEC", "Melanocytes",
    "Pericyte/SMC", "Plasma", "Schwann", "Sweat gland",
    "Cornified keratinocytes",
}

MAST_ALIASES = {"mast", "mast cells", "mast cell", "mc", "mastcells"}


class ContractViolation(AssertionError):
    """Raised when an analysis breaches CLAUDE.md."""


def require_mast_subject(subject: str) -> str:
    """§2 — the subject of every analysis must be mast cells.

    Control populations are legitimate *comparators*, never the subject.
    """
    if str(subject).strip().lower() not in MAST_ALIASES:
        raise ContractViolation(
            f"CLAUDE.md §2: subject of analysis must be mast cells, got {subject!r}. "
            "Control populations appear only to demonstrate specificity."
        )
    return subject


def assert_control_role(populations) -> None:
    """§2 — non-mast populations may only be carried as controls."""
    for p in populations:
        if str(p).strip().lower() in MAST_ALIASES:
            continue
        if p not in CONTROL_POPULATIONS:
            raise ContractViolation(
                f"CLAUDE.md §2: {p!r} is neither mast cells nor a declared control population."
            )


# --------------------------------------------------------------------------
# §3  Comparison hierarchy
# --------------------------------------------------------------------------
PRIMARY = "Healthy_vs_AD"
SECONDARY = ("Healthy_vs_AD_NL", "Healthy_vs_AD_LS")
REFINEMENT = "AD_LS_vs_AD_NL"

CANONICAL_ORDER = (PRIMARY,) + SECONDARY + (REFINEMENT,)


def comparison_order(comparisons):
    """§3 — enforce the reporting order; refuse to lead with lesional vs non-lesional."""
    comps = list(comparisons)
    unknown = [c for c in comps if c not in CANONICAL_ORDER]
    if unknown:
        raise ContractViolation(
            f"CLAUDE.md §3: unrecognised comparison(s) {unknown}; "
            f"expected a subset of {list(CANONICAL_ORDER)}"
        )
    if comps and comps[0] == REFINEMENT:
        raise ContractViolation(
            "CLAUDE.md §3: must not lead with AD lesional vs non-lesional. "
            "The primary comparison is Healthy vs AD; report it first even if confounded."
        )
    if PRIMARY in comps and comps[0] != PRIMARY:
        raise ContractViolation(
            f"CLAUDE.md §3: Healthy vs AD must be reported first, got {comps[0]!r}."
        )
    return sorted(comps, key=CANONICAL_ORDER.index)


# --------------------------------------------------------------------------
# §6  Mandatory controls on every claim
# --------------------------------------------------------------------------
DEPTH_TERM_RE = re.compile(r"log10[_(]?(depth|total_counts|umi|lib)", re.I)


def check_model_terms(terms, *, grouping=None) -> None:
    """§6 — every group comparison adjusts for library size and clusters by donor."""
    joined = " + ".join(map(str, terms))
    if not DEPTH_TERM_RE.search(joined):
        raise ContractViolation(
            "CLAUDE.md §6: model must adjust for library size — no log10(depth) term "
            f"in {joined!r}. TNFRSF9 detection is depth-dominated."
        )
    if grouping is None:
        raise ContractViolation(
            "CLAUDE.md §6: cell-level tests pseudo-replicate; declare the donor-level "
            "grouping (aggregate by donor, or cluster SEs on donor)."
        )
    if not re.search(r"donor|subject|patient|ind_id", str(grouping), re.I):
        raise ContractViolation(
            f"CLAUDE.md §6: grouping must be donor-level, got {grouping!r}."
        )


MAST_QC_MARKERS = ("TPSAB1", "TPSB2", "CPA3")


def marker_qc_mask(counts_by_gene, *, min_markers: int = 2):
    """§6 — mast cells must express >=2 of TPSAB1/TPSB2/CPA3; deposited labels are unreliable."""
    missing = [g for g in MAST_QC_MARKERS if g not in counts_by_gene]
    if missing:
        raise ContractViolation(
            f"CLAUDE.md §6: mast marker QC impossible, absent from matrix: {missing}"
        )
    import numpy as np
    stack = np.vstack([np.asarray(counts_by_gene[g]) > 0 for g in MAST_QC_MARKERS])
    return stack.sum(axis=0) >= min_markers


# --------------------------------------------------------------------------
# Repo audit
# --------------------------------------------------------------------------
# §2 — banned framings: control populations written up as findings.
BANNED_PATTERNS = [
    (r"\bthe signal comes from T[- ]cells\b", "§2 T cells reported as the result"),
    (r"\bT[- ]cell activation state", "§2 T-cell activation described as a finding"),
    (r"\bregress\w*\s+(tissue\s+)?signal\s+on(to)?\s+T[- ]cell", "§2 regression onto T-cell content"),
    (r"\bwe (identify|describe|report) a novel (T[- ]cell|macrophage)", "§2 control population as finding"),
]

# §8 — bare P values.
BARE_P_RE = re.compile(r"\bP\s*[<=]\s*0?\.\d+")

COMPARISON_HINT = re.compile(
    r"healthy.{0,12}vs.{0,12}ad|ad.{0,12}vs.{0,12}healthy|lesional.{0,20}vs.{0,20}non[- ]lesional",
    re.I,
)


def _audit_scripts(problems: list[str]) -> None:
    """Every script making a group comparison must use the guardrail API (§9).

    Two tiers, because a plotting script and a modelling script owe different
    things: anything that *orders* comparisons must validate that order, but only
    something that actually *fits* a model can declare its terms.
    """
    for path in sorted(ANALYSIS.glob("s*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        makes_comparison = bool(COMPARISON_HINT.search(text)) or "GROUPS" in text
        if not makes_comparison:
            continue
        if "guardrails" not in text:
            problems.append(f"{path.name}: makes a group comparison but does not import guardrails (§9)")
            continue
        required = ["require_mast_subject", "comparison_order"]
        fits_model = bool(re.search(r"sm\.GLM|statsmodels|\.fit\(", text))
        if fits_model:
            required.append("check_model_terms")
        for fn in required:
            if fn not in text:
                problems.append(f"{path.name}: does not call {fn}() (§9)")


def _audit_figures(problems: list[str]) -> None:
    """Figure conventions that have regressed repeatedly.

    1. Panel titles must come from palette.panel_label() — bold letter, plain
       description. Passing fontweight="bold" to set_title bolds the description
       too and produces inconsistent labelling across panels.
    2. §2 — control populations must not be drawn as co-equal series in a MAIN
       figure. They belong in the supplementary controls figure, the tables and
       the report text. The subject of every main panel is mast cells.
    """
    CONTROL_NAMES = ("Fibroblast", "Keratinocyte", "T/NK", "Macrophage")
    for path in sorted(ANALYSIS.glob("s3*.py")):
        text = path.read_text(encoding="utf-8", errors="replace")
        if re.search(r'set_title\([^)]*fontweight\s*=\s*["\']bold', text):
            problems.append(
                f"{path.name}: set_title(..., fontweight='bold') — use "
                f"palette.panel_label() so only the panel letter is bold")
        if "figS" in path.name:      # the controls figure may plot controls
            continue
        for name in CONTROL_NAMES:
            if re.search(rf'label\s*=\s*["\'][^"\']*{re.escape(name)}', text):
                problems.append(
                    f"{path.name}: plots {name!r} as a labelled series in a main "
                    f"figure (§2: control populations are not co-equal findings)")


def _audit_prose(problems: list[str]) -> None:
    """§2 and §8 checks over written output."""
    targets = list((ROOT / "results").rglob("*.md")) + list((ROOT / "docs").rglob("*.md"))
    targets += [p for p in ROOT.glob("*.md") if p.name != "CLAUDE.md"]
    for path in targets:
        text = path.read_text(encoding="utf-8", errors="replace")
        for pattern, why in BANNED_PATTERNS:
            if re.search(pattern, text, re.I):
                problems.append(f"{path.relative_to(ROOT)}: banned framing — {why}")
        # §8: a P value must travel with an effect size and an n somewhere in the same block.
        for block in re.split(r"\n\s*\n", text):
            if BARE_P_RE.search(block):
                # accept Unicode subscripts (log₂FC) as well as ASCII
                # A group median or rate reported per arm IS the effect size for a
                # rank test, so it counts; a lone P value still does not.
                has_effect = re.search(
                    r"log[2₂]?\s?FC|logFC|fold|OR\s*=|odds ratio|β|beta|rate ratio|"
                    r"RR\s*=|estimate|%|prevalence|Cliff|d\s*=|median|rho|per 10",
                    block, re.I)
                has_n = re.search(r"\bn\s*=|\bN\s*=|donors?\b|cells\b|samples?\b|spots\b",
                                  block, re.I)
                if not (has_effect and has_n):
                    snippet = " ".join(block.split())[:90]
                    problems.append(
                        f"{path.relative_to(ROOT)}: bare P value without effect size + n (§8): {snippet!r}")


def _audit_ordering(problems: list[str]) -> None:
    """§3 — the report must not lead with the lesional/non-lesional refinement."""
    for path in (ROOT / "results").rglob("*.md"):
        text = path.read_text(encoding="utf-8", errors="replace")
        low = text.lower()
        i_ls_nl = low.find("lesional versus non-lesional")
        if i_ls_nl == -1:
            i_ls_nl = low.find("lesional vs non-lesional")
        i_primary = min([i for i in (low.find("healthy versus ad"), low.find("healthy vs ad")) if i != -1] or [-1])
        if i_ls_nl != -1 and i_primary != -1 and i_ls_nl < i_primary:
            problems.append(
                f"{path.relative_to(ROOT)}: leads with lesional vs non-lesional before Healthy vs AD (§3)")


def _audit_selfconsistency(problems: list[str]) -> None:
    """Contract file must be present and the guardrail API intact (§9)."""
    if not (ROOT / "CLAUDE.md").exists():
        problems.append("CLAUDE.md missing from repo root")
    if not (ANALYSIS / "palette.py").exists():
        problems.append("analysis/palette.py missing (§8 figure standard)")


def main() -> int:
    problems: list[str] = []
    _audit_selfconsistency(problems)
    _audit_scripts(problems)
    _audit_figures(problems)
    _audit_prose(problems)
    _audit_ordering(problems)

    if problems:
        print("CONTRACT VIOLATIONS", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        print(f"\n{len(problems)} violation(s).", file=sys.stderr)
        return 1
    print("guardrails: OK — no CLAUDE.md violations detected.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
