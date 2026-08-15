"""Figure 2 — quantification of Figure 1d, with statistics.

Four numbers, per arm, per sample, nothing else:

  a  mast cells recovered per sample        (absolute)
  b  mast cells as a percentage of all cells (abundance, capture-independent)
  c  TNFRSF9+ mast cells per sample          (absolute)
  d  TNFRSF9+ mast cells as a percentage of mast cells   <- Figure 1d, quantified

Each point is one sample (one donor x arm). Statistics are donor-level
Mann-Whitney against the healthy arm, plus a pooled Fisher exact test on the raw
TNFRSF9+ cell counts for panel d — both reported on the figure.

No ratios of ratios, no products: the CLAUDE.md §4 decomposition is recorded in
the tables and in Section 2.1 of the report, but the figure states the measured
quantities directly.

Subject: MAST CELLS (CLAUDE.md §2).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, panel_label, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY, *G.SECONDARY, G.REFINEMENT])

ORDER = ["Healthy", "AD_NL", "AD_LS"]
LAB = {"Healthy": "Healthy", "AD_NL": "AD NL", "AD_LS": "AD LS"}
COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"]}


def pstar(p: float) -> str:
    if not np.isfinite(p):
        return "n.s."
    return f"P = {p:.3f}" if p >= 0.001 else "P < 0.001"


def strip_panel(ax, s: pd.DataFrame, col: str, ylab: str, letter: str,
                title: str, tests: dict | None = None, pooled: str | None = None):
    rng = np.random.default_rng(11)
    ymax = 0.0
    for i, a in enumerate(ORDER):
        v = s.loc[s.arm == a, col].astype(float).values
        v = v[np.isfinite(v)]
        ax.scatter(i + rng.uniform(-.14, .14, len(v)), v, s=16, lw=.35,
                   facecolor=COL[a], edgecolor="#333333", alpha=.95, zorder=3)
        if len(v):
            ax.plot([i - .27, i + .27], [np.median(v)] * 2, color="#1F1F1F",
                    lw=1.3, zorder=4, solid_capstyle="butt")
            ymax = max(ymax, np.nanmax(v))
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LAB[a] for a in ORDER], fontsize=6.4)
    ax.set_ylabel(ylab, fontsize=6.4)
    ax.set_title(rf"$\bf{{{letter}}}$   {title}", loc="left", fontsize=6.8, pad=4)
    ax.set_ylim(bottom=-0.04 * max(ymax, 1e-9))

    # significance brackets against the healthy arm
    if tests:
        span = ymax if ymax > 0 else 1.0
        for j, (arm, p) in enumerate(tests.items()):
            i = ORDER.index(arm)
            y = span * (1.10 + 0.16 * j)
            ax.plot([0, 0, i, i], [y - span * .03, y, y, y - span * .03],
                    lw=.7, color="#5A5A5A")
            ax.text((0 + i) / 2, y + span * .015, pstar(p), ha="center",
                    va="bottom", fontsize=5.6, color="#2B2B2B")
        ax.set_ylim(top=span * (1.10 + 0.16 * len(tests) + 0.13))
    if pooled:
        ax.text(.5, -.30, pooled, transform=ax.transAxes, ha="center",
                fontsize=5.6, color="#4A4A4A")


def main() -> int:
    set_style()
    s = pd.read_csv(TAB / "sc_sample_level.csv")
    s["pct_mast"] = 100 * s.n_mast / s.n_cells
    s["pct_t9pos"] = 100 * s.n_mast_t9pos / s.n_mast.replace(0, np.nan)

    def mwu(col: str, arm: str) -> float:
        h = s.loc[s.arm == "Healthy", col].astype(float).dropna()
        a = s.loc[s.arm == arm, col].astype(float).dropna()
        if len(h) < 2 or len(a) < 2:
            return np.nan
        return float(stats.mannwhitneyu(a, h, alternative="two-sided").pvalue)

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.40))
    gs = fig.add_gridspec(1, 4, wspace=0.52)

    # ---- a. absolute mast cells per sample -------------------------------
    strip_panel(fig.add_subplot(gs[0, 0]), s, "n_mast",
                "mast cells per sample (n)", "a", "mast cells recovered",
                tests={"AD_NL": mwu("n_mast", "AD_NL"), "AD_LS": mwu("n_mast", "AD_LS")},
                pooled=f"total {int(s.n_mast.sum()):,} mast cells")

    # ---- b. mast cells as % of all cells ---------------------------------
    strip_panel(fig.add_subplot(gs[0, 1]), s, "pct_mast",
                "mast cells (% of all cells)", "b", "mast-cell abundance",
                tests={"AD_NL": mwu("pct_mast", "AD_NL"), "AD_LS": mwu("pct_mast", "AD_LS")},
                pooled=None)

    # ---- c. absolute TNFRSF9+ mast cells ---------------------------------
    strip_panel(fig.add_subplot(gs[0, 2]), s, "n_mast_t9pos",
                f"{TARGET}$^+$ mast cells per sample (n)", "c",
                f"{TARGET}$^+$ mast cells",
                tests={"AD_NL": mwu("n_mast_t9pos", "AD_NL"),
                       "AD_LS": mwu("n_mast_t9pos", "AD_LS")},
                pooled=f"total {int(s.n_mast_t9pos.sum())} {TARGET}$^+$ cells")

    # ---- d. % of mast cells that are TNFRSF9+  (= Fig 1d, quantified) ----
    ax = fig.add_subplot(gs[0, 3])
    # pooled Fisher exact, healthy vs AD (lesional + non-lesional)
    h_pos = int(s.loc[s.arm == "Healthy", "n_mast_t9pos"].sum())
    h_tot = int(s.loc[s.arm == "Healthy", "n_mast"].sum())
    a_pos = int(s.loc[s.arm != "Healthy", "n_mast_t9pos"].sum())
    a_tot = int(s.loc[s.arm != "Healthy", "n_mast"].sum())
    orr, pf = stats.fisher_exact([[a_pos, a_tot - a_pos], [h_pos, h_tot - h_pos]])
    strip_panel(ax, s, "pct_t9pos", f"{TARGET}$^+$ (% of mast cells)", "d",
                f"{TARGET}$^+$ fraction",
                tests={"AD_NL": mwu("pct_t9pos", "AD_NL"),
                       "AD_LS": mwu("pct_t9pos", "AD_LS")},
                pooled=(f"pooled {h_pos}/{h_tot} vs {a_pos}/{a_tot} cells\n"
                        f"Fisher OR {orr:.1f}, {pstar(pf)}"))

    fig.savefig(FIG / "fig2_quantification.pdf")
    fig.savefig(FIG / "fig2_quantification.png")

    # ---- the same numbers as a table -------------------------------------
    out = s.groupby("arm").apply(lambda g: pd.Series({
        "samples": len(g), "donors": g.donor.nunique(),
        "cells_total": int(g.n_cells.sum()),
        "mast_total": int(g.n_mast.sum()),
        "mast_per_sample_median": g.n_mast.median(),
        "mast_pct_median": g.pct_mast.median(),
        "t9pos_total": int(g.n_mast_t9pos.sum()),
        "t9pos_per_sample_median": g.n_mast_t9pos.median(),
        "t9pos_pct_pooled": 100 * g.n_mast_t9pos.sum() / g.n_mast.sum(),
        "t9pos_pct_median": g.pct_t9pos.median()})).reindex(ORDER)
    out.to_csv(TAB / "fig2_quantification.csv")
    print(out.round(3).to_string())
    print(f"\npooled Fisher (Healthy vs AD): {h_pos}/{h_tot} vs {a_pos}/{a_tot}, "
          f"OR {orr:.2f}, P = {pf:.2e}")
    for c in ("n_mast", "pct_mast", "n_mast_t9pos", "pct_t9pos"):
        print(f"  MWU vs healthy — {c:14s} AD_NL P={mwu(c,'AD_NL'):.4f}  "
              f"AD_LS P={mwu(c,'AD_LS'):.4f}")
    print("\nwrote fig2_quantification.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
