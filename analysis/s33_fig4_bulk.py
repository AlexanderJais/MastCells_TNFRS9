"""Figure 4 — whole-skin bulk RNA-seq (GSE121212), with psoriasis as a
disease-specificity control.

CLAUDE.md §4 is the governing caveat: in bulk tissue mast-cell abundance and
per-mast-cell expression cannot be separated at all, and neither can the
mast-cell contribution be separated from any other source. What bulk can say is
whether whole-skin TNFRSF9 moves, and whether mast-cell content moves with it.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY, *G.SECONDARY, G.REFINEMENT])
ORDER = ["Healthy", "AD_NL", "AD_LS", "PSO_NL", "PSO_LS"]
LAB = {"Healthy": "Healthy", "AD_NL": "AD NL", "AD_LS": "AD LS",
       "PSO_NL": "Pso NL", "PSO_LS": "Pso LS"}
COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"],
       "PSO_NL": "#B9B0C2", "PSO_LS": "#8E8098"}


def strip(ax, s, col, ylab, title, logy=True):
    rng = np.random.default_rng(7)
    for i, a in enumerate(ORDER):
        v = s.loc[s.arm == a, col].astype(float).values
        v = v[np.isfinite(v)]
        ax.scatter(i + rng.uniform(-.16, .16, len(v)), v, s=8, lw=.25,
                   facecolor=COL[a], edgecolor="#3A3A3A", alpha=.85, zorder=3)
        ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="#2B2B2B", lw=1.1, zorder=4)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LAB[a] for a in ORDER], fontsize=5.8, rotation=30)
    ax.set_ylabel(ylab, fontsize=6)
    ax.set_title(title, loc="left", fontweight="bold")
    if logy:
        ax.set_yscale("log")


def main() -> int:
    set_style()
    s = pd.read_csv(TAB / "bulk_GSE121212_per_sample.csv", index_col=0)
    R = pd.read_csv(TAB / "bulk_GSE121212_glm.csv")

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.42))
    gs = fig.add_gridspec(1, 4, wspace=0.62)

    strip(fig.add_subplot(gs[0, 0]), s, "TNFRSF9_cpm", "TNFRSF9 (CPM)",
          "a   whole-skin TNFRSF9")
    strip(fig.add_subplot(gs[0, 1]), s, "mast_content_cpm",
          "TPSAB1+TPSB2+CPA3 (CPM)", "b   mast-cell content")

    # ---- c. the two together --------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    for a in ORDER:
        d = s[s.arm == a]
        ax.scatter(d.mast_content_cpm, d.TNFRSF9_cpm, s=9, lw=.25,
                   facecolor=COL[a], edgecolor="#3A3A3A", alpha=.85, label=LAB[a])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mast content (CPM)", fontsize=6)
    ax.set_ylabel("TNFRSF9 (CPM)", fontsize=6)
    ax.legend(fontsize=4.8, loc="upper left", ncol=2, handletextpad=.1,
              columnspacing=.6, borderpad=.15)
    ax.set_title("c   §4 not separable in bulk", loc="left", fontweight="bold")

    # ---- d. forest -------------------------------------------------------
    ax = fig.add_subplot(gs[0, 3])
    comps = [(G.PRIMARY, TARGET, "TNFRSF9 · AD"),
             ("Healthy_vs_AD_LS", TARGET, "TNFRSF9 · AD LS"),
             ("Healthy_vs_PSO", TARGET, "TNFRSF9 · psoriasis"),
             (G.PRIMARY, "TPSAB1", "TPSAB1 · AD"),
             (G.PRIMARY, "TPSB2", "TPSB2 · AD"),
             (G.PRIMARY, "COL1A1", "COL1A1 · AD"),
             (G.PRIMARY, "KRT14", "KRT14 · AD")]
    for i, (c, g, lab) in enumerate(comps):
        r = R[(R.comparison == c) & (R.gene == g)].iloc[0]
        col = MUTED["accent"] if g == TARGET else ("#7A5C7B" if g in ("TPSAB1", "TPSB2") else "#9A9A9A")
        ax.plot([r.log2FC - r.se_log2, r.log2FC + r.se_log2], [i, i],
                color="#6B6B6B", lw=1.1)   # +/- 1 SEM
        ax.plot(r.log2FC, i, "o", ms=4, color=col, zorder=5)
    ax.axvline(0, color="#B0B0B0", lw=.6, ls="--")
    ax.set_yticks(range(len(comps)))
    ax.set_yticklabels([c[2] for c in comps], fontsize=5.2)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change  (± SEM)", fontsize=6)
    ax.set_title("d   vs healthy, library-size-adjusted", loc="left", fontweight="bold")

    fig.savefig(FIG / "fig4_bulk.pdf")
    fig.savefig(FIG / "fig4_bulk.png")
    print("wrote fig4_bulk.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
