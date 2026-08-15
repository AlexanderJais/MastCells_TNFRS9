"""Figure 4 — whole-skin bulk RNA-seq (GSE121212, 147 samples), with psoriasis
as a disease-specificity control.

This arm exists for three things the single-cell cohorts cannot deliver, all of
which need a large, independent, disease-controlled tissue cohort:

  1. a well-powered dose-response for whole-skin TNFRSF9
     (healthy -> AD non-lesional -> AD lesional);
  2. an independent test of the single-cell abundance result — mast-cell content
     stays flat while TNFRSF9 rises 4-fold;
  3. the specificity check — psoriasis raises whole-skin TNFRSF9 as much as AD
     does, so the tissue-level signal marks inflamed skin, not AD.

The previous version carried a 147-point five-colour TNFRSF9-vs-mast-content
scatter. It was unreadable and the relationship it was meant to show is a single
number: TNFRSF9 does not track mast-cell content within any arm (Spearman rho
+0.12 to +0.26, all P > 0.18). That number is annotated on panel b instead.

CLAUDE.md §4 caveat: in bulk tissue, abundance and per-cell expression are not
separable, and the mast-cell contribution is not separable from any other source.
Error bars are +/- 1 SEM.
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

ORDER = ["Healthy", "AD_NL", "AD_LS", "PSO_NL", "PSO_LS"]
LAB = {"Healthy": "Healthy", "AD_NL": "AD NL", "AD_LS": "AD LS",
       "PSO_NL": "Pso NL", "PSO_LS": "Pso LS"}
COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"],
       "PSO_NL": "#C0B6C7", "PSO_LS": "#8A7C96"}


def strip(ax, s, col, ylab, letter, title, note=None):
    rng = np.random.default_rng(7)
    for i, a in enumerate(ORDER):
        v = s.loc[s.arm == a, col].astype(float).values
        v = v[np.isfinite(v) & (v > 0)]
        ax.scatter(i + rng.uniform(-.17, .17, len(v)), v, s=9, lw=.28,
                   facecolor=COL[a], edgecolor="#333333", alpha=.9, zorder=3)
        ax.plot([i - .3, i + .3], [np.median(v)] * 2, color="#1F1F1F",
                lw=1.3, zorder=4, solid_capstyle="butt")
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([LAB[a] for a in ORDER], fontsize=5.8, rotation=30,
                       ha="right")
    ax.set_ylabel(ylab, fontsize=6.4)
    ax.set_yscale("log")
    ax.set_title(rf"$\bf{{{letter}}}$   {title}", loc="left", fontsize=6.8, pad=4)
    if note:
        ax.text(.02, .03, note, transform=ax.transAxes, fontsize=5.2,
                color="#4A4A4A", va="bottom")


def main() -> int:
    set_style()
    s = pd.read_csv(TAB / "bulk_GSE121212_per_sample.csv", index_col=0)
    R = pd.read_csv(TAB / "bulk_GSE121212_glm.csv")

    fig = plt.figure(figsize=(DOUBLE_COL * 0.82, DOUBLE_COL * 0.36))
    gs = fig.add_gridspec(1, 3, wspace=0.58, width_ratios=[1, 1, 1.25])

    # ---- a. the dose-response, and psoriasis alongside it -----------------
    strip(fig.add_subplot(gs[0, 0]), s, "TNFRSF9_cpm", "TNFRSF9 (CPM)",
          "a", "whole-skin TNFRSF9")

    # ---- b. mast content over the same samples ----------------------------
    rho_txt = []
    for a in ("Healthy", "AD_LS"):
        g = s[s.arm == a]
        r, _ = stats.spearmanr(g.TNFRSF9_cpm, g.mast_content_cpm)
        rho_txt.append(f"{LAB[a]} ρ={r:+.2f}")
    strip(fig.add_subplot(gs[0, 1]), s, "mast_content_cpm",
          "TPSAB1+TPSB2+CPA3 (CPM)", "b", "mast-cell content",
          note="TNFRSF9 vs mast content:\n" + ", ".join(rho_txt))

    # ---- c. the two diseases side by side, target vs mast content ---------
    ax = fig.add_subplot(gs[0, 2])
    rows = [("Healthy_vs_AD", TARGET, "AD"), ("Healthy_vs_PSO", TARGET, "Psoriasis")]
    mast_rows = [("Healthy_vs_AD", "TPSAB1", "AD"), ("Healthy_vs_PSO", "TPSAB1", "Psoriasis")]
    ys = np.arange(2)
    w = 0.34
    for k, (grp, colr, lbl) in enumerate([(rows, MUTED["accent"], "TNFRSF9"),
                                          (mast_rows, MUTED["mast"], "TPSAB1 (mast content)")]):
        off = (k - 0.5) * w
        for i, (comp, gene, _) in enumerate(grp):
            r = R[(R.comparison == comp) & (R.gene == gene)].iloc[0]
            ax.barh(i + off, r.log2FC, height=w * .9, color=colr, edgecolor="none",
                    label=lbl if i == 0 else None)
            ax.errorbar(r.log2FC, i + off, xerr=r.se_log2, fmt="none",
                        ecolor="#4A4A4A", elinewidth=.8, capsize=1.6)
    ax.axvline(0, color="#B0B0B0", lw=.6, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels(["AD\n(n = 38 vs 54)", "Psoriasis\n(n = 38 vs 55)"], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change vs healthy  (± SEM)", fontsize=6.2)
    ax.legend(fontsize=5.4, loc="upper left", frameon=False,
              handlelength=1.2, handletextpad=.4)
    ax.set_title(r"$\bf{c}$   target up, mast content flat — in both diseases",
                 loc="left", fontsize=6.8, pad=4)

    fig.savefig(FIG / "fig4_bulk.pdf")
    fig.savefig(FIG / "fig4_bulk.png")

    print("median CPM by arm:")
    print(s.groupby("arm")[["TNFRSF9_cpm", "mast_content_cpm"]].median().round(2).to_string())
    print("\nSpearman(TNFRSF9, mast content) within arm:")
    for a, g in s.groupby("arm"):
        r, p = stats.spearmanr(g.TNFRSF9_cpm, g.mast_content_cpm)
        print(f"  {a:9s} n={len(g):3d}  rho={r:+.3f}  P={p:.3f}")
    print("\nwrote fig4_bulk.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
