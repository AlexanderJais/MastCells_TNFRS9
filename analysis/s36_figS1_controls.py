"""Supplementary Figure 1 — the CLAUDE.md §6 controls and the calibration.

These are the checks that decide how much the main result is worth: the mandated
control populations, the empirical calibration against abundance-matched genes,
and the dissociation-stress burden per arm. The measured quantities themselves
are quantified directly in Figure 2 (s32_fig2_quantification.py).
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
ARM_COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"]}
ARM_LAB = {"Healthy": "Healthy", "AD_NL": "AD NL", "AD_LS": "AD LS"}
ORDER = ["Healthy", "AD_NL", "AD_LS"]


def strip(ax, s: pd.DataFrame, col: str, ylab: str, title: str, logy=False):
    rng = np.random.default_rng(3)
    for i, a in enumerate(ORDER):
        v = s.loc[s.arm == a, col].astype(float).values
        v = v[np.isfinite(v)]
        x = i + rng.uniform(-.13, .13, len(v))
        ax.scatter(x, v, s=11, lw=.3, facecolor=ARM_COL[a], edgecolor="#3A3A3A",
                   alpha=.9, zorder=3)
        if len(v):
            ax.plot([i - .28, i + .28], [np.median(v)] * 2, color="#2B2B2B", lw=1.1, zorder=4)
    ax.set_xticks(range(len(ORDER)))
    ax.set_xticklabels([ARM_LAB[a] for a in ORDER], fontsize=6)
    ax.set_ylabel(ylab, fontsize=6)
    ax.set_title(title, loc="left", fontweight="bold")
    if logy:
        ax.set_yscale("symlog", linthresh=1e-2)


def main() -> int:
    set_style()
    s = pd.read_csv(TAB / "sc_sample_level.csv")
    # §4 decomposition in per-CELL units: the identity F1 x F2 = product holds
    # exactly, and per-cell is the report's primary estimand (Section 7.1).
    s["f2_per100"] = 100 * s.t9_mast / s.n_mast.replace(0, np.nan)
    s["prod_per1k"] = 1000 * s.t9_mast / s.n_cells

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.34))
    gs = fig.add_gridspec(1, 3, wspace=0.46)

    # ---- §6 control populations ------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    C = pd.read_csv(TAB / "sc_control_populations.csv")
    C = C[C.comparison == G.PRIMARY].set_index("population")
    pops = ["Mast (SUBJECT)", "Fibroblasts", "Keratinocytes", "T/NK", "Macrophages"]
    for i, p in enumerate(pops):
        r = C.loc[p]
        col = MUTED["mast"] if "SUBJECT" in p else "#9A9A9A"
        ax.plot([r.log2FC - r.se_log2, r.log2FC + r.se_log2], [i, i],
                color="#6B6B6B", lw=1.1)   # +/- 1 SEM
        ax.plot(r.log2FC, i, "o", ms=4.5, color=col, zorder=5)
    ax.axvline(0, color="#B0B0B0", lw=.6, ls="--")
    ax.set_yticks(range(len(pops)))
    ax.set_yticklabels(["Mast (subject)", "Fibroblasts", "Keratinocytes",
                        "T/NK", "Macrophages"], fontsize=5.8)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ FC per UMI (adjusted), Healthy vs AD  (± SEM)", fontsize=6)
    ax.set_title("a   §6 control populations", loc="left", fontweight="bold")

    # ---- calibration against abundance-matched genes ---------------------
    ax = fig.add_subplot(gs[0, 1])
    N = pd.read_csv(TAB / "sc_calibration_matched_genes.csv")
    ax.hist(N.log2FC.clip(-6, 6), bins=60, color="#D5D5D5", edgecolor="none")
    t9 = N.loc[N.gene == TARGET, "log2FC"].iat[0]
    ax.axvline(t9, color=MUTED["accent"], lw=1.2)
    ax.text(t9, ax.get_ylim()[1] * .92, " TNFRSF9", fontsize=6,
            color=MUTED["accent"], ha="left", va="top")
    ax.set_xlabel("log$_2$ FC (Healthy vs AD), mast cells", fontsize=6)
    ax.set_ylabel(f"abundance-matched genes (n={len(N)})", fontsize=5.8)
    ax.set_title("b   TNFRSF9 vs the background shift", loc="left", fontweight="bold")

    # ---- §6 dissociation stress ------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    strip(ax, s, "stress_mast", "HSP+IEG per 10k mast UMI",
          "c   §6 dissociation stress")
    ax.text(.02, .96, "healthy arm is the most stressed:\na stressed arm reads low",
            transform=ax.transAxes, fontsize=5.2, va="top", color="#5A5A5A")

    fig.savefig(FIG / "figS1_controls.pdf")
    fig.savefig(FIG / "figS1_controls.png")
    print("wrote figS1_controls.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
