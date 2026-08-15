"""Figure 6 — TNFRSF9 is an activation-induced gene in primary human mast cells.

Subject: MAST CELLS (CLAUDE.md §2). This is the controlled experiment the tissue
data cannot be: resting versus stimulated primary mast cells, no dissociation,
no ambient RNA, no doublets, no depth confound.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402

import guardrails as G  # noqa: E402
from genes import TARGET, TARGET_LIGAND  # noqa: E402
from palette import DOUBLE_COL, MUTED, panel_label, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")

ORDER = ["UN", "TSLP", "IL25", "IgE", "IL33", "IL33IL25TSLP",
         "IgETSLP", "IgEIL25", "IgEIL33", "IgEIL33IL25TSLP"]
NICE = {"UN": "resting", "TSLP": "TSLP", "IL25": "IL-25", "IgE": "IgE/Ag",
        "IL33": "IL-33", "IL33IL25TSLP": "IL-33+IL-25\n+TSLP",
        "IgETSLP": "IgE+TSLP", "IgEIL25": "IgE+IL-25", "IgEIL33": "IgE+IL-33",
        "IgEIL33IL25TSLP": "IgE+all three"}
TEFF = ["Ctl", "LT", "LT_B", "IgE", "IL33"]
TEFF_NICE = {"Ctl": "resting", "LT": "resting\nT cells", "LT_B": "activated\nT cells",
             "IgE": "IgE/Ag", "IL33": "IL-33"}


def main() -> int:
    set_style()
    R = pd.read_csv(TAB / "invitro_GSE196862_fpkm.csv", index_col=0)
    F = pd.read_csv(TAB / "invitro_GSE196862_log2fc.csv", index_col=0)
    T = pd.read_csv(TAB / "invitro_GSE235240.csv", index_col=0)

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.46))
    gs = fig.add_gridspec(1, 3, wspace=0.62, width_ratios=[1.5, 1.15, 0.95])

    # ---- a. TNFRSF9 FPKM by condition ------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    order = [o for o in ORDER if o in R.columns]
    vals = R.loc[TARGET, order].astype(float).values
    colors = [MUTED["healthy"] if o == "UN" else
              ("#C7C7C7" if o in ("TSLP", "IL25") else MUTED["ad_ls"]) for o in order]
    ax.bar(np.arange(len(order)), vals, color=colors, edgecolor="none", width=.72)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([NICE[o] for o in order], rotation=45, ha="right", fontsize=5.2)
    ax.set_ylabel(f"{TARGET} (FPKM)", fontsize=6.2)
    panel_label(ax, "a", "primary human skin mast cells")
    ax.text(.02, .95, "resting FPKM 0.14\nIgE+IL-33 → 8.7  (59×)", transform=ax.transAxes,
            fontsize=5.4, va="top", color="#4A4A4A")

    # ---- b. specificity: log2FC heat ------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    rows = [TARGET, TARGET_LIGAND, "TNFRSF18", "TNFRSF4", "IL13", "TNF",
            "TPSAB1", "TPSB2", "B2M", "RPL13A"]
    rows = [r for r in rows if r in F.index]
    cols = [c for c in order if c != "UN"]
    M = F.loc[rows, cols].astype(float).values
    im = ax.imshow(M, cmap="RdBu_r", vmin=-6, vmax=6, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([NICE[c] for c in cols], rotation=45, ha="right", fontsize=4.8)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([r if r != TARGET else f"$\\bf{{{r}}}$" for r in rows], fontsize=5.4)
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.02)
    cb.set_label("log$_2$ FC vs resting", fontsize=5.0, labelpad=1)
    cb.ax.tick_params(labelsize=4.6)
    panel_label(ax, "b", "specificity")

    # ---- c. paired donors, second cohort --------------------------------
    ax = fig.add_subplot(gs[0, 2])
    conds = [c for c in TEFF if c in T.columns]
    y = [T.loc[TARGET, c] for c in conds]
    ax.bar(np.arange(len(conds)), y,
           color=[MUTED["healthy"] if c == "Ctl" else
                  ("#C7C7C7" if c == "LT" else MUTED["ad_ls"]) for c in conds],
           edgecolor="none", width=.7)
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(conds)))
    ax.set_xticklabels([TEFF_NICE[c] for c in conds], rotation=45, ha="right", fontsize=5.2)
    ax.set_ylabel(f"{TARGET} (normalised counts)", fontsize=6.0, labelpad=1)
    panel_label(ax, "c", "4 donors, paired")
    for i, c in enumerate(conds):
        p = T.get(f"p_{c}")
        if p is not None and np.isfinite(T.loc[TARGET, f"p_{c}"]):
            ax.text(i, y[i] * 1.35, f"P={T.loc[TARGET, f'p_{c}']:.3f}",
                    ha="center", fontsize=4.4, color="#4A4A4A")

    fig.savefig(FIG / "fig6_invitro.pdf")
    fig.savefig(FIG / "fig6_invitro.png")
    print("wrote fig6_invitro.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
