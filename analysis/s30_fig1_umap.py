"""Figure 1 — the mast-cell compartment and TNFRSF9 in the discovery atlas.

Subject: MAST CELLS (CLAUDE.md §2). Other populations are drawn only to place
the mast-cell island in context.

Panels are restricted to the cell map and the TNFRSF9 detection it supports.
Quantities that depend on per-cell UMI are NOT shown here: per-cell UMI is
sampling effort (library reads x capture x the cell's own mRNA content), not
sequencing depth, and the analysis establishing that lives in
s08_depth_decomposition.py and Section 7.1 of the report. Per-donor rates are in
Fig. 2a-c.
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
import matplotlib.patheffects as pe  # noqa: E402

import guardrails as G  # noqa: E402
from genes import MAST_QC, TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")
ARMS = ["Healthy", "AD_NL", "AD_LS"]
ARM_LAB = {"Healthy": "Healthy", "AD_NL": "AD non-lesional", "AD_LS": "AD lesional"}
ARM_COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"]}

CT_COL = {
    "Keratinocytes": "#C9BBA8", "Cornified keratinocytes": "#DCD2C4",
    "Fibroblasts": "#A8968A", "Pericyte/SMC": "#B3A7B8", "VEC": "#93A7C4",
    "LEC": "#A9BFD4", "Melanocytes": "#C3B0C6", "Schwann": "#BFC8B4",
    "Sweat gland": "#D0C6B0", "DC": "#8FA9AE", "Macrophages": "#7FA1A8",
    "T/NK": "#8FA36B", "B cells": "#A7B98C", "Plasma": "#B6C79E",
    "Neutrophils": "#C7C39A", "Mast": MUTED["mast"],
}
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]


def panel_label(ax, letter: str, text: str = "") -> None:
    """Bold panel letter, plain descriptive text — identical treatment everywhere.

    Mathtext keeps the letter bold and the description plain inside one title
    string, so the gap between them cannot collapse on a narrow panel.
    """
    ax.set_title(rf"$\bf{{{letter}}}$   {text}" if text else rf"$\bf{{{letter}}}$",
                 loc="left", fontsize=6.6, fontweight="normal", pad=4)


def bare(ax) -> None:
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")          # UMAPs must not be stretched
    for s in ax.spines.values():
        s.set_visible(False)


def load() -> pd.DataFrame:
    o = pd.read_parquet(PROC / "obs_all.parquet")
    u = pd.read_parquet(PROC / "umap_all.parquet")
    d = o.join(u[["umap1", "umap2"]], how="inner")
    return d[d.arm.isin(ARMS)].copy()


def main() -> int:
    set_style()
    d = load()
    d["ct"] = d["Cell type"].astype(str)
    d["mast_score"] = d[[f"g_{g}" for g in MAST_QC]].sum(axis=1)
    d["t9pos"] = d[f"g_{TARGET}"] > 0

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.60))
    gs = fig.add_gridspec(2, 12, height_ratios=[1.0, 0.95],
                          hspace=0.34, wspace=0.60)

    xlim = np.percentile(d.umap1, [0.02, 99.98])
    ylim = np.percentile(d.umap2, [0.02, 99.98])

    # ---------------- a. atlas by cell type -------------------------------
    ax = fig.add_subplot(gs[0, 0:6])
    rest = d[d.ct != "Mast"]
    ax.scatter(rest.umap1, rest.umap2, s=.6, lw=0, rasterized=True,
               c=[CT_COL.get(c, "#CCCCCC") for c in rest.ct], alpha=.55)
    mm = d[d.ct == "Mast"]
    ax.scatter(mm.umap1, mm.umap2, s=3.0, lw=0, rasterized=True,
               c=MUTED["mast"], alpha=.95)
    for c in CT_COL:
        sub = d[d.ct == c]
        if len(sub) < 200:
            continue
        # The mast label must never be drawn in the mast colour on the mast
        # island — it becomes invisible. Dark text, white halo, always on top.
        ax.text(sub.umap1.median(), sub.umap2.median(), c,
                fontsize=6.2 if c == "Mast" else 5.4, ha="center", va="center",
                fontweight="bold" if c == "Mast" else "normal",
                color="#2B2B2B", zorder=10, path_effects=HALO)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    bare(ax)
    ax.set_xlabel("UMAP 1", fontsize=6, labelpad=1)
    ax.set_ylabel("UMAP 2", fontsize=6, labelpad=1)
    panel_label(ax, "a", f"{len(d):,} cells · {d.donor.nunique()} donors · "
                         f"{d.sample_label.nunique()} samples")

    # ---------------- b. mast marker score --------------------------------
    ax = fig.add_subplot(gs[0, 6:9])
    o = d.sort_values("mast_score")
    sc = ax.scatter(o.umap1, o.umap2, s=.6, lw=0, rasterized=True,
                    c=np.log1p(o.mast_score), cmap="BuPu", vmin=0, vmax=3)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    bare(ax)
    panel_label(ax, "b", "tryptase / CPA3")
    cb = fig.colorbar(sc, ax=ax, fraction=.045, pad=.03)
    cb.set_label("log(1+counts)", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)

    # ---------------- c. mast island, by arm ------------------------------
    mast = d[d.mast_strict]
    x0, x1 = mast.umap1.quantile([.005, .995])
    y0, y1 = mast.umap2.quantile([.005, .995])
    pad = 1.2
    near = d[(d.umap1.between(x0 - pad, x1 + pad)) & (d.umap2.between(y0 - pad, y1 + pad))]

    ax = fig.add_subplot(gs[0, 9:12])
    ax.scatter(near.umap1, near.umap2, s=1.2, lw=0, c="#E2E2E2", rasterized=True)
    for a in ARMS:
        s = mast[mast.arm == a]
        ax.scatter(s.umap1, s.umap2, s=2.0, lw=0, c=ARM_COL[a], alpha=.85,
                   rasterized=True, label=ARM_LAB[a])
    ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + pad)
    bare(ax)
    panel_label(ax, "c", "mast cells")
    ax.legend(loc="upper left", ncol=1, handletextpad=.2, borderpad=.25,
              labelspacing=.22, markerscale=2.0, fontsize=5.2,
              framealpha=.85, facecolor="white", edgecolor="none")

    # ---------------- d. TNFRSF9+ mast cells per arm ----------------------
    spans = [(0, 4), (4, 8), (8, 12)]
    for i, (a, (c0, c1)) in enumerate(zip(ARMS, spans)):
        ax = fig.add_subplot(gs[1, c0:c1])
        sub = mast[mast.arm == a]
        ax.scatter(near.umap1, near.umap2, s=1.0, lw=0, c="#EDEDED", rasterized=True)
        ax.scatter(sub.umap1, sub.umap2, s=1.8, lw=0, c="#C9C9C9", rasterized=True)
        pos = sub[sub.t9pos]
        ax.scatter(pos.umap1, pos.umap2, s=5.5, lw=.25, marker="o",
                   facecolor=ARM_COL[a], edgecolor="#2B2B2B", zorder=5)
        ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + pad)
        bare(ax)
        n, k = len(sub), int(sub.t9pos.sum())
        if i == 0:
            panel_label(ax, "d", f"TNFRSF9$^+$ mast cells · {ARM_LAB[a]}")
        else:
            ax.set_title(ARM_LAB[a], loc="left", fontsize=6.6, pad=4)
        ax.text(.5, .01, f"{k}/{n}  ({100*k/n:.1f}%)", fontsize=6,
                transform=ax.transAxes, ha="center", color="#2B2B2B")

    fig.savefig(FIG / "fig1_umap_atlas.pdf")
    fig.savefig(FIG / "fig1_umap_atlas.png")
    print("wrote fig1_umap_atlas.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
