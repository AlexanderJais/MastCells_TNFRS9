"""Figure 1 — the mast-cell compartment and TNFRSF9 in the discovery atlas.

Subject: MAST CELLS (CLAUDE.md §2). Other populations are drawn only to place
the mast-cell island in context and to show the depth control.
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
from matplotlib.lines import Line2D  # noqa: E402

import guardrails as G  # noqa: E402
from genes import MAST_QC, TARGET  # noqa: E402
from palette import DOUBLE_COL, GROUP_COLORS, MUTED, set_style  # noqa: E402

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

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.82))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.20, 1.0, 0.9],
                          hspace=0.75, wspace=0.42)

    # ---------------- a. atlas by cell type -------------------------------
    ax = fig.add_subplot(gs[0, :2])
    order = [c for c in CT_COL if c in set(d.ct)]
    rest = d[d.ct != "Mast"]
    ax.scatter(rest.umap1, rest.umap2, s=.6, lw=0, rasterized=True,
               c=[CT_COL.get(c, "#CCCCCC") for c in rest.ct], alpha=.55)
    mm = d[d.ct == "Mast"]
    ax.scatter(mm.umap1, mm.umap2, s=2.6, lw=0, rasterized=True,
               c=MUTED["mast"], alpha=.95)
    for c in order:
        sub = d[d.ct == c]
        if len(sub) < 200:
            continue
        ax.text(sub.umap1.median(), sub.umap2.median(),
                "Mast" if c == "Mast" else c, fontsize=5.4, ha="center", va="center",
                fontweight="bold" if c == "Mast" else "normal",
                color="#2B2B2B" if c != "Mast" else MUTED["mast"])
    ax.set_xlim(*np.percentile(d.umap1, [0.02, 99.98]))
    ax.set_ylim(*np.percentile(d.umap2, [0.02, 99.98]))
    ax.set_title(f"a   {len(d):,} cells, 27 donors", loc="left", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    for s in ax.spines.values():
        s.set_visible(False)

    # ---------------- b. mast marker score --------------------------------
    ax = fig.add_subplot(gs[0, 2])
    o = d.sort_values("mast_score")
    sc = ax.scatter(o.umap1, o.umap2, s=.6, lw=0, rasterized=True,
                    c=np.log1p(o.mast_score), cmap="BuPu", vmin=0, vmax=3)
    ax.set_title("b   tryptase/CPA3", loc="left", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(sc, ax=ax, fraction=.04, pad=.02)
    cb.set_label("log(1+counts)", fontsize=5.5)
    cb.ax.tick_params(labelsize=5)

    # ---------------- c. mast island, by arm ------------------------------
    mast = d[d.mast_strict]
    x0, x1 = mast.umap1.quantile([.005, .995])
    y0, y1 = mast.umap2.quantile([.005, .995])
    pad = 1.2
    ax = fig.add_subplot(gs[0, 3])
    near = d[(d.umap1.between(x0 - pad, x1 + pad)) & (d.umap2.between(y0 - pad, y1 + pad))]
    ax.scatter(near.umap1, near.umap2, s=1.2, lw=0, c="#E2E2E2", rasterized=True)
    for a in ARMS:
        s = mast[mast.arm == a]
        ax.scatter(s.umap1, s.umap2, s=2.2, lw=0, c=ARM_COL[a], alpha=.85,
                   rasterized=True, label=ARM_LAB[a])
    ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + pad)
    ax.set_title("c   mast cells", loc="left", fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    ax.legend(loc="lower center", bbox_to_anchor=(.5, -.3), ncol=1,
              handletextpad=.2, borderpad=.1, labelspacing=.25, markerscale=2.2)
    for s in ax.spines.values():
        s.set_visible(False)

    # ---------------- d. TNFRSF9+ mast cells per arm ----------------------
    for i, a in enumerate(ARMS):
        ax = fig.add_subplot(gs[1, i])
        sub = mast[mast.arm == a]
        ax.scatter(near.umap1, near.umap2, s=1.0, lw=0, c="#EDEDED", rasterized=True)
        ax.scatter(sub.umap1, sub.umap2, s=2.0, lw=0, c="#C9C9C9", rasterized=True)
        pos = sub[sub.t9pos]
        ax.scatter(pos.umap1, pos.umap2, s=13, lw=.35, marker="o",
                   facecolor=ARM_COL[a], edgecolor="#2B2B2B", zorder=5)
        ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + pad)
        n, k = len(sub), int(sub.t9pos.sum())
        ax.set_title(f"{'d' if i == 0 else ''}   {ARM_LAB[a]}",
                     loc="left", fontweight="bold" if i == 0 else "normal")
        ax.text(.5, -.09, f"{k}/{n} TNFRSF9$^+$ ({100*k/n:.1f}%)", fontsize=6,
                transform=ax.transAxes, ha="center", color="#2B2B2B")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)

    # ---------------- e. the depth control (§6) ---------------------------
    ax = fig.add_subplot(gs[1, 3])
    cts = ["Mast", "T/NK", "Macrophages", "Fibroblasts", "Keratinocytes"]
    data = [d.loc[d.ct == c, "depth_retained"].values for c in cts]
    bp = ax.boxplot(data, vert=False, widths=.6, showfliers=False,
                    patch_artist=True, medianprops=dict(color="#2B2B2B", lw=.9))
    for patch, c in zip(bp["boxes"], cts):
        patch.set_facecolor(MUTED["mast"] if c == "Mast" else "#D8D8D8")
        patch.set_edgecolor("#6B6B6B"); patch.set_linewidth(.5)
    ax.set_yticks(range(1, len(cts) + 1))
    ax.set_yticklabels(["Mast" if c == "Mast" else c for c in cts], fontsize=6)
    ax.set_xscale("log")
    ax.set_xlabel("UMI per cell", labelpad=1)
    ax.set_title("e   sequencing depth", loc="left", fontweight="bold")

    # ---------------- f. detection vs depth -------------------------------
    ax = fig.add_subplot(gs[2, :2])
    bins = [0, 200, 400, 800, 1600, 3200, np.inf]
    lab = ["<200", "200–400", "400–800", "800–1600", "1600–3200", ">3200"]
    mast2 = mast.copy()
    mast2["b"] = pd.cut(mast2.depth_retained, bins, labels=lab)
    w = .38
    xs = np.arange(len(lab))
    for j, (a, off) in enumerate([("Healthy", -w / 2), ("AD", w / 2)]):
        sub = mast2[(mast2.arm == "Healthy") if a == "Healthy" else (mast2.arm != "Healthy")]
        vals, ns = [], []
        for b in lab:
            s = sub[sub.b == b]
            vals.append(1e4 * s[f"g_{TARGET}"].sum() / max(s.depth_retained.sum(), 1))
            ns.append(len(s))
        ax.bar(xs + off, vals, w, color=MUTED["healthy"] if a == "Healthy" else MUTED["ad_all"],
               edgecolor="none", label="Healthy" if a == "Healthy" else "AD")
        for x, v, n in zip(xs + off, vals, ns):
            ax.text(x, v + .02, f"{n}", fontsize=4.6, ha="center", color="#6B6B6B")
    ax.set_xticks(xs); ax.set_xticklabels(lab, fontsize=5.6)
    ax.set_xlabel("mast-cell sequencing depth (UMI)")
    ax.set_ylabel("TNFRSF9 per 10k UMI")
    ax.legend(ncol=2, loc="upper left")
    ax.set_title("f   depth-matched", loc="left", fontweight="bold")

    # ---------------- g. per-donor rates ----------------------------------
    ax = fig.add_subplot(gs[2, 2:])
    don = pd.read_csv(ROOT / "results" / "tables" / "sc_donor_level.csv")
    don = don.sort_values(["arm", "cp10k"], ascending=[True, True])
    cols = [MUTED["healthy"] if a == "Healthy" else MUTED["ad_all"] for a in don.arm]
    ax.bar(np.arange(len(don)), don.cp10k, color=cols, edgecolor="none", width=.72)
    ax.set_xticks(np.arange(len(don)))
    ax.set_xticklabels(don.donor, rotation=90, fontsize=4.8)
    ax.set_ylabel("TNFRSF9 per 10k mast UMI")
    ax.set_title("g   per donor", loc="left", fontweight="bold")
    ax.legend(handles=[Line2D([0], [0], marker="s", lw=0, markersize=4,
                              markerfacecolor=MUTED["healthy"], markeredgewidth=0, label="Healthy"),
                       Line2D([0], [0], marker="s", lw=0, markersize=4,
                              markerfacecolor=MUTED["ad_all"], markeredgewidth=0, label="AD")],
              loc="upper left", ncol=2)

    fig.savefig(FIG / "fig1_umap_atlas.pdf")
    fig.savefig(FIG / "fig1_umap_atlas.png")
    print("wrote fig1_umap_atlas.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
