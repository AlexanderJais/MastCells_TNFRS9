"""Figure 5 — replication (CLAUDE.md §7): what reproduces and what does not.

Subject: MAST CELLS (§2). UMAPs are the authors'-free embeddings computed here
for each replication cohort; mast cells are called by the rule calibrated on the
discovery cohort (s41).
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
from genes import MAST_QC, TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REP = ROOT / "data" / "processed" / "replication"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])
THR = 50.0

COH = [("REP1 · GSE222840+GSE173205", "rep1_obs.parquet", "rep1_umap.parquet"),
       ("REP2 · GSE153760 biopsies", "rep2_obs.parquet", "rep2_umap.parquet")]


def main() -> int:
    set_style()
    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.52))
    gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.45,
                          width_ratios=[1, 1, 1, 1.35])

    for row, (title, fo, fu) in enumerate(COH):
        o = pd.read_parquet(REP / fo)
        u = pd.read_parquet(REP / fu)
        # obs already carries umap1/umap2 from build(); drop before joining
        o = o.drop(columns=[c for c in ("umap1", "umap2") if c in o.columns])
        d = o.join(u[["umap1", "umap2"]], how="inner")
        d["tryp"] = 1e4 * d[[f"g_{g}" for g in MAST_QC]].sum(axis=1) / d.depth_retained
        d["mast"] = d.mast_qc_pass & (d.tryp >= THR)
        d["t9pos"] = d[f"g_{TARGET}"] > 0

        # a: cluster map with mast highlighted
        ax = fig.add_subplot(gs[row, 0])
        ax.scatter(d.umap1, d.umap2, s=.5, lw=0, c="#DFDFDF", rasterized=True)
        m = d[d.mast]
        ax.scatter(m.umap1, m.umap2, s=3.5, lw=0, c=MUTED["mast"], rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title(f"{'a' if row == 0 else 'd'}   {title}", loc="left",
                     fontweight="bold", fontsize=6.2)
        ax.text(.02, .03, f"{len(m):,} mast / {len(d):,} cells "
                          f"({100*len(m)/len(d):.2f}%)",
                transform=ax.transAxes, fontsize=5.2, color="#4A4A4A")

        # b: tryptase overlay
        ax = fig.add_subplot(gs[row, 1])
        s_ = d.sort_values("tryp")
        sc = ax.scatter(s_.umap1, s_.umap2, s=.5, lw=0, rasterized=True,
                        c=np.log1p(s_.tryp), cmap="BuPu", vmin=0, vmax=6)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_title("tryptase/CPA3", loc="left", fontsize=5.8)
        cb = fig.colorbar(sc, ax=ax, fraction=.04, pad=.02)
        cb.ax.tick_params(labelsize=4.4)

        # c: mast cells by arm, TNFRSF9+ marked
        ax = fig.add_subplot(gs[row, 2])
        ax.scatter(m.umap1, m.umap2, s=3.5, lw=0, c="#CFCFCF", rasterized=True)
        for arm, col in (("Healthy", MUTED["healthy"]), ("AD", MUTED["ad_all"])):
            p = m[(m.arm == arm) & m.t9pos]
            ax.scatter(p.umap1, p.umap2, s=15, lw=.35, facecolor=col,
                       edgecolor="#1F1F1F", zorder=5, label=arm)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        nh = int(m[(m.arm == "Healthy") & m.t9pos].shape[0])
        na = int(m[(m.arm == "AD") & m.t9pos].shape[0])
        ax.set_title(f"TNFRSF9$^+$: {nh} healthy, {na} AD", loc="left", fontsize=6)
        if row == 0:
            ax.legend(fontsize=5, loc="lower right", markerscale=.8,
                      handletextpad=.15, borderpad=.15)

    # ---- forest plot across cohorts --------------------------------------
    ax = fig.add_subplot(gs[:, 3])
    M = pd.read_csv(TAB / "meta_analysis.csv")
    labels, xs, los, his, cols = [], [], [], [], []
    for r in M.itertuples():
        short = (r.cohort.replace(" (GSE204762, 3')", "").replace(" (GSE222840+GSE173205, 5')", "")
                 .replace(" (GSE153760 biopsies, 3' v3)", ""))
        labels.append(short)
        xs.append(r.log2FC); los.append(r.ci_lo); his.append(r.ci_hi)
        cols.append(MUTED["accent"] if "POOLED" in r.cohort else MUTED["mast"])
    y = np.arange(len(labels))
    for i in range(len(labels)):
        ax.plot([los[i], his[i]], [i, i], color="#6B6B6B", lw=1)
        ax.plot(xs[i], i, "D" if "POOLED" in labels[i] else "o",
                ms=4.5 if "POOLED" in labels[i] else 4, color=cols[i], zorder=5)
    ax.axvline(0, color="#B0B0B0", lw=.6, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ TNFRSF9 per mast UMI, AD vs healthy", fontsize=5.8)
    ax.set_title("g   direction consistent, magnitude not", loc="left",
                 fontweight="bold", fontsize=6.4)
    ax.text(.02, .02, "I² = 58%", transform=ax.transAxes, fontsize=5.6,
            color="#4A4A4A")

    fig.savefig(FIG / "fig5_replication.pdf")
    fig.savefig(FIG / "fig5_replication.png")
    print("wrote fig5_replication.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
