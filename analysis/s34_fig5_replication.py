"""Figure 5 — replication in an independent cohort (CLAUDE.md §7).

Panels mirror Figures 1 and 2 so the replication is read against the discovery
cohort on identical axes: the mast-cell map, the TNFRSF9+ mast cells on it, and
the TNFRSF9+ fraction with statistics.

Panel d records mast-cell recovery per cohort. It is the reason GSE222840 +
GSE173205 does not enter the analysis: 0.11% of cells recovered as mast cells
against 1.31% and 4.41% elsewhere, giving 27 mast cells in the healthy arm. A
cohort that cannot recover the cell type cannot test a gene in it.

Subject: MAST CELLS (§2).
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
from genes import MAST_QC, TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
REP = ROOT / "data" / "processed" / "replication"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])
THR = 50.0

ARM_COL = {"Healthy": MUTED["healthy"], "AD": MUTED["ad_all"]}


def load_rep2() -> pd.DataFrame:
    o = pd.read_parquet(REP / "rep2_obs.parquet")
    u = pd.read_parquet(REP / "rep2_umap.parquet")
    o = o.drop(columns=[c for c in ("umap1", "umap2") if c in o.columns])
    d = o.join(u[["umap1", "umap2"]], how="inner")
    d["tryp"] = 1e4 * d[[f"g_{g}" for g in MAST_QC]].sum(axis=1) / d.depth_retained
    d["mast"] = d.mast_qc_pass & (d.tryp >= THR)
    d["t9pos"] = d[f"g_{TARGET}"] > 0
    return d


def bare(ax):
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_aspect("equal")
    for s in ax.spines.values():
        s.set_visible(False)


def main() -> int:
    set_style()
    d = load_rep2()
    mast = d[d.mast]

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.34))
    gs = fig.add_gridspec(1, 4, wspace=0.62, width_ratios=[1, 1, 1.2, 1.0])

    # ---- a. where the mast cells are -------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(d.umap1, d.umap2, s=.6, lw=0, c="#DFDFDF", rasterized=True)
    ax.scatter(mast.umap1, mast.umap2, s=2.4, lw=0, c=MUTED["mast"], rasterized=True)
    bare(ax)
    ax.set_title(r"$\bf{a}$   GSE153760 · " f"{len(mast):,} mast cells",
                 loc="left", fontsize=6.6, pad=4)

    # ---- b. TNFRSF9+ mast cells on the mast island -----------------------
    ax = fig.add_subplot(gs[0, 1])
    x0, x1 = mast.umap1.quantile([.005, .995]); y0, y1 = mast.umap2.quantile([.005, .995])
    pad = 1.0
    ax.scatter(mast.umap1, mast.umap2, s=2.4, lw=0, c="#CFCFCF", rasterized=True)
    for arm in ("Healthy", "AD"):
        p = mast[(mast.arm == arm) & mast.t9pos]
        ax.scatter(p.umap1, p.umap2, s=7, lw=.3, facecolor=ARM_COL[arm],
                   edgecolor="#1F1F1F", zorder=5, label=arm)
    ax.set_xlim(x0 - pad, x1 + pad); ax.set_ylim(y0 - pad, y1 + pad)
    bare(ax)
    ax.set_title(r"$\bf{b}$   TNFRSF9$^+$ mast cells", loc="left", fontsize=6.6, pad=4)
    ax.legend(fontsize=5.4, loc="lower center", bbox_to_anchor=(.5, -.14), ncol=2,
              frameon=False, markerscale=1.2, handletextpad=.15, columnspacing=.8)

    # ---- c. the TNFRSF9+ fraction, both cohorts, per sample --------------
    ax = fig.add_subplot(gs[0, 2])
    disc = pd.read_csv(TAB / "sc_sample_level.csv")
    disc["pct"] = 100 * disc.n_mast_t9pos / disc.n_mast.replace(0, np.nan)
    disc["armp"] = np.where(disc.arm == "Healthy", "Healthy", "AD")
    r2 = pd.read_csv(TAB / "replication_donor_level.csv")
    r2 = r2[r2.cohort.str.startswith("REP2")].copy()
    r2["pct"] = 100 * r2.n_mast_t9pos / r2.n_mast.replace(0, np.nan)
    r2["armp"] = r2.arm

    rng = np.random.default_rng(5)
    pos, labels, notes = [], [], []
    for k, (name, df) in enumerate([("Discovery", disc), ("GSE153760", r2)]):
        for j, arm in enumerate(("Healthy", "AD")):
            i = k * 2.6 + j
            v = df.loc[df.armp == arm, "pct"].astype(float).dropna().values
            ax.scatter(i + rng.uniform(-.16, .16, len(v)), v, s=14, lw=.3,
                       facecolor=ARM_COL[arm], edgecolor="#333333", zorder=3)
            ax.plot([i - .28, i + .28], [np.median(v)] * 2, color="#1F1F1F", lw=1.2, zorder=4)
            pos.append(i); labels.append(arm)
        # pooled Fisher within cohort
        h = df[df.armp == "Healthy"]; a = df[df.armp == "AD"]
        orr, pf = stats.fisher_exact(
            [[int(a.n_mast_t9pos.sum()), int(a.n_mast.sum() - a.n_mast_t9pos.sum())],
             [int(h.n_mast_t9pos.sum()), int(h.n_mast.sum() - h.n_mast_t9pos.sum())]])
        stat = (f"OR {orr:.1f}, P = {pf:.3f}" if pf >= .001
                else f"OR {orr:.1f}, P < 0.001")
        notes.append((k, name, stat))
    ax.set_ylim(-0.4, 12.6)
    for k, name, stat in notes:
        xa = (k * 2.6 + .5) / 3.6 + .06        # data -> axes fraction
        ax.text(xa, 1.005, name, transform=ax.transAxes, ha="center",
                va="bottom", fontsize=5.6, color="#2B2B2B")
        ax.text(xa, .95, stat, transform=ax.transAxes, ha="center",
                va="top", fontsize=5.2, color="#4A4A4A")
    ax.set_xticks(pos); ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("TNFRSF9$^+$ (% of mast cells)", fontsize=6.2)
    ax.set_title(r"$\bf{c}$", loc="left", fontsize=6.6, pad=12)

    # ---- d. mast-cell recovery, and why REP1 is excluded ------------------
    ax = fig.add_subplot(gs[0, 3])
    rec = pd.DataFrame([
        ("GSE153760", 4.41, True), ("GSE204762", 1.31, True),
        ("GSE222840", 0.106, False)],
        columns=["cohort", "pct", "ok"])
    cols = [MUTED["mast"] if o else "#C9C9C9" for o in rec.ok]
    ax.barh(np.arange(len(rec)), rec.pct, color=cols, edgecolor="none", height=.62)
    ax.axvline(0.5, color="#8A8A8A", lw=.7, ls="--")
    ax.text(0.55, 2.62, "threshold", fontsize=5.0, color="#5A5A5A",
            va="top", ha="left")
    for i, r in enumerate(rec.itertuples()):
        ax.text(r.pct * 1.12, i, f"{r.pct:.2f}%", va="center", fontsize=5.4)
    ax.set_yticks(np.arange(len(rec)))
    ax.set_yticklabels(rec.cohort, fontsize=5.8)
    ax.set_ylim(-0.6, 2.9)
    ax.set_xscale("log")
    ax.set_xlabel("mast cells (% of all cells)", fontsize=6.2)
    ax.set_title(r"$\bf{d}$   mast-cell recovery", loc="left", fontsize=6.6, pad=4)

    fig.savefig(FIG / "fig5_replication.pdf")
    fig.savefig(FIG / "fig5_replication.png")
    print("wrote fig5_replication.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
