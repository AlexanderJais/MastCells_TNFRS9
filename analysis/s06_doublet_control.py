"""Doublet control: are TNFRSF9+ "mast cells" really mast-T doublets?

T/NK cells carry TNFRSF9 at a rate comparable to AD mast cells, so a mast-cell/
T-cell doublet would present exactly as a tryptase-positive, TNFRSF9-positive
cell. This is the single most dangerous artefact for the result and it is tested
three ways:

  1. Do TNFRSF9+ mast cells carry more T-cell transcripts than TNFRSF9- mast
     cells *at matched sequencing depth*? (Depth alone raises detection of every
     gene, so the raw comparison is confounded.)
  2. Does the primary Healthy-vs-AD result survive deleting every mast cell that
     carries any T-cell transcript at all — the strictest possible purge?
  3. Does it survive deleting cells with any monocyte/DC transcript too?

Subject: MAST CELLS (CLAUDE.md §2). T-cell markers appear here purely as a
contamination diagnostic.
"""
from __future__ import annotations

import sys
import warnings
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])

TCELL = ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "IL7R"]
MYELOID = ["LYZ", "CD68", "CD163", "AIF1", "C1QA", "ITGAX", "CD14"]


def exact_perm(counts, expo, alt):
    nd, k = len(counts), int(alt.sum())
    masks = np.zeros((comb(nd, k), nd), bool)
    for i, c in enumerate(combinations(range(nd), k)):
        masks[i, list(c)] = True
    ea = masks @ expo
    er = expo.sum() - ea
    ca = masks @ counts
    cr = counts.sum() - ca
    S = np.log((ca + .5) / ea) - np.log((cr + .5) / er)
    i0 = int(np.where((masks == alt[None, :]).all(axis=1))[0][0])
    return S[i0] / np.log(2), float((np.abs(S) >= abs(S[i0]) - 1e-12).mean())


def primary(m: pd.DataFrame) -> dict:
    """Donor-level exact permutation, Healthy vs AD, on the surviving cells."""
    g = m.groupby("donor").agg(
        arm=("arm", lambda v: "Healthy" if (v == "Healthy").all() else "AD"),
        umi=("depth_retained", "sum"), t9=(f"g_{TARGET}", "sum"),
        n=("depth_retained", "size")).reset_index()
    g["arm"] = np.where(g.donor.str.startswith("Healthy"), "Healthy", "AD")
    alt = (g.arm == "AD").values
    lf, p = exact_perm(g.t9.values.astype(float), g.umi.values.astype(float), alt)
    return dict(log2FC=lf, p_perm=p, n_mast=int(g.n.sum()),
                t9_healthy=int(g.loc[~alt, "t9"].sum()), t9_ad=int(g.loc[alt, "t9"].sum()),
                cells_healthy=int(g.loc[~alt, "n"].sum()), cells_ad=int(g.loc[alt, "n"].sum()),
                rate_healthy=1e4 * g.loc[~alt, "t9"].sum() / g.loc[~alt, "umi"].sum(),
                rate_ad=1e4 * g.loc[alt, "t9"].sum() / g.loc[alt, "umi"].sum())


def main() -> int:
    o = pd.read_parquet(PROC / "obs_all.parquet")
    m = o[o.mast_strict & o.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    m["t9pos"] = m[f"g_{TARGET}"] > 0
    m["t_counts"] = m[[f"g_{g}" for g in TCELL if f"g_{g}" in m.columns]].sum(axis=1)
    m["mye_counts"] = m[[f"g_{g}" for g in MYELOID if f"g_{g}" in m.columns]].sum(axis=1)

    # ---- 1. depth-matched contamination check ----------------------------
    print("=" * 92)
    print("1. T-CELL TRANSCRIPT CARRIAGE IN TNFRSF9+ vs TNFRSF9- MAST CELLS")
    bins = [0, 400, 600, 900, 1400, 2200, np.inf]
    m["dbin"] = pd.cut(m.depth_retained, bins)
    rows = []
    for b, g in m.groupby("dbin", observed=True):
        a, c = g[g.t9pos], g[~g.t9pos]
        if len(a) == 0:
            continue
        rows.append(dict(depth_bin=str(b), n_t9pos=len(a), n_t9neg=len(c),
                         pct_Tpos_in_t9pos=100 * (a.t_counts > 0).mean(),
                         pct_Tpos_in_t9neg=100 * (c.t_counts > 0).mean()))
    D = pd.DataFrame(rows)
    print(D.round(2).to_string(index=False))
    # Cochran-Mantel-Haenszel style pooled OR across depth strata
    num = den = 0.0
    for b, g in m.groupby("dbin", observed=True):
        a, c = g[g.t9pos], g[~g.t9pos]
        if len(a) == 0 or len(c) == 0:
            continue
        n = len(g)
        a1, a0 = (a.t_counts > 0).sum(), (a.t_counts == 0).sum()
        c1, c0 = (c.t_counts > 0).sum(), (c.t_counts == 0).sum()
        num += a1 * c0 / n
        den += a0 * c1 / n
    print(f"  depth-stratified OR (T-cell transcript | TNFRSF9+): {num/den if den else np.inf:.2f}")
    D.to_csv(TAB / "sc_doublet_depth_strata.csv", index=False)

    # ---- 2/3. purge and re-test ------------------------------------------
    print("\n" + "=" * 92)
    print("2. DOES THE PRIMARY RESULT SURVIVE PURGING SUSPECT CELLS?")
    variants = {
        "all marker-QC mast cells": m,
        "excluding any T-cell transcript": m[m.t_counts == 0],
        "excluding any T-cell or myeloid transcript": m[(m.t_counts == 0) & (m.mye_counts == 0)],
        "excluding scrublet-flagged": m[m.get("scrublet", pd.Series(False, index=m.index)).astype(str) != "True"],
    }
    out = []
    for name, sub in variants.items():
        r = primary(sub)
        r["cell_set"] = name
        out.append(r)
    R = pd.DataFrame(out)[["cell_set", "n_mast", "cells_healthy", "cells_ad",
                           "t9_healthy", "t9_ad", "rate_healthy", "rate_ad",
                           "log2FC", "p_perm"]]
    R.to_csv(TAB / "sc_doublet_purge.csv", index=False)
    print(R.round(4).to_string(index=False))
    print("\n  If the effect persists after deleting every mast cell carrying a single")
    print("  T-cell transcript, mast-T doublets cannot be its source.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
