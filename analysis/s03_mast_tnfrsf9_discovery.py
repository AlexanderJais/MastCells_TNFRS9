"""TNFRSF9 in cutaneous mast cells — discovery cohort GSE204762.

THE QUESTION (CLAUDE.md §1): is TNFRSF9 expressed in cutaneous mast cells in
atopic dermatitis, and does it change between healthy and AD skin?

THE SUBJECT IS MAST CELLS (§2). Fibroblasts, keratinocytes, T/NK and macrophages
are carried only as control populations, to show that a mast-cell result is
specific and not a batch or global effect.

§4 — mast-cell-derived TNFRSF9 is ONE quantity with two factors, always reported
together and as their product:

    mast_TNFRSF9_fraction_of_tissue  =  mast_UMI_share  x  TNFRSF9_rate_in_mast
                                         (abundance)      (per-cell expression)

The identity is exact, which is why the abundance factor is expressed as a UMI
share as well as a cell fraction. The per-cell factor is conditional on which
mast cells survived dissociation (§4) — the Visium arm exists to test that.

§6 controls on every claim: log-depth offset; donor-level aggregation with
donor-clustered SEs; identical test in fibroblasts and keratinocytes; per-group
median depth; marker-QC mast cells; dissociation-stress burden per group.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import HSP_GENES, IEG_GENES, TARGET, TARGET_LIGAND  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")
CONTROLS = ["Fibroblasts", "Keratinocytes", "T/NK", "Macrophages"]

SPECS = {
    G.PRIMARY:          ("arm_pooled", ("Healthy", "AD")),
    "Healthy_vs_AD_NL": ("arm", ("Healthy", "AD_NL")),
    "Healthy_vs_AD_LS": ("arm", ("Healthy", "AD_LS")),
    G.REFINEMENT:       ("arm", ("AD_NL", "AD_LS")),
}


def poisson_fit(d: pd.DataFrame, ycol: str, offcol: str, key: str,
                ref: str, alt: str, cluster: str = "donor") -> dict:
    d = d[d[key].isin([ref, alt])].copy()
    d = d[d[offcol] > 0]
    y = d[ycol].astype(float).values
    x = (d[key] == alt).astype(float).values
    off = np.log(d[offcol].astype(float).values)
    X = sm.add_constant(np.column_stack([x]), has_constant="add")
    G.check_model_terms(["intercept", key, f"log10_depth(offset=log {offcol})"],
                        grouping=cluster)
    out = dict(n_ref=int((~(d[key] == alt)).sum()), n_alt=int((d[key] == alt).sum()),
               donors_ref=int(d.loc[d[key] == ref, cluster].nunique()),
               donors_alt=int(d.loc[d[key] == alt, cluster].nunique()))
    r_ref = d.loc[d[key] == ref, ycol].sum() / d.loc[d[key] == ref, offcol].sum()
    r_alt = d.loc[d[key] == alt, ycol].sum() / d.loc[d[key] == alt, offcol].sum()
    out.update(rate_ref=r_ref, rate_alt=r_alt)
    try:
        res = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(
            cov_type="cluster", cov_kwds={"groups": d[cluster].values})
        b, se = float(res.params[1]), float(res.bse[1])
        out.update(log2FC=b / np.log(2), se_log2=se / np.log(2),
                   ci_lo=(b - 1.96 * se) / np.log(2),
                   ci_hi=(b + 1.96 * se) / np.log(2), p=float(res.pvalues[1]))
    except Exception as e:
        out.update(log2FC=np.nan, p=np.nan, note=str(e)[:60])
    return out


def binom_fit(d: pd.DataFrame, succ: str, tot: str, key: str, ref: str, alt: str,
              cluster: str = "donor") -> dict:
    d = d[d[key].isin([ref, alt])].copy()
    y = np.column_stack([d[succ].values, (d[tot] - d[succ]).values]).astype(float)
    x = (d[key] == alt).astype(float).values
    X = sm.add_constant(np.column_stack([x]), has_constant="add")
    G.check_model_terms(["intercept", key, "log10_depth(binomial denominator=all cells)"],
                        grouping=cluster)
    out = dict(n_ref=int((~(d[key] == alt)).sum()), n_alt=int((d[key] == alt).sum()),
               donors_ref=int(d.loc[d[key] == ref, cluster].nunique()),
               donors_alt=int(d.loc[d[key] == alt, cluster].nunique()),
               pct_ref=100 * d.loc[d[key] == ref, succ].sum() / d.loc[d[key] == ref, tot].sum(),
               pct_alt=100 * d.loc[d[key] == alt, succ].sum() / d.loc[d[key] == alt, tot].sum())
    try:
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d[cluster].values})
        b, se = float(res.params[1]), float(res.bse[1])
        out.update(log2OR=b / np.log(2), ci_lo=(b - 1.96 * se) / np.log(2),
                   ci_hi=(b + 1.96 * se) / np.log(2), p=float(res.pvalues[1]))
    except Exception as e:
        out.update(log2OR=np.nan, p=np.nan, note=str(e)[:60])
    return out


def build_sample_table(obs: pd.DataFrame) -> pd.DataFrame:
    """One row per donor x arm: the §4 factors and their product."""
    obs = obs.copy()
    obs["is_mast"] = obs["mast_strict"].values
    rows = []
    for (donor, arm), g in obs.groupby(["donor", "arm"], observed=True):
        m = g[g.is_mast]
        rec = dict(
            donor=donor, arm=arm,
            arm_pooled="Healthy" if arm == "Healthy" else "AD",
            n_cells=len(g), n_mast=len(m),
            tissue_umi=int(g.depth_retained.sum()),
            mast_umi=int(m.depth_retained.sum()),
            mast_median_depth=float(m.depth_retained.median()) if len(m) else np.nan,
            tissue_median_depth=float(g.depth_retained.median()),
            t9_mast=int(m[f"g_{TARGET}"].sum()) if len(m) else 0,
            t9_tissue=int(g[f"g_{TARGET}"].sum()),
            t9l_mast=int(m[f"g_{TARGET_LIGAND}"].sum()) if len(m) else 0,
            t9l_tissue=int(g[f"g_{TARGET_LIGAND}"].sum()),
            n_mast_t9pos=int((m[f"g_{TARGET}"] > 0).sum()) if len(m) else 0,
            chemistry=g["chemistry_10x"].astype(str).mode().iat[0] if "chemistry_10x" in g else "?",
        )
        # control populations: same quantities, same test (§6)
        for c in CONTROLS:
            cc = g[g["Cell type"].astype(str) == c]
            rec[f"n_{c}"] = len(cc)
            rec[f"umi_{c}"] = int(cc.depth_retained.sum())
            rec[f"t9_{c}"] = int(cc[f"g_{TARGET}"].sum()) if len(cc) else 0
        # §6 dissociation-stress burden, mast cells and whole sample
        hsp = [f"g_{x}" for x in HSP_GENES if f"g_{x}" in g.columns]
        ieg = [f"g_{x}" for x in IEG_GENES if f"g_{x}" in g.columns]
        rec["stress_mast"] = float(m[hsp + ieg].sum().sum() / max(m.depth_retained.sum(), 1) * 1e4) if len(m) else np.nan
        rec["stress_tissue"] = float(g[hsp + ieg].sum().sum() / max(g.depth_retained.sum(), 1) * 1e4)
        rows.append(rec)
    s = pd.DataFrame(rows)

    # ---- §4: the two factors and their exact product ----------------------
    s["f1_mast_umi_share"] = s.mast_umi / s.tissue_umi
    s["f2_t9_rate_in_mast"] = s.t9_mast / s.mast_umi.replace(0, np.nan)
    s["product_mast_t9_share"] = s.t9_mast / s.tissue_umi          # == f1 * f2
    s["mast_pct_cells"] = 100 * s.n_mast / s.n_cells
    s["t9_per_mast_cell"] = s.t9_mast / s.n_mast.replace(0, np.nan)
    s["pct_mast_t9pos"] = 100 * s.n_mast_t9pos / s.n_mast.replace(0, np.nan)
    s["mast_share_of_tissue_t9"] = s.t9_mast / s.t9_tissue.replace(0, np.nan)
    return s


def main() -> int:
    obs = pd.read_parquet(PROC / "obs_all.parquet")
    obs = obs[obs.arm.isin(["Healthy", "AD_NL", "AD_LS"])]
    s = build_sample_table(obs)
    s.to_csv(TAB / "sc_sample_level.csv", index=False)

    order = G.comparison_order(list(SPECS))

    print("=" * 100)
    print("COHORT (GSE204762)")
    print(s.groupby("arm", observed=True).agg(
        donors=("donor", "nunique"), samples=("arm", "size"),
        cells=("n_cells", "sum"), mast=("n_mast", "sum"),
        mast_pct=("mast_pct_cells", "median"),
        mast_depth=("mast_median_depth", "median"),
        tissue_depth=("tissue_median_depth", "median")).to_string())

    # ---------------- is TNFRSF9 expressed in mast cells at all? -----------
    print("\n" + "=" * 100)
    print("Q1  IS TNFRSF9 EXPRESSED IN CUTANEOUS MAST CELLS?")
    tot_mast = s.n_mast.sum()
    tot_pos = s.n_mast_t9pos.sum()
    print(f"  mast cells (marker-QC): {tot_mast}")
    print(f"  TNFRSF9+ mast cells:    {tot_pos} ({100*tot_pos/tot_mast:.2f}%)")
    print(f"  TNFRSF9 counts in mast: {s.t9_mast.sum()}  over {s.mast_umi.sum():,} mast UMI")
    print(f"  rate: {1e4*s.t9_mast.sum()/s.mast_umi.sum():.3f} counts per 10k mast UMI")
    by = s.groupby("arm", observed=True).apply(lambda d: pd.Series({
        "mast_cells": d.n_mast.sum(), "t9+_cells": d.n_mast_t9pos.sum(),
        "pct_t9+": 100 * d.n_mast_t9pos.sum() / max(d.n_mast.sum(), 1),
        "t9_counts": d.t9_mast.sum(),
        "cp10k": 1e4 * d.t9_mast.sum() / max(d.mast_umi.sum(), 1),
        "median_mast_depth": d.mast_median_depth.median()}))
    print(by.to_string())

    # ---------------- §4 factors and product, in contract order ------------
    print("\n" + "=" * 100)
    print("Q2  HEALTHY vs AD  — §4 abundance x per-cell expression = product")
    rows = []
    for comp in order:
        key, (ref, alt) = SPECS[comp]
        a = binom_fit(s, "n_mast", "n_cells", key, ref, alt)
        a.update(comparison=comp, quantity="F1 mast abundance (% of cells)", metric="log2OR")
        b = poisson_fit(s, "t9_mast", "mast_umi", key, ref, alt)
        b.update(comparison=comp, quantity="F2 TNFRSF9 per mast UMI", metric="log2FC")
        c = poisson_fit(s, "t9_mast", "tissue_umi", key, ref, alt)
        c.update(comparison=comp, quantity="PRODUCT mast TNFRSF9 / tissue UMI", metric="log2FC")
        rows += [a, b, c]
    F = pd.DataFrame(rows)
    F.to_csv(TAB / "sc_factors_product.csv", index=False)
    cols = ["comparison", "quantity", "log2OR", "log2FC", "ci_lo", "ci_hi", "p",
            "pct_ref", "pct_alt", "rate_ref", "rate_alt", "donors_ref", "donors_alt"]
    print(F.reindex(columns=cols).to_string(index=False))

    # ---------------- §6 control populations, identical test ---------------
    print("\n" + "=" * 100)
    print("§6 CONTROL POPULATIONS — same test, TNFRSF9 per population UMI")
    rows = []
    for comp in order:
        key, (ref, alt) = SPECS[comp]
        r = poisson_fit(s, "t9_mast", "mast_umi", key, ref, alt)
        r.update(comparison=comp, population="Mast (SUBJECT)")
        rows.append(r)
        for c in CONTROLS:
            r = poisson_fit(s, f"t9_{c}", f"umi_{c}", key, ref, alt)
            r.update(comparison=comp, population=c)
            rows.append(r)
    C = pd.DataFrame(rows)
    C.to_csv(TAB / "sc_control_populations.csv", index=False)
    print(C[["comparison", "population", "log2FC", "ci_lo", "ci_hi", "p",
             "rate_ref", "rate_alt"]].to_string(index=False))

    # ---------------- §6 stress burden -------------------------------------
    print("\n" + "=" * 100)
    print("§6 DISSOCIATION-STRESS BURDEN (HSP+IEG counts per 10k UMI)")
    print(s.groupby("arm", observed=True)[["stress_mast", "stress_tissue"]]
          .median().round(1).to_string())

    # ---------------- mast share of the tissue TNFRSF9 pool ----------------
    print("\n" + "=" * 100)
    print("MAST-CELL SHARE OF ALL TNFRSF9 IN THE TISSUE")
    sh = s.groupby("arm", observed=True).apply(lambda d: pd.Series({
        "t9_mast": d.t9_mast.sum(), "t9_tissue": d.t9_tissue.sum(),
        "mast_share_pct": 100 * d.t9_mast.sum() / max(d.t9_tissue.sum(), 1),
        "mast_umi_share_pct": 100 * d.mast_umi.sum() / d.tissue_umi.sum()}))
    print(sh.to_string())
    sh.to_csv(TAB / "sc_mast_share_of_tissue_t9.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
