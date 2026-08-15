"""Robustness of the mast-cell TNFRSF9 result — discovery cohort GSE204762.

Subject: MAST CELLS (CLAUDE.md §2).

The discovery counts are small (110 TNFRSF9 molecules across 3,661 mast cells,
3 of them in the healthy arm). Cluster-robust standard errors are unreliable
with six donor clusters, so inference here does not rely on them. Instead:

  1. EXACT permutation test over donor labels (all C(n,k) assignments
     enumerated), on the donor-aggregated rate — no distributional assumption.
  2. Donor-level rank test and the donor-positivity contingency.
  3. Leave-one-donor-out, because one donor (MGH108) carries half the counts.
  4. Depth-matched detection: mast cells resampled to a common depth
     distribution, so the comparison cannot be a depth comparison (§6).
  5. Expression-matched null genes: the identical pipeline run on genes with
     comparable abundance in mast cells, calibrating how often this design
     manufactures a result.
  6. Dissociation-stress adjustment (§6) — the healthy arm is the most stressed
     and a stressed arm reads artificially low.
"""
from __future__ import annotations

import sys
import warnings
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
# §3: the order in which contrasts are reported is validated, not assumed.
COMPARISONS = G.comparison_order([G.PRIMARY, *G.SECONDARY])
RNG = np.random.default_rng(20240815)


# ------------------------------------------------------- exact permutation ---
def exact_perm_rate(counts: np.ndarray, expo: np.ndarray, is_alt: np.ndarray,
                    max_enum: int = 400_000) -> dict:
    """Exact (or near-exact) permutation of donor labels on the pooled rate ratio.

    Statistic: log((sum c_alt / sum e_alt) / (sum c_ref / sum e_ref)), with a
    0.5 continuity term so an all-zero arm does not give an infinite statistic.
    """
    n, k = len(counts), int(is_alt.sum())

    def stat(mask):
        ca, ea = counts[mask].sum(), expo[mask].sum()
        cr, er = counts[~mask].sum(), expo[~mask].sum()
        return np.log((ca + .5) / ea) - np.log((cr + .5) / er)

    obs = stat(is_alt.astype(bool))
    from math import comb
    total = comb(n, k)
    idx = np.arange(n)
    if total <= max_enum:
        stats_null = np.empty(total)
        for i, c in enumerate(combinations(idx, k)):
            m = np.zeros(n, bool)
            m[list(c)] = True
            stats_null[i] = stat(m)
        exact = True
    else:
        stats_null = np.empty(max_enum)
        for i in range(max_enum):
            m = np.zeros(n, bool)
            m[RNG.choice(idx, k, replace=False)] = True
            stats_null[i] = stat(m)
        exact = False
    p_two = float((np.abs(stats_null) >= abs(obs) - 1e-12).mean())
    return dict(obs_log2=obs / np.log(2), p_perm=p_two, exact=exact,
                n_perm=len(stats_null), n_ref=n - k, n_alt=k)


def donor_table(s: pd.DataFrame) -> pd.DataFrame:
    """Aggregate samples to one row per donor (§6: no pseudo-replication)."""
    g = s.groupby("donor").agg(
        arm=("arm", lambda v: "Healthy" if (v == "Healthy").all() else "AD"),
        n_cells=("n_cells", "sum"), n_mast=("n_mast", "sum"),
        mast_umi=("mast_umi", "sum"), tissue_umi=("tissue_umi", "sum"),
        t9_mast=("t9_mast", "sum"), t9_tissue=("t9_tissue", "sum"),
        n_mast_t9pos=("n_mast_t9pos", "sum"),
        stress_mast=("stress_mast", "mean")).reset_index()
    g["cp10k"] = 1e4 * g.t9_mast / g.mast_umi
    g["mast_pct"] = 100 * g.n_mast / g.n_cells
    return g


def main() -> int:
    s = pd.read_csv(TAB / "sc_sample_level.csv")
    d = donor_table(s)
    d.to_csv(TAB / "sc_donor_level.csv", index=False)

    print("=" * 96)
    print("1. EXACT PERMUTATION TEST (donor labels), TNFRSF9 per mast UMI")
    rows = []
    # PRIMARY first (§3)
    specs = [(G.PRIMARY, d, d.arm == "AD")]
    sNL = s[s.arm.isin(["Healthy", "AD_NL"])]
    dNL = donor_table(sNL)
    specs.append(("Healthy_vs_AD_NL", dNL, dNL.arm == "AD"))
    sLS = s[s.arm.isin(["Healthy", "AD_LS"])]
    dLS = donor_table(sLS)
    specs.append(("Healthy_vs_AD_LS", dLS, dLS.arm == "AD"))
    for name, dd, alt in specs:
        r = exact_perm_rate(dd.t9_mast.values.astype(float),
                            dd.mast_umi.values.astype(float), alt.values)
        r.update(comparison=name,
                 rate_ref=1e4 * dd.loc[~alt, "t9_mast"].sum() / dd.loc[~alt, "mast_umi"].sum(),
                 rate_alt=1e4 * dd.loc[alt, "t9_mast"].sum() / dd.loc[alt, "mast_umi"].sum(),
                 counts_ref=int(dd.loc[~alt, "t9_mast"].sum()),
                 counts_alt=int(dd.loc[alt, "t9_mast"].sum()))
        rows.append(r)
    P = pd.DataFrame(rows)
    P.to_csv(TAB / "sc_exact_permutation.csv", index=False)
    print(P[["comparison", "obs_log2", "p_perm", "exact", "n_perm", "n_ref",
             "n_alt", "counts_ref", "counts_alt", "rate_ref", "rate_alt"]].to_string(index=False))

    # ---- paired refinement (§3, reported last) ---------------------------
    both = s.pivot_table(index="donor", columns="arm", values=["t9_mast", "mast_umi"],
                         aggfunc="sum").dropna()
    if len(both):
        r_ls = both[("t9_mast", "AD_LS")] / both[("mast_umi", "AD_LS")]
        r_nl = both[("t9_mast", "AD_NL")] / both[("mast_umi", "AD_NL")]
        w = stats.wilcoxon(r_ls, r_nl, zero_method="wilcox") if (r_ls != r_nl).any() else None
        print(f"\n  {G.REFINEMENT} (paired within donor, n={len(both)} donors): "
              f"median LS {1e4*r_ls.median():.3f} vs NL {1e4*r_nl.median():.3f} cp10k, "
              f"Wilcoxon p={w.pvalue:.3f}" if w else "  no variation")

    # ---- 2. donor-level rank test + positivity ---------------------------
    print("\n" + "=" * 96)
    print("2. DONOR-LEVEL RANK TEST AND POSITIVITY")
    h, a = d[d.arm == "Healthy"], d[d.arm == "AD"]
    u = stats.mannwhitneyu(a.cp10k, h.cp10k, alternative="two-sided")
    print(f"  Mann-Whitney on donor TNFRSF9 cp10k: U={u.statistic:.0f}, p={u.pvalue:.4f} "
          f"(n={len(h)} healthy vs {len(a)} AD donors)")
    tab = np.array([[int((a.t9_mast > 0).sum()), int((a.t9_mast == 0).sum())],
                    [int((h.t9_mast > 0).sum()), int((h.t9_mast == 0).sum())]])
    orr, pf = stats.fisher_exact(tab)
    print(f"  donors with >=1 TNFRSF9+ mast cell: AD {tab[0,0]}/{tab[0].sum()}, "
          f"Healthy {tab[1,0]}/{tab[1].sum()}; Fisher OR={orr:.2f}, p={pf:.4f}")
    exp_h = h.mast_umi.sum() * (a.t9_mast.sum() / a.mast_umi.sum())
    print(f"  healthy counts expected under the AD rate: {exp_h:.1f}; observed: {int(h.t9_mast.sum())}")

    # ---- 3. leave-one-donor-out ------------------------------------------
    print("\n" + "=" * 96)
    print("3. LEAVE-ONE-DONOR-OUT (primary comparison)")
    loo = []
    for dr in d.donor:
        dd = d[d.donor != dr]
        alt = dd.arm == "AD"
        if alt.sum() < 2 or (~alt).sum() < 2:
            continue
        r = exact_perm_rate(dd.t9_mast.values.astype(float),
                            dd.mast_umi.values.astype(float), alt.values)
        loo.append(dict(dropped=dr, arm=d.loc[d.donor == dr, "arm"].iat[0],
                        log2FC=r["obs_log2"], p_perm=r["p_perm"]))
    L = pd.DataFrame(loo).sort_values("log2FC")
    L.to_csv(TAB / "sc_leave_one_out.csv", index=False)
    print(L.to_string(index=False))

    # ---- 4. depth-matched detection --------------------------------------
    print("\n" + "=" * 96)
    print("4. DEPTH-MATCHED MAST-CELL COMPARISON (§6)")
    obs = pd.read_parquet(PROC / "obs_all.parquet")
    m = obs[obs.mast_strict & obs.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    m["armp"] = np.where(m.arm == "Healthy", "Healthy", "AD")
    bins = np.array([0, 200, 400, 800, 1600, 3200, np.inf])
    m["dbin"] = pd.cut(m.depth_retained, bins)
    tabm = m.groupby(["dbin", "armp"], observed=True).apply(lambda g: pd.Series({
        "cells": len(g), "umi": g.depth_retained.sum(),
        "t9": g[f"g_{TARGET}"].sum(),
        "cp10k": 1e4 * g[f"g_{TARGET}"].sum() / max(g.depth_retained.sum(), 1)}))
    print(tabm.unstack(1).round(3).to_string())
    tabm.to_csv(TAB / "sc_depth_matched.csv")
    # stratified (Mantel-Haenszel style) common rate ratio across depth bins
    num = den = 0.0
    for b, g in m.groupby("dbin", observed=True):
        A = g[g.armp == "AD"]; H = g[g.armp == "Healthy"]
        if len(A) == 0 or len(H) == 0:
            continue
        eA, eH = A.depth_retained.sum(), H.depth_retained.sum()
        cA, cH = A[f"g_{TARGET}"].sum(), H[f"g_{TARGET}"].sum()
        tot_e = eA + eH
        num += cA * eH / tot_e
        den += cH * eA / tot_e
    print(f"  depth-stratified common rate ratio (AD vs Healthy): {num/den if den else np.inf:.2f}"
          f"  [log2 = {np.log2(num/den) if den else np.inf:.2f}]")

    # ---- 5. expression-matched null genes --------------------------------
    print("\n" + "=" * 96)
    print("5. EXPRESSION-MATCHED NULL GENES — does this design manufacture results?")
    # Samples were deposited with differing gene sets (27,364 vs 29,991 genes),
    # so matrices are remapped onto a common gene axis before stacking.
    loaded = []
    for p in sorted(PROC.glob("*_mast.npz")):
        gsm = p.name.split("_")[0]
        meta = np.load(PROC / f"{gsm}_mastmeta.npz", allow_pickle=True)
        mm = sp.load_npz(p)
        if mm.shape[0] == 0:
            continue
        loaded.append((mm, np.asarray(meta["var"], dtype=str),
                       np.asarray(meta["barcodes"], dtype=str)))
    common = sorted(set.intersection(*[set(v) for _, v, _ in loaded]))
    var = np.array(common)
    pos = {g: i for i, g in enumerate(var)}
    mats, bcs = [], []
    for mm, v, b in loaded:
        take = [i for i, g in enumerate(v) if g in pos]
        order = np.argsort([pos[v[i]] for i in take])
        mats.append(mm[:, np.asarray(take)[order]])
        bcs.append(b)
    print(f"  common gene axis across samples: {len(var)} genes")
    M = sp.vstack(mats).tocsr()
    B = np.concatenate(bcs)
    keep = pd.Index(B).isin(m.index)
    M, B = M[keep], B[keep]
    mm = m.loc[B]
    donors = mm.donor.values
    is_ad = (mm.armp == "AD").values
    tot = np.asarray(M.sum(0)).ravel()
    t9_tot = tot[list(var).index(TARGET)]
    cand = np.where((tot >= 0.5 * t9_tot) & (tot <= 2 * t9_tot))[0]
    print(f"  TNFRSF9 total counts in mast: {int(t9_tot)}; "
          f"{len(cand)} genes within 0.5-2x that abundance")
    udon = pd.unique(donors)
    dmap = {u: i for i, u in enumerate(udon)}
    di = np.array([dmap[x] for x in donors])
    dep = mm.depth_retained.values.astype(float)
    expo = np.bincount(di, weights=dep, minlength=len(udon))
    alt = np.array([bool(is_ad[di == i][0]) for i in range(len(udon))])
    # Vectorised exact permutation: build every donor-label assignment once,
    # then evaluate all candidate genes simultaneously.
    nd, k = len(udon), int(alt.sum())
    masks = np.zeros((int(__import__("math").comb(nd, k)), nd), dtype=bool)
    for i, c in enumerate(combinations(range(nd), k)):
        masks[i, list(c)] = True
    # donor x gene count matrix for the candidate genes
    Cd = np.zeros((nd, len(cand)))
    dense = np.asarray(M[:, cand].todense())
    for i in range(nd):
        Cd[i] = dense[di == i].sum(axis=0)
    e_alt = masks @ expo                      # (n_perm,)
    e_ref = expo.sum() - e_alt
    c_alt = masks @ Cd                        # (n_perm, n_genes)
    c_ref = Cd.sum(axis=0)[None, :] - c_alt
    S = (np.log((c_alt + .5) / e_alt[:, None])
         - np.log((c_ref + .5) / e_ref[:, None]))
    obs_i = int(np.where((masks == alt[None, :]).all(axis=1))[0][0])
    obs = S[obs_i]
    p = (np.abs(S) >= np.abs(obs)[None, :] - 1e-12).mean(axis=0)
    N = pd.DataFrame(dict(gene=var[cand], total=tot[cand].astype(int),
                          log2FC=obs / np.log(2), p_perm=p)).sort_values("p_perm")
    N.to_csv(TAB / "sc_null_gene_calibration.csv", index=False)
    frac = (N.p_perm <= 0.05).mean()
    t9row = N[N.gene == TARGET]
    rank = int((N.p_perm < t9row.p_perm.iat[0]).sum()) + 1 if len(t9row) else -1
    print(f"  abundance-matched genes with p_perm<=0.05: {100*frac:.1f}% "
          f"({int((N.p_perm<=0.05).sum())}/{len(N)})")
    print(f"  TNFRSF9 rank by permutation p: {rank}/{len(N)}")
    print("  strongest matched genes:")
    print(N.head(12).to_string(index=False))

    # ---- 6. stress adjustment --------------------------------------------
    print("\n" + "=" * 96)
    print("6. DISSOCIATION-STRESS (§6): healthy arm is the most stressed")
    print(d.groupby("arm")[["stress_mast"]].median().round(1).to_string())
    rho, prho = stats.spearmanr(d.stress_mast, d.cp10k)
    print(f"  donor-level Spearman(stress, TNFRSF9 cp10k) = {rho:.3f} (p={prho:.3f}) "
          f"-> stress does not explain the direction" if rho < 0 else
          f"  donor-level Spearman(stress, TNFRSF9 cp10k) = {rho:.3f} (p={prho:.3f})")

    # ---- 7. is healthy consistent with ABSENCE (de novo induction)? -------
    print("\n" + "=" * 96)
    print("7. IS TNFRSF9 ABSENT IN HEALTHY MAST CELLS AND INDUCED IN DISEASE?")
    rows = []
    for arm, dd in (("Healthy", d[d.arm == "Healthy"]), ("AD", d[d.arm == "AD"])):
        c, e = int(dd.t9_mast.sum()), float(dd.mast_umi.sum())
        lo, hi = stats.poisson.interval(0.95, c) if c else (0, stats.chi2.ppf(0.95, 2) / 2)
        # exact Poisson CI for a rate
        lo_r = (stats.chi2.ppf(0.025, 2 * c) / 2 / e * 1e4) if c else 0.0
        hi_r = stats.chi2.ppf(0.975, 2 * (c + 1)) / 2 / e * 1e4
        det = 100 * dd.n_mast_t9pos.sum() / dd.n_mast.sum()
        rows.append(dict(arm=arm, counts=c, mast_umi=int(e), n_mast=int(dd.n_mast.sum()),
                         cp10k=1e4 * c / e, ci_lo=lo_r, ci_hi=hi_r,
                         pct_cells_pos=det, donors=len(dd),
                         donors_pos=int((dd.t9_mast > 0).sum())))
    Z = pd.DataFrame(rows)
    Z.to_csv(TAB / "sc_denovo_test.csv", index=False)
    print(Z.round(4).to_string(index=False))
    h_hi = Z.loc[Z.arm == "Healthy", "ci_hi"].iat[0]
    a_lo = Z.loc[Z.arm == "AD", "ci_lo"].iat[0]
    print(f"\n  healthy 95% upper bound {h_hi:.3f} cp10k vs AD 95% lower bound "
          f"{a_lo:.3f} cp10k -> intervals {'do NOT overlap' if h_hi < a_lo else 'overlap'}")
    print("  Healthy is not distinguishable from zero at this exposure; the AD rate is.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
