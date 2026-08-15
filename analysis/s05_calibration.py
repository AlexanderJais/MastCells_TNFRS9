"""Calibration: is the discovery test well behaved, or does the design invent results?

Section 5 of s04 found that 17.8% of abundance-matched genes reach a nominal
permutation p<=0.05 for Healthy vs AD in mast cells — 3.6x the 5% expected. Two
explanations, with opposite consequences for the TNFRSF9 claim:

  (i)  the test is miscalibrated (donor-level overdispersion, unequal arm sizes),
       in which case a nominal p means very little; or
  (ii) the mast-cell transcriptome genuinely differs between healthy and AD skin
       across many genes, in which case the test is fine but TNFRSF9 must be
       judged against that broad background rather than against 5%.

This script separates them with a TRUE null: arm labels permuted *within the AD
donors only* (5 vs 6), where no disease contrast exists by construction. The
same statistic and the same abundance-matched gene set are used.
"""
from __future__ import annotations

import sys
import warnings
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])


def load_mast():
    loaded = []
    for p in sorted(PROC.glob("*_mast.npz")):
        gsm = p.name.split("_")[0]
        meta = np.load(PROC / f"{gsm}_mastmeta.npz", allow_pickle=True)
        m = sp.load_npz(p)
        if m.shape[0] == 0:
            continue
        loaded.append((m, np.asarray(meta["var"], dtype=str),
                       np.asarray(meta["barcodes"], dtype=str)))
    common = sorted(set.intersection(*[set(v) for _, v, _ in loaded]))
    pos = {g: i for i, g in enumerate(common)}
    mats, bcs = [], []
    for m, v, b in loaded:
        take = [i for i, g in enumerate(v) if g in pos]
        order = np.argsort([pos[v[i]] for i in take])
        mats.append(m[:, np.asarray(take)[order]])
        bcs.append(b)
    return sp.vstack(mats).tocsr(), np.concatenate(bcs), np.array(common)


def perm_all(Cd: np.ndarray, expo: np.ndarray, alt: np.ndarray):
    """Vectorised exact permutation over donors for every gene column."""
    nd, k = len(expo), int(alt.sum())
    masks = np.zeros((comb(nd, k), nd), dtype=bool)
    for i, c in enumerate(combinations(range(nd), k)):
        masks[i, list(c)] = True
    e_alt = masks @ expo
    e_ref = expo.sum() - e_alt
    c_alt = masks @ Cd
    c_ref = Cd.sum(axis=0)[None, :] - c_alt
    S = np.log((c_alt + .5) / e_alt[:, None]) - np.log((c_ref + .5) / e_ref[:, None])
    obs_i = int(np.where((masks == alt[None, :]).all(axis=1))[0][0])
    obs = S[obs_i]
    p = (np.abs(S) >= np.abs(obs)[None, :] - 1e-12).mean(axis=0)
    return obs, p


def main() -> int:
    obs = pd.read_parquet(PROC / "obs_all.parquet")
    m = obs[obs.mast_strict & obs.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    M, B, var = load_mast()
    keep = pd.Index(B).isin(m.index)
    M, B = M[keep], B[keep]
    mm = m.loc[B]

    donors = mm.donor.values
    udon = pd.unique(donors)
    dmap = {u: i for i, u in enumerate(udon)}
    di = np.array([dmap[x] for x in donors])
    expo = np.bincount(di, weights=mm.depth_retained.values.astype(float),
                       minlength=len(udon))
    is_h = np.array([str(u).startswith("Healthy") for u in udon])

    tot = np.asarray(M.sum(0)).ravel()
    t9 = tot[list(var).index(TARGET)]
    cand = np.where((tot >= 0.5 * t9) & (tot <= 2 * t9))[0]
    dense = np.asarray(M[:, cand].todense())
    Cd = np.zeros((len(udon), len(cand)))
    for i in range(len(udon)):
        Cd[i] = dense[di == i].sum(axis=0)
    gname = var[cand]
    i9 = int(np.where(gname == TARGET)[0][0])

    print(f"abundance-matched gene set: {len(cand)} genes "
          f"(TNFRSF9 = {int(t9)} counts in mast cells)")

    # ---- observed contrast: Healthy vs AD --------------------------------
    obs_o, p_o = perm_all(Cd, expo, ~is_h)
    frac_o = float((p_o <= .05).mean())
    rank9 = int((p_o < p_o[i9]).sum()) + 1
    print("\nOBSERVED  Healthy vs AD (6 vs 11 donors)")
    print(f"  genes with p<=0.05: {100*frac_o:.1f}%  ({int((p_o<=.05).sum())}/{len(cand)})")
    print(f"  TNFRSF9: log2FC {obs_o[i9]/np.log(2):.2f}, p_perm {p_o[i9]:.4f}, "
          f"rank {rank9}/{len(cand)} (empirical p = {rank9/len(cand):.3f})")

    # ---- TRUE null: split the AD donors, no disease contrast -------------
    # A single arbitrary split is not a null estimate: because genes are highly
    # correlated, the fraction of "significant" genes varies enormously from one
    # split to the next. The null is the distribution over many random splits.
    ad = np.where(~is_h)[0]
    Cd_ad, expo_ad = Cd[ad], expo[ad]
    rng = np.random.default_rng(20240815)
    fracs, seen = [], set()
    for _ in range(40):
        k = int(rng.integers(5, 7))
        sel = tuple(sorted(rng.choice(len(ad), k, replace=False)))
        if sel in seen:
            continue
        seen.add(sel)
        alt = np.zeros(len(ad), bool)
        alt[list(sel)] = True
        _, p_n = perm_all(Cd_ad, expo_ad, alt)
        fracs.append(float((p_n <= .05).mean()))
    fracs = np.array(fracs)
    print(f"\nTRUE NULL  AD donors split against themselves ({len(fracs)} random splits,"
          f" no disease contrast)")
    print(f"  genes with p<=0.05: median {100*np.median(fracs):.1f}%  "
          f"(IQR {100*np.percentile(fracs,25):.1f}-{100*np.percentile(fracs,75):.1f}%, "
          f"range {100*fracs.min():.1f}-{100*fracs.max():.1f}%)")
    null_frac = float(np.median(fracs))
    print(f"  observed contrast sits at the {100*(fracs < frac_o).mean():.0f}th percentile "
          f"of this null distribution")

    # ---- healthy-only null is not possible (6 donors, 1 informative) -----
    print(f"\nINTERPRETATION")
    print(f"  observed 'significant' fraction : {100*frac_o:.1f}%")
    print(f"  true-null fraction              : {100*null_frac:.1f}%")
    if null_frac < 0.08:
        print("  -> the test is approximately calibrated; the excess in the observed")
        print("     contrast reflects genuine widespread Healthy-vs-AD differences in")
        print("     the mast-cell transcriptome, not a broken test.")
    else:
        print("  -> the test itself is anticonservative; nominal p values are not")
        print("     interpretable and only the empirical rank should be used.")
    print(f"  TNFRSF9 must therefore be read against that background: empirical")
    print(f"  p = {rank9/len(cand):.3f} among abundance-matched genes.")

    out = pd.DataFrame(dict(gene=gname, log2FC=obs_o / np.log(2), p_perm=p_o))
    out["empirical_p"] = out.p_perm.rank(method="min") / len(out)
    out.sort_values("p_perm").to_csv(TAB / "sc_calibration_matched_genes.csv", index=False)
    pd.DataFrame(dict(quantity=["observed_frac_p05", "true_null_frac_p05",
                                "TNFRSF9_log2FC", "TNFRSF9_p_perm",
                                "TNFRSF9_rank", "TNFRSF9_empirical_p", "n_genes"],
                      value=[frac_o, null_frac, obs_o[i9] / np.log(2), p_o[i9],
                             rank9, rank9 / len(cand), len(cand)])
                 ).to_csv(TAB / "sc_calibration_summary.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
