"""Replication test (CLAUDE.md §7), using a mast-cell rule calibrated on the
discovery cohort.

Neither replication deposit carries author cell-type labels, and detection-based
marker QC alone is not enough in whole-skin data: ambient tryptase released by
lysed mast cells contaminates keratinocyte droplets, so ">=2 of TPSAB1/TPSB2/
CPA3 detected" tags hundreds of keratinocytes. The rule used here therefore adds
a magnitude requirement, and the threshold is calibrated against the deposited
labels of the discovery cohort rather than chosen by eye:

    mast cell  :=  (>=2 of TPSAB1/TPSB2/CPA3 detected)          [§6 marker QC]
                   AND (tryptase+CPA3 >= 50 transcripts per 10k UMI)

Against the 4,313 author-labelled mast cells of GSE204762 this rule has
precision 96.0% and sensitivity 79.3%; loosening or tightening the threshold
trades sensitivity for very little precision (see calibration table printed
below).

The same donor-level exact permutation test as the discovery cohort is applied.
Subject: MAST CELLS (§2).
"""
from __future__ import annotations

import sys
import warnings
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import MAST_QC, TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DISC = ROOT / "data" / "processed" / "GSE204762"
REP = ROOT / "data" / "processed" / "replication"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])

TRYP_THRESHOLD = 50.0   # transcripts per 10k UMI, calibrated below


def tryp_rate(df: pd.DataFrame) -> pd.Series:
    return 1e4 * df[[f"g_{g}" for g in MAST_QC]].sum(axis=1) / df.depth_retained


def calibrate() -> pd.DataFrame:
    o = pd.read_parquet(DISC / "obs_all.parquet")
    o = o[o.arm.isin(["Healthy", "AD_NL", "AD_LS"])]
    lab = o["Cell type"].astype(str) == "Mast"
    r = tryp_rate(o)
    rows = []
    for thr in (0, 25, 50, 100, 200):
        pred = (r >= thr) & o.mast_qc_pass
        tp = int((pred & lab).sum()); fp = int((pred & ~lab).sum())
        fn = int(((~pred) & lab).sum())
        rows.append(dict(threshold_per10k=thr, n_called=int(pred.sum()),
                         sensitivity=100 * tp / max(tp + fn, 1),
                         precision=100 * tp / max(tp + fp, 1)))
    return pd.DataFrame(rows)


def exact_perm(counts, expo, alt):
    nd, k = len(counts), int(alt.sum())
    if k == 0 or k == nd:
        return np.nan, np.nan
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


def rule_sensitivity() -> pd.DataFrame:
    """Does the discovery result depend on WHICH mast-cell rule is applied?

    The discovery cohort carries author labels and uses `label AND >=2 markers`;
    the replication cohorts have no labels and use `>=2 markers AND tryptase
    magnitude`. A replication is only interpretable if the two arms of the
    comparison are not selected by different rules, so the discovery cohort is
    re-analysed here under the replication's rule.
    """
    o = pd.read_parquet(DISC / "obs_all.parquet")
    o = o[o.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    o["armp"] = np.where(o.arm == "Healthy", "Healthy", "AD")
    o["mast_cal"] = o.mast_qc_pass & (tryp_rate(o) >= TRYP_THRESHOLD)
    rules = {
        "discovery rule (label AND >=2 markers)": o.mast_strict.values.astype(bool),
        "replication rule (>=2 markers AND tryptase)": o.mast_cal.values.astype(bool),
        ">=2 markers only": o.mast_qc_pass.values.astype(bool),
    }
    from scipy import stats as sps
    rows = []
    for name, mask in rules.items():
        m = o[mask]
        rec = dict(rule=name)
        for arm in ("Healthy", "AD"):
            g = m[m.armp == arm]
            rec[f"mast_{arm}"] = len(g)
            rec[f"t9pos_{arm}"] = int((g[f"g_{TARGET}"] > 0).sum())
            rec[f"pct_{arm}"] = 100 * (g[f"g_{TARGET}"] > 0).sum() / max(len(g), 1)
        orr, p = sps.fisher_exact(
            [[rec["t9pos_AD"], rec["mast_AD"] - rec["t9pos_AD"]],
             [rec["t9pos_Healthy"], rec["mast_Healthy"] - rec["t9pos_Healthy"]]])
        rec.update(odds_ratio=orr, p_fisher=p)
        rows.append(rec)
    return pd.DataFrame(rows)


def run(obs: pd.DataFrame, cohort: str) -> tuple[pd.DataFrame, dict]:
    obs = obs.copy()
    obs["tryp_rate"] = tryp_rate(obs)
    obs["mast_cal"] = obs.mast_qc_pass & (obs.tryp_rate >= TRYP_THRESHOLD)
    per = []
    for donor, g in obs.groupby("donor", observed=True):
        m = g[g.mast_cal]
        per.append(dict(cohort=cohort, donor=donor, arm=g.arm.iat[0],
                        n_cells=len(g), n_mast=len(m),
                        mast_pct=100 * len(m) / len(g),
                        mast_umi=int(m.depth_retained.sum()),
                        t9_mast=int(m[f"g_{TARGET}"].sum()) if len(m) else 0,
                        n_mast_t9pos=int((m[f"g_{TARGET}"] > 0).sum()) if len(m) else 0,
                        mast_median_depth=float(m.depth_retained.median()) if len(m) else np.nan))
    P = pd.DataFrame(per)
    alt = (P.arm == "AD").values
    usable = P.mast_umi.sum() > 0 and P.loc[alt, "mast_umi"].sum() > 0 and P.loc[~alt, "mast_umi"].sum() > 0
    if usable:
        lf, p = exact_perm(P.t9_mast.values.astype(float),
                           P.mast_umi.values.astype(float), alt)
    else:
        lf, p = np.nan, np.nan
    summ = dict(
        cohort=cohort, donors_healthy=int((~alt).sum()), donors_ad=int(alt.sum()),
        mast_healthy=int(P.loc[~alt, "n_mast"].sum()), mast_ad=int(P.loc[alt, "n_mast"].sum()),
        umi_healthy=int(P.loc[~alt, "mast_umi"].sum()), umi_ad=int(P.loc[alt, "mast_umi"].sum()),
        t9_healthy=int(P.loc[~alt, "t9_mast"].sum()), t9_ad=int(P.loc[alt, "t9_mast"].sum()),
        rate_healthy=1e4 * P.loc[~alt, "t9_mast"].sum() / max(P.loc[~alt, "mast_umi"].sum(), 1),
        rate_ad=1e4 * P.loc[alt, "t9_mast"].sum() / max(P.loc[alt, "mast_umi"].sum(), 1),
        log2FC=lf, p_perm=p)
    return P, summ


def main() -> int:
    print("CALIBRATION of the mast-cell rule on the discovery cohort (author labels):")
    cal = calibrate()
    print(cal.round(1).to_string(index=False))
    cal.to_csv(TAB / "replication_mast_rule_calibration.csv", index=False)
    print(f"\nusing threshold {TRYP_THRESHOLD:.0f} transcripts per 10k UMI\n")

    print("SENSITIVITY — the discovery cohort under each mast-cell rule")
    print("(the replication cannot be read as a failure to replicate if the two")
    print(" cohorts select their mast cells by different rules)")
    RS = rule_sensitivity()
    RS.to_csv(TAB / "replication_rule_sensitivity.csv", index=False)
    print(RS.round(3).to_string(index=False))
    print()

    frames, summaries = [], []
    for tag, f in (("REP1 (GSE222840+GSE173205, 5')", "rep1_obs.parquet"),
                   ("REP2 (GSE153760 biopsies, 3' v3)", "rep2_obs.parquet")):
        obs = pd.read_parquet(REP / f)
        P, s = run(obs, tag)
        frames.append(P); summaries.append(s)
        print("=" * 96)
        print(tag)
        print(P[["donor", "arm", "n_cells", "n_mast", "mast_pct", "mast_umi",
                 "t9_mast", "n_mast_t9pos", "mast_median_depth"]].round(2).to_string(index=False))

    S = pd.DataFrame(summaries)
    pd.concat(frames).to_csv(TAB / "replication_donor_level.csv", index=False)
    S.to_csv(TAB / "replication_summary.csv", index=False)
    print("\n" + "=" * 96)
    print("REPLICATION SUMMARY — same test as discovery (donor-level exact permutation)")
    print(S.round(4).to_string(index=False))

    print("\nDiscovery for comparison: 6 vs 11 donors, 538 vs 3,123 mast cells,")
    print("  3 vs 107 transcripts, rate 0.057 vs 0.469 per 10k, log2FC +2.81, p=0.012")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
