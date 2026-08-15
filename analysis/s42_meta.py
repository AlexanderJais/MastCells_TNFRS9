"""Meta-analysis across the three single-cell cohorts (CLAUDE.md §7).

The discovery cohort gives a large effect, and neither replication cohort
reproduces its magnitude. Rather than declare a winner, all three are combined
on the same scale so that the pooled estimate and its heterogeneity are visible
together.

Model: log rate ratio of TNFRSF9 per mast UMI, healthy versus AD, per cohort,
with Poisson variance on the counts (1/c_healthy + 1/c_AD, Haldane-Anscombe 0.5
correction for zero cells), combined by inverse-variance weighting under both
fixed- and random-effects (DerSimonian-Laird) assumptions. Heterogeneity is
reported as I-squared, because it is the point.

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

import guardrails as G  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])


def log_rr(c_h, e_h, c_a, e_a, corr=0.5):
    """log rate ratio (AD vs healthy) with Haldane-Anscombe correction."""
    ch, ca = c_h + corr, c_a + corr
    lrr = np.log(ca / e_a) - np.log(ch / e_h)
    se = np.sqrt(1 / ch + 1 / ca)
    return lrr, se


def main() -> int:
    disc = pd.read_csv(TAB / "sc_denovo_test.csv").set_index("arm")
    rep = pd.read_csv(TAB / "replication_summary.csv")

    rows = [dict(cohort="Discovery (GSE204762, 3')",
                 donors_h=int(disc.loc["Healthy", "donors"]),
                 donors_a=int(disc.loc["AD", "donors"]),
                 mast_h=int(disc.loc["Healthy", "n_mast"]),
                 mast_a=int(disc.loc["AD", "n_mast"]),
                 c_h=int(disc.loc["Healthy", "counts"]),
                 e_h=float(disc.loc["Healthy", "mast_umi"]),
                 c_a=int(disc.loc["AD", "counts"]),
                 e_a=float(disc.loc["AD", "mast_umi"]))]
    for r in rep.itertuples():
        rows.append(dict(cohort=r.cohort, donors_h=r.donors_healthy,
                         donors_a=r.donors_ad, mast_h=r.mast_healthy,
                         mast_a=r.mast_ad, c_h=r.t9_healthy, e_h=float(r.umi_healthy),
                         c_a=r.t9_ad, e_a=float(r.umi_ad)))
    M = pd.DataFrame(rows)
    M["rate_h"] = 1e4 * M.c_h / M.e_h
    M["rate_a"] = 1e4 * M.c_a / M.e_a
    lrr, se = zip(*[log_rr(r.c_h, r.e_h, r.c_a, r.e_a) for r in M.itertuples()])
    M["log2FC"] = np.array(lrr) / np.log(2)
    M["se_log2"] = np.array(se) / np.log(2)
    M["ci_lo"] = M.log2FC - 1.96 * M.se_log2
    M["ci_hi"] = M.log2FC + 1.96 * M.se_log2

    y, v = np.array(lrr), np.array(se) ** 2
    w = 1 / v
    fe = float((w * y).sum() / w.sum())
    fe_se = float(np.sqrt(1 / w.sum()))
    Q = float((w * (y - fe) ** 2).sum())
    df = len(y) - 1
    I2 = max(0.0, 100 * (Q - df) / Q) if Q > 0 else 0.0
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum()))
    wr = 1 / (v + tau2)
    re = float((wr * y).sum() / wr.sum())
    re_se = float(np.sqrt(1 / wr.sum()))

    print("PER-COHORT (TNFRSF9 per 10k mast UMI, healthy vs AD)")
    print(M[["cohort", "donors_h", "donors_a", "mast_h", "mast_a", "c_h", "c_a",
             "rate_h", "rate_a", "log2FC", "ci_lo", "ci_hi"]].round(3).to_string(index=False))

    print("\nPOOLED")
    print(f"  fixed effects  : log2FC {fe/np.log(2):+.2f} "
          f"(95% CI {(fe-1.96*fe_se)/np.log(2):+.2f} to {(fe+1.96*fe_se)/np.log(2):+.2f}), "
          f"P = {2*stats.norm.sf(abs(fe/fe_se)):.4f}")
    print(f"  random effects : log2FC {re/np.log(2):+.2f} "
          f"(95% CI {(re-1.96*re_se)/np.log(2):+.2f} to {(re+1.96*re_se)/np.log(2):+.2f}), "
          f"P = {2*stats.norm.sf(abs(re/re_se)):.4f}")
    print(f"  heterogeneity  : Q = {Q:.2f} (df {df}), I2 = {I2:.0f}%, "
          f"P_Q = {stats.chi2.sf(Q, df):.3f}")

    M.loc[len(M)] = {**{c: np.nan for c in M.columns}, "cohort": "POOLED (fixed)",
                     "log2FC": fe / np.log(2), "se_log2": fe_se / np.log(2),
                     "ci_lo": (fe - 1.96 * fe_se) / np.log(2),
                     "ci_hi": (fe + 1.96 * fe_se) / np.log(2)}
    M.loc[len(M)] = {**{c: np.nan for c in M.columns}, "cohort": "POOLED (random)",
                     "log2FC": re / np.log(2), "se_log2": re_se / np.log(2),
                     "ci_lo": (re - 1.96 * re_se) / np.log(2),
                     "ci_hi": (re + 1.96 * re_se) / np.log(2)}
    M.to_csv(TAB / "meta_analysis.csv", index=False)

    print("\nMAST-CELL RECOVERY, the limiting factor")
    for r in M.dropna(subset=["mast_h"]).itertuples():
        print(f"  {r.cohort:36s} healthy {int(r.mast_h):5d} / AD {int(r.mast_a):5d} mast cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
