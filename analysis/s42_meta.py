"""Meta-analysis across the single-cell cohorts that can measure the target
(CLAUDE.md §7).

QUALIFICATION CRITERION, applied before any TNFRSF9 number is looked at: a cohort
contributes only if it recovered at least 100 marker-QC mast cells in EACH arm
and at least 0.5% of its cells as mast cells. This is a mast-cell QC threshold
judged on tryptase/CPA3 marker data alone; it is independent of TNFRSF9 and of
the outcome.

  Discovery GSE204762  538 / 3,123 mast cells, 1.31% recovery   -> qualifies
  REP2 GSE153760       278 / 1,101 mast cells, 4.41% recovery   -> qualifies
  REP1 GSE222840+GSE173205  27 / 91 mast cells, 0.106% recovery -> EXCLUDED

REP1 recovered mast cells 12-40x less efficiently than the other two and yielded
2 TNFRSF9 transcripts in 111,370 cells. Pooling it would be pooling a failed
experiment: it cannot detect the target, so its contribution is noise, not
evidence. It is reported as an attempted replication that failed on mast-cell
recovery, not as a negative result about TNFRSF9.

Model: log rate ratio of TNFRSF9 per mast UMI, healthy versus AD, per cohort,
with Poisson variance on the counts (1/c_healthy + 1/c_AD, Haldane-Anscombe 0.5
correction for zero cells), combined by inverse-variance weighting under both
fixed- and random-effects (DerSimonian-Laird) assumptions. With two qualifying
cohorts an I-squared estimate is unstable and is reported for completeness only.

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
    disc_s = pd.read_csv(TAB / "sc_sample_level.csv")
    rep = pd.read_csv(TAB / "replication_summary.csv")
    rep_d = pd.read_csv(TAB / "replication_donor_level.csv")
    rep_cells = rep_d.groupby("cohort").n_cells.sum()

    # Both halves of the qualification criterion are applied in code, not only
    # stated in prose: a cohort must clear an absolute count in each arm AND a
    # recovery fraction over the whole cohort. Both are judged on tryptase/CPA3
    # marker data alone, before any TNFRSF9 number is looked at.
    MIN_MAST_PER_ARM = 100
    MIN_MAST_RECOVERY_PCT = 0.5

    rows = [dict(cohort="Discovery (GSE204762, 3')",
                 donors_h=int(disc.loc["Healthy", "donors"]),
                 donors_a=int(disc.loc["AD", "donors"]),
                 mast_h=int(disc.loc["Healthy", "n_mast"]),
                 mast_a=int(disc.loc["AD", "n_mast"]),
                 cells_total=int(disc_s.n_cells.sum()),
                 c_h=int(disc.loc["Healthy", "counts"]),
                 e_h=float(disc.loc["Healthy", "mast_umi"]),
                 c_a=int(disc.loc["AD", "counts"]),
                 e_a=float(disc.loc["AD", "mast_umi"]))]
    for r in rep.itertuples():
        rows.append(dict(cohort=r.cohort, donors_h=r.donors_healthy,
                         donors_a=r.donors_ad, mast_h=r.mast_healthy,
                         mast_a=r.mast_ad, cells_total=int(rep_cells[r.cohort]),
                         c_h=r.t9_healthy, e_h=float(r.umi_healthy),
                         c_a=r.t9_ad, e_a=float(r.umi_ad)))
    M = pd.DataFrame(rows)
    M["recovery_pct"] = 100 * (M.mast_h + M.mast_a) / M.cells_total
    M["qualifies"] = ((M.mast_h >= MIN_MAST_PER_ARM) & (M.mast_a >= MIN_MAST_PER_ARM)
                      & (M.recovery_pct >= MIN_MAST_RECOVERY_PCT))
    excluded = M[~M.qualifies].copy()
    M = M[M.qualifies].reset_index(drop=True)
    for r in excluded.itertuples():
        print(f"EXCLUDED before analysis: {r.cohort} — {int(r.mast_h)}/{int(r.mast_a)} "
              f"mast cells recovered (threshold {MIN_MAST_PER_ARM} per arm) at "
              f"{r.recovery_pct:.3f}% recovery (threshold {MIN_MAST_RECOVERY_PCT}%)")
    # Second exposure: mast CELLS rather than mast UMI. This is the unadjusted
    # ("absolute") estimand — TNFRSF9 molecules captured per mast cell — and it
    # is the conservative one wherever healthy mast cells are the deeper arm.
    M["rate_h"] = 1e4 * M.c_h / M.e_h
    M["rate_a"] = 1e4 * M.c_a / M.e_a
    lrr_c, se_c = zip(*[log_rr(r.c_h, r.mast_h, r.c_a, r.mast_a) for r in M.itertuples()])
    M["log2FC_perCELL"] = np.array(lrr_c) / np.log(2)
    M["ci_lo_perCELL"] = (np.array(lrr_c) - 1.96 * np.array(se_c)) / np.log(2)
    M["ci_hi_perCELL"] = (np.array(lrr_c) + 1.96 * np.array(se_c)) / np.log(2)
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

    # ---- the same pooling on the unadjusted (per-cell) exposure ----------
    # This is the PRIMARY estimand (§7.1 of the report), so it is written into
    # the table alongside the per-UMI pool rather than left in stdout.
    yc, vc = np.array(lrr_c), np.array(se_c) ** 2
    wc = 1 / vc
    fec = float((wc * yc).sum() / wc.sum())
    fec_se = float(np.sqrt(1 / wc.sum()))
    Qc = float((wc * (yc - fec) ** 2).sum())
    I2c = max(0.0, 100 * (Qc - df) / Qc) if Qc > 0 else 0.0
    tau2c = max(0.0, (Qc - df) / (wc.sum() - (wc ** 2).sum() / wc.sum()))
    wrc = 1 / (vc + tau2c)
    rec = float((wrc * yc).sum() / wrc.sum())
    rec_se = float(np.sqrt(1 / wrc.sum()))

    def pooled_row(label, b, b_se, bc, bc_se, i2, i2c):
        return {**{c: np.nan for c in M.columns}, "cohort": label,
                "log2FC": b / np.log(2), "se_log2": b_se / np.log(2),
                "ci_lo": (b - 1.96 * b_se) / np.log(2),
                "ci_hi": (b + 1.96 * b_se) / np.log(2),
                "p": float(2 * stats.norm.sf(abs(b / b_se))),
                "log2FC_perCELL": bc / np.log(2),
                "ci_lo_perCELL": (bc - 1.96 * bc_se) / np.log(2),
                "ci_hi_perCELL": (bc + 1.96 * bc_se) / np.log(2),
                "p_perCELL": float(2 * stats.norm.sf(abs(bc / bc_se))),
                "I2_pct": i2, "I2_pct_perCELL": i2c}

    M["p"] = np.nan
    M["p_perCELL"] = np.nan
    M["I2_pct"] = np.nan
    M["I2_pct_perCELL"] = np.nan
    M.loc[len(M)] = pooled_row("POOLED (fixed)", fe, fe_se, fec, fec_se, I2, I2c)
    M.loc[len(M)] = pooled_row("POOLED (random)", re, re_se, rec, rec_se, I2, I2c)
    M.to_csv(TAB / "meta_analysis.csv", index=False)

    print("\nPOOLED — UNADJUSTED ESTIMAND (TNFRSF9 per mast CELL)")
    print(f"  per-cohort log2FC: " +
          ", ".join(f"{r.cohort.split(' (')[0]} {r.log2FC_perCELL:+.2f}"
                    for r in M.itertuples() if not r.cohort.startswith("POOLED")))
    print(f"  fixed effects  : log2FC {fec/np.log(2):+.2f} "
          f"(95% CI {(fec-1.96*fec_se)/np.log(2):+.2f} to {(fec+1.96*fec_se)/np.log(2):+.2f}), "
          f"P = {2*stats.norm.sf(abs(fec/fec_se)):.4f}")
    print(f"  random effects : log2FC {rec/np.log(2):+.2f} "
          f"(95% CI {(rec-1.96*rec_se)/np.log(2):+.2f} to {(rec+1.96*rec_se)/np.log(2):+.2f}), "
          f"P = {2*stats.norm.sf(abs(rec/rec_se)):.4f}")
    print(f"  heterogeneity  : I2 = {I2c:.0f}%, P_Q = {stats.chi2.sf(Qc, df):.3f}")

    print("\nMAST-CELL DEPTH BY ARM (does the exposure choice matter per cohort?)")
    print("  ratio = mean UMI per mast cell, AD / healthy; <1 means healthy is deeper")
    for r in M.dropna(subset=["mast_h"]).itertuples():
        ratio = (r.e_a / r.mast_a) / (r.e_h / r.mast_h)
        print(f"  {r.cohort:36s} {ratio:.2f}  "
              f"({'healthy deeper -> per-cell is conservative' if ratio < 1 else 'AD deeper -> per-cell is anticonservative'})")

    print("\nMAST-CELL RECOVERY, the limiting factor")
    for r in M.dropna(subset=["mast_h"]).itertuples():
        print(f"  {r.cohort:36s} healthy {int(r.mast_h):5d} / AD {int(r.mast_a):5d} mast cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
