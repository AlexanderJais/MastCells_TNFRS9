"""In-vitro test: is TNFRSF9 induced when primary human mast cells are given an
AD-like inflammatory stimulus?

The tissue analyses cannot separate "mast cells transcribe more TNFRSF9 in AD"
from "AD skin contains a different mast-cell population". A controlled
stimulation experiment can, and two public ones do exactly the proposed
experiment:

  GSE196862 — primary human SKIN mast cells (HSMCs), the correct cell type and
      tissue, stimulated for the AD milieu: IgE receptor cross-linking, and the
      epithelial alarmins IL-33, TSLP and IL-25, alone and in combination.
      Cufflinks FPKM, 24 libraries.

  GSE235240 — primary human mast cells co-cultured with autologous effector/
      memory CD4+ T cells (resting or anti-CD3/CD28-activated), with IL-33 and
      IgE/antigen as reference stimuli. 4 donors, paired design, normalised
      counts.

Positive controls (genes that MUST move if the stimulation worked) and
housekeeping negative controls are tested alongside TNFRSF9, so that a null
result can be distinguished from a failed experiment.

Subject: MAST CELLS (CLAUDE.md §2).
"""
from __future__ import annotations

import gzip
import re
import sys
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import TARGET, TARGET_LIGAND  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")

# genes that must respond if an activation stimulus worked
POSITIVE = ["TNF", "IL6", "CXCL8", "CCL1", "NR4A1", "EGR1", "IL13", "CSF2",
            "PTGS2", "IL1RL1"]
NEGATIVE = ["ACTB", "GAPDH", "B2M", "RPL13A", "TMSB4X"]
MASTID = ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"]
PANEL = [TARGET, TARGET_LIGAND, "TNFRSF4", "TNFRSF18"] + POSITIVE + NEGATIVE + MASTID


# ------------------------------------------------------------- GSE196862 ---
def load_hsmc() -> pd.DataFrame:
    tar = RAW / "GSE196862_RAW.tar"
    out = {}
    with tarfile.open(tar) as t:
        for m in t.getmembers():
            if not m.name.endswith(".fpkm_tracking.gz"):
                continue
            label = re.sub(r"^GSM\d+_HSMCs_", "", m.name).replace(".fpkm_tracking.gz", "")
            with gzip.open(t.extractfile(m), "rt") as fh:
                df = pd.read_csv(fh, sep="\t")
            gcol = "gene_short_name" if "gene_short_name" in df.columns else "tracking_id"
            s = df.groupby(gcol)["FPKM"].max()
            out[label] = s
    X = pd.DataFrame(out)
    X.columns = [c for c in X.columns]
    return X


def hsmc_report(X: pd.DataFrame) -> pd.DataFrame:
    cond = pd.Series({c: re.sub(r"_\d+$", "", c) for c in X.columns})
    order = ["UN", "IgE", "IL33", "TSLP", "IL25", "IL33IL25TSLP",
             "IgEIL33", "IgETSLP", "IgEIL25", "IgEIL33IL25TSLP"]
    order = [o for o in order if (cond == o).any()]
    rows = []
    for g in PANEL:
        if g not in X.index:
            continue
        r = {"gene": g}
        for o in order:
            v = X.loc[g, cond[cond == o].index].astype(float)
            r[o] = float(v.mean())
        rows.append(r)
    R = pd.DataFrame(rows).set_index("gene")
    un = R["UN"].replace(0, np.nan)
    F = np.log2((R.add(0.01)).div(un + 0.01, axis=0))
    return R, F, order


# ------------------------------------------------------------- GSE235240 ---
def load_teff() -> pd.DataFrame:
    return pd.read_csv(RAW / "GSE235240_counts.txt.gz", sep="\t", index_col=0)


def teff_report(X: pd.DataFrame):
    cond = pd.Series({c: c.rsplit("_", 1)[0] for c in X.columns})
    donor = pd.Series({c: c.rsplit("_", 1)[1] for c in X.columns})
    conds = ["Ctl", "LT", "LT_B", "IL33", "IgE"]
    conds = [c for c in conds if (cond == c).any()]
    rows = []
    for g in PANEL:
        if g not in X.index:
            continue
        base = X.loc[g, cond[cond == "Ctl"].index]
        base.index = [donor[c] for c in base.index]
        r = {"gene": g, "Ctl": float(base.mean())}
        for c in conds:
            if c == "Ctl":
                continue
            v = X.loc[g, cond[cond == c].index]
            v.index = [donor[i] for i in v.index]
            common = base.index.intersection(v.index)
            lf = np.log2((v[common] + 1) / (base[common] + 1))
            r[c] = float(v.mean())
            r[f"log2FC_{c}"] = float(lf.mean())
            r[f"p_{c}"] = float(stats.ttest_rel(np.log2(v[common] + 1),
                                                np.log2(base[common] + 1)).pvalue) \
                if len(common) > 1 else np.nan
            r[f"n_{c}"] = int(len(common))
        rows.append(r)
    return pd.DataFrame(rows).set_index("gene")


def main() -> int:
    print("=" * 100)
    print("GSE196862 — PRIMARY HUMAN SKIN MAST CELLS + AD ALARMINS / IgE")
    X = load_hsmc()
    print(f"  {X.shape[1]} libraries, {X.shape[0]} genes")
    R, F, order = hsmc_report(X)
    print("\n  mean FPKM by condition:")
    print(R.reindex([TARGET, TARGET_LIGAND] + POSITIVE[:5] + MASTID[:3] + NEGATIVE[:2])
          .round(2).to_string())
    print("\n  log2 fold change vs unstimulated:")
    print(F.reindex([TARGET, TARGET_LIGAND] + POSITIVE[:5] + MASTID[:3] + NEGATIVE[:2])
          .drop(columns=["UN"]).round(2).to_string())
    R.to_csv(TAB / "invitro_GSE196862_fpkm.csv")
    F.to_csv(TAB / "invitro_GSE196862_log2fc.csv")

    # per-replicate test for the target: unstimulated vs each stimulus
    cond = pd.Series({c: re.sub(r"_\d+$", "", c) for c in X.columns})
    un = X.loc[TARGET, cond[cond == "UN"].index].astype(float)
    print(f"\n  {TARGET} per replicate:")
    for o in order:
        v = X.loc[TARGET, cond[cond == o].index].astype(float)
        print(f"    {o:18s} n={len(v)}  FPKM {np.round(v.values,3)}  mean {v.mean():.3f}")
    best = F.loc[TARGET].drop("UN").astype(float)
    print(f"\n  strongest induction: {best.idxmax()} "
          f"(log2FC {best.max():+.2f}, {2**best.max():.1f}x)")

    print("\n" + "=" * 100)
    print("GSE235240 — PRIMARY HUMAN MAST CELLS + ACTIVATED CD4+ T CELLS (4 donors, paired)")
    Y = load_teff()
    print(f"  {Y.shape[1]} libraries, {Y.shape[0]} genes")
    T = teff_report(Y)
    cols = ["Ctl"] + [c for c in ["LT", "LT_B", "IL33", "IgE"] if c in T.columns]
    show = [c for c in T.columns if c.startswith("log2FC_") or c.startswith("p_")]
    print("\n  normalised counts by condition:")
    print(T.reindex([TARGET, TARGET_LIGAND] + POSITIVE[:4] + MASTID[:3])[cols]
          .round(1).to_string())
    print("\n  paired log2FC vs control (n donors in n_ columns):")
    print(T.reindex([TARGET, TARGET_LIGAND] + POSITIVE[:4] + MASTID[:3])[show]
          .round(3).to_string())
    T.to_csv(TAB / "invitro_GSE235240.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
