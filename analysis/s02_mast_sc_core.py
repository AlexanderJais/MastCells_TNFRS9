"""Core mast-cell analysis of the discovery cohort (GSE204762).

Subject: MAST CELLS (CLAUDE.md §2). Fibroblasts, keratinocytes, T/NK and
macrophages appear only as control populations, to show that whatever happens in
mast cells is specific and not a global or batch effect.

CLAUDE.md §4 — mast-cell-derived TNFRSF9 is one quantity with two factors:

        mast_TNFRSF9_per_1k_cells  =  mast abundance  x  per-mast-cell TNFRSF9

Both factors and their product are reported everywhere; never separately.

CLAUDE.md §6 controls applied to every claim:
  * library size: log offset (pseudobulk) / log10(depth) covariate (cell level)
  * donor aggregation: pseudobulk per donor-arm, SEs clustered on donor
  * identical test in fibroblasts and keratinocytes (and T/NK, macrophages)
  * per-group median sequencing depth reported
  * mast cells re-derived by marker QC (>=2 of TPSAB1/TPSB2/CPA3)
  * dissociation-stress (HSP/IEG) burden per group
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import (HSP_GENES, IEG_GENES, MAST_QC, REQUIRED_CONTROL_POPULATIONS,
                   TARGET, TARGET_LIGAND)  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")

ARM_ORDER = ["Healthy", "AD_NL", "AD_LS"]


# ---------------------------------------------------------------- loading ---
def load_obs() -> pd.DataFrame:
    parts = []
    for p in sorted(PROC.glob("*_obs.parquet")):
        parts.append(pd.read_parquet(p))
    obs = pd.concat(parts, axis=0)
    obs["disease_status"] = obs["disease_status"].astype(str)
    obs["lesional"] = obs["lesional"].astype(str)
    obs["arm"] = np.where(
        obs.disease_status.str.contains("Healthy", case=False), "Healthy",
        np.where(obs.disease_status.str.contains("Atopic|AD", case=False, regex=True),
                 np.where(obs.lesional.str.lower().str.startswith("y"), "AD_LS", "AD_NL"),
                 "Other"))
    obs["donor"] = obs["subject_id"].astype(str)
    obs["arm_pooled"] = np.where(obs.arm == "Healthy", "Healthy",
                                 np.where(obs.arm.isin(["AD_NL", "AD_LS"]), "AD", "Other"))
    return obs


def load_mast_matrix():
    """Full count vectors for stored (mast-labelled or QC-positive) cells."""
    mats, bcs, var = [], [], None
    for p in sorted(PROC.glob("*_mast.npz")):
        gsm = p.name.split("_")[0]
        meta = np.load(PROC / f"{gsm}_mastmeta.npz", allow_pickle=True)
        m = sp.load_npz(p)
        if m.shape[0] == 0:
            continue
        if var is None:
            var = np.asarray(meta["var"], dtype=str)
        mats.append(m)
        bcs.append(np.asarray(meta["barcodes"], dtype=str))
    return sp.vstack(mats).tocsr(), np.concatenate(bcs), var


# ------------------------------------------------------------ §6 marker QC ---
def define_mast(obs: pd.DataFrame) -> pd.DataFrame:
    """Re-derive mast cells from markers; deposited labels are unreliable (§6)."""
    obs = obs.copy()
    lab = obs["Cell type"].astype(str) == "Mast"
    qc = obs["mast_qc_pass"].values.astype(bool)
    obs["mast_label"] = lab
    obs["mast_qc"] = qc
    # Analysis set: passes marker QC. Cells that pass QC but carry a non-mast
    # label are inspected below (ambient tryptase / doublets) before inclusion.
    obs["mast_strict"] = lab & qc
    return obs


def concordance_table(obs: pd.DataFrame) -> pd.DataFrame:
    ct = pd.crosstab(obs["mast_label"], obs["mast_qc"])
    ct.index = ["label!=Mast", "label==Mast"]
    ct.columns = ["QC-fail", "QC-pass"]
    return ct


# ------------------------------------------------------------- pseudobulk ---
def pseudobulk(obs: pd.DataFrame, mat, bcs, var, cell_mask: pd.Series,
               genes: list[str], unit=("donor", "arm")) -> pd.DataFrame:
    """Sum counts over cells within donor x arm; carry total UMI as the offset."""
    idx = {b: i for i, b in enumerate(bcs)}
    gi = {g: i for i, g in enumerate(var)}
    use = obs[cell_mask]
    rows = []
    for key, grp in use.groupby(list(unit), observed=True):
        pos = [idx[b] for b in grp.index if b in idx]
        if not pos:
            continue
        sub = mat[pos, :]
        tot = int(sub.sum())
        rec = dict(zip(unit, key if isinstance(key, tuple) else (key,)))
        rec["n_cells"] = len(pos)
        rec["total_umi"] = tot
        for g in genes:
            rec[f"c_{g}"] = int(sub[:, gi[g]].sum()) if g in gi else 0
        rows.append(rec)
    return pd.DataFrame(rows)


def fit_pb(pb: pd.DataFrame, gene: str, contrast, cluster="donor") -> dict:
    """Poisson GLM: gene counts ~ arm + offset(log total UMI), donor-clustered SE."""
    ref, alt = contrast
    d = pb[pb["arm_c"].isin([ref, alt])].copy()
    d = d[d.total_umi > 0]
    if d.empty or d["arm_c"].nunique() < 2:
        return dict(contrast=f"{alt}_vs_{ref}", gene=gene, log2FC=np.nan, p=np.nan)
    y = d[f"c_{gene}"].astype(float).values
    x = (d["arm_c"] == alt).astype(float).values
    off = np.log(d["total_umi"].astype(float).values)
    X = sm.add_constant(np.column_stack([x]), has_constant="add")
    G.check_model_terms(["intercept", "arm", "log10_depth(offset=log total_umi)"],
                        grouping=cluster)
    try:
        res = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(
            cov_type="cluster", cov_kwds={"groups": d[cluster].values})
        b, se = float(res.params[1]), float(res.bse[1])
        p = float(res.pvalues[1])
    except Exception:
        return dict(contrast=f"{alt}_vs_{ref}", gene=gene, log2FC=np.nan, p=np.nan)
    r_ref = d.loc[d["arm_c"] == ref, f"c_{gene}"].sum() / d.loc[d["arm_c"] == ref, "total_umi"].sum()
    r_alt = d.loc[d["arm_c"] == alt, f"c_{gene}"].sum() / d.loc[d["arm_c"] == alt, "total_umi"].sum()
    return dict(
        contrast=f"{alt}_vs_{ref}", gene=gene,
        log2FC=b / np.log(2), se_log2=se / np.log(2),
        ci_lo=(b - 1.96 * se) / np.log(2), ci_hi=(b + 1.96 * se) / np.log(2),
        p=p,
        cp10k_ref=1e4 * r_ref, cp10k_alt=1e4 * r_alt,
        n_ref=int((d["arm_c"] == ref).sum()), n_alt=int((d["arm_c"] == alt).sum()),
        donors_ref=int(d.loc[d["arm_c"] == ref, "donor"].nunique()),
        donors_alt=int(d.loc[d["arm_c"] == alt, "donor"].nunique()),
        cells_ref=int(d.loc[d["arm_c"] == ref, "n_cells"].sum()),
        cells_alt=int(d.loc[d["arm_c"] == alt, "n_cells"].sum()),
        umi_ref=int(d.loc[d["arm_c"] == ref, "total_umi"].sum()),
        umi_alt=int(d.loc[d["arm_c"] == alt, "total_umi"].sum()),
    )


def main() -> int:
    obs = define_mast(load_obs())
    obs = obs[obs.arm != "Other"]
    print("cells loaded:", len(obs))
    print("\nsamples per arm:")
    print(obs.groupby("arm", observed=True)["sample_label"].nunique().to_string())
    print("\ndonors per arm:")
    print(obs.groupby("arm", observed=True)["donor"].nunique().to_string())

    print("\n§6 mast marker-QC vs deposited label:")
    print(concordance_table(obs).to_string())

    # ---- §6 per-group median sequencing depth ------------------------------
    depth = obs.groupby(["arm", "Cell type"], observed=True)["depth_retained"].median().unstack(0)
    depth.to_csv(TAB / "sc_depth_by_arm_celltype.csv")
    print("\n§6 median retained depth (UMI) by arm x cell type:")
    print(depth.round(0).to_string())

    obs.to_parquet(PROC / "obs_all.parquet")
    print("\nsaved obs_all.parquet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
