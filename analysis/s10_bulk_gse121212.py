"""Bulk arm — GSE121212 (Tsoi et al.): whole-skin RNA-seq, one study, one pipeline.

38 healthy controls, 21 AD lesional and 27 AD non-lesional biopsies, with
psoriasis lesional/non-lesional carried as a disease-specificity control.

CLAUDE.md §4 governs the interpretation: in bulk tissue, mast-cell abundance and
per-mast-cell expression are *not separable*. Whole-skin TNFRSF9 is the sum over
every cell type of (abundance x per-cell expression). This script therefore
reports three linked quantities and never presents them as separate questions:

  1. whole-skin TNFRSF9
  2. mast-cell content (tryptase/CPA3 transcripts, the bulk proxy for abundance)
  3. TNFRSF9 per unit mast-cell content

Model: Poisson GLM on raw counts with log(library size) offset (§6 depth
control) and cluster-robust standard errors on patient (§6 no pseudo-
replication: AD lesional and non-lesional come from the same individual).
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import MAST_QC  # noqa: E402
from palette import DOUBLE_COL, GROUP_COLORS, MUTED, SINGLE_COL, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE121212_readcount.txt.gz"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")

# Genes reported alongside the target: mast content, and the §6 control
# populations that must show the same test.
CONTROL_GENES = {
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN"],
    "Keratinocytes": ["KRT14", "KRT5", "KRT10"],
    "T/NK": ["CD3D", "CD3E", "TRAC"],
    "Macrophages": ["CD68", "CD163", "AIF1"],
}
MAST_GENES = MAST_QC + ["MS4A2", "CMA1", "KIT", "HDC"]


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(RAW, sep="\t", index_col=0)
    df = df.loc[:, ~df.columns.duplicated()]
    meta = []
    for c in df.columns:
        m = re.match(r"^(AD|PSO|CTRL)_(\d+)_(.*)$", c)
        if not m:
            continue
        dis, pid, tis = m.group(1), m.group(2), m.group(3).strip()
        if dis == "CTRL":
            grp, les = "Healthy", "Healthy"
        elif dis == "AD":
            les = "LS" if "lesion" in tis and "non" not in tis else "NL"
            grp = "AD"
        else:
            les = "LS" if "lesion" in tis and "non" not in tis else "NL"
            grp = "PSO"
        meta.append(dict(sample=c, group=grp, lesional=les,
                         patient=f"{dis}_{pid}", raw_tissue=tis))
    meta = pd.DataFrame(meta).set_index("sample")
    df = df[meta.index]
    return df, meta


def fit(counts: pd.Series, meta: pd.DataFrame, lib: pd.Series,
        contrast: tuple[str, str], key: str = "arm") -> dict:
    """Poisson GLM, log-library-size offset, cluster-robust SE on patient."""
    ref, alt = contrast
    keep = meta[key].isin([ref, alt])
    m = meta[keep]
    y = counts[keep.index[keep]].astype(float).values
    x = (m[key] == alt).astype(float).values
    off = np.log(lib[keep.index[keep]].astype(float).values)
    X = sm.add_constant(np.column_stack([x]), has_constant="add")
    G.check_model_terms(["intercept", key, "log10_depth(offset log library size)"],
                        grouping="patient")
    res = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(
        cov_type="cluster", cov_kwds={"groups": m["patient"].values})
    beta = float(res.params[1])
    se = float(res.bse[1])
    return dict(
        contrast=f"{alt}_vs_{ref}", log2FC=beta / np.log(2),
        se_log2=se / np.log(2), z=beta / se if se else np.nan,
        p=float(res.pvalues[1]),
        ci_lo=(beta - 1.96 * se) / np.log(2), ci_hi=(beta + 1.96 * se) / np.log(2),
        n_ref=int((~(m[key] == alt)).sum()), n_alt=int((m[key] == alt).sum()),
        n_patients=int(m["patient"].nunique()),
        cpm_ref=float(np.median(1e6 * counts[m.index[m[key] == ref]] / lib[m.index[m[key] == ref]])),
        cpm_alt=float(np.median(1e6 * counts[m.index[m[key] == alt]] / lib[m.index[m[key] == alt]])),
    )


def main() -> int:
    df, meta = load()
    lib = df.sum(axis=0)
    meta["lib"] = lib

    # arms for the contract's comparison hierarchy (§3)
    meta["arm"] = np.where(
        meta.group == "Healthy", "Healthy",
        np.where(meta.group == "AD", "AD_" + meta.lesional, "PSO_" + meta.lesional))
    meta["arm_pooled"] = np.where(meta.group == "Healthy", "Healthy",
                                  np.where(meta.group == "AD", "AD", "PSO"))

    print("cohort composition:")
    print(meta.groupby(["group", "lesional"], observed=True).size().to_string())
    print("\nper-group median library size (§6):")
    depth = meta.groupby("arm_pooled", observed=True)["lib"].median()
    print(depth.to_string())
    depth.rename("median_library_size").to_frame().to_csv(TAB / "bulk_GSE121212_depth.csv")

    # ---- §3 comparison hierarchy, in order --------------------------------
    order = G.comparison_order([G.PRIMARY, *G.SECONDARY, G.REFINEMENT])
    specs = {
        G.PRIMARY:      ("arm_pooled", ("Healthy", "AD")),
        "Healthy_vs_AD_NL": ("arm", ("Healthy", "AD_NL")),
        "Healthy_vs_AD_LS": ("arm", ("Healthy", "AD_LS")),
        G.REFINEMENT:   ("arm", ("AD_NL", "AD_LS")),
    }
    # disease-specificity control (psoriasis) — reported, never the headline
    extra = {"Healthy_vs_PSO": ("arm_pooled", ("Healthy", "PSO"))}

    genes = ["TNFRSF9", "TNFSF9"] + MAST_GENES + \
            [g for v in CONTROL_GENES.values() for g in v]
    genes = [g for g in dict.fromkeys(genes) if g in df.index]

    rows = []
    for comp in order + list(extra):
        key, contrast = {**specs, **extra}[comp]
        for g in genes:
            r = fit(df.loc[g], meta, lib, contrast, key=key)
            r.update(gene=g, comparison=comp)
            rows.append(r)
    res = pd.DataFrame(rows)

    # BH across genes within each comparison
    from statsmodels.stats.multitest import multipletests
    res["fdr"] = np.nan
    for comp, idx in res.groupby("comparison").groups.items():
        res.loc[idx, "fdr"] = multipletests(res.loc[idx, "p"], method="fdr_bh")[1]
    res.to_csv(TAB / "bulk_GSE121212_glm.csv", index=False)

    show = res[res.gene.isin(["TNFRSF9", "TPSAB1", "TPSB2", "CPA3", "COL1A1",
                              "KRT14", "CD3D", "CD68"])]
    print("\n" + "=" * 100)
    print(show[["comparison", "gene", "log2FC", "ci_lo", "ci_hi", "p", "fdr",
                "cpm_ref", "cpm_alt", "n_ref", "n_alt"]].to_string(index=False))

    # ---- §4 abundance x per-cell expression, as one quantity ---------------
    cpm = 1e6 * df / lib
    mast_content = cpm.loc[[g for g in MAST_QC if g in cpm.index]].sum(axis=0)
    t9 = cpm.loc["TNFRSF9"]
    out = meta.copy()
    out["TNFRSF9_cpm"] = t9
    out["mast_content_cpm"] = mast_content
    out["TNFRSF9_per_mast"] = t9 / mast_content.replace(0, np.nan)
    out.to_csv(TAB / "bulk_GSE121212_per_sample.csv")
    print("\n§4 decomposition (median per arm):")
    print(out.groupby("arm", observed=True)[
        ["TNFRSF9_cpm", "mast_content_cpm", "TNFRSF9_per_mast"]].median().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
