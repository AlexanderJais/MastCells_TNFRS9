"""Replication (CLAUDE.md §7) — two independent AD cohorts.

A finding from one cohort is a hypothesis. The discovery result (GSE204762) is
tested here in cohorts that share no samples with it, and what fails to
replicate is reported in the same place as what succeeds.

  REP1  "Vienna 5'"  — GSE222840 (5 AD lesional biopsies) + GSE173205
        (4 healthy control biopsies). Same laboratory, same 10x 5' protocol,
        but the healthy arm comes from a companion series: a cross-series
        comparator, disclosed here as CLAUDE.md §5 requires.
  REP2  "GSE153760"  — 4 AD and 2 healthy *biopsies* from a single study and a
        single 3' v3 run. Blister samples from the same study are excluded:
        the authors themselves report that suction blisters do not recover mast
        cells, which is the population under study.

Neither deposit carries author cell-type labels, so mast cells are derived here:
Leiden clustering, cluster annotation by lineage markers, then the §6 marker QC
(>=2 of TPSAB1/TPSB2/CPA3) applied within the mast cluster.
"""
from __future__ import annotations

import gzip
import shutil
import sys
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import (CONTROL_MARKERS, HSP_GENES, IEG_GENES, MAST_IDENTITY,
                   MAST_QC, PANEL, TARGET)  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed" / "replication"
TAB = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")

REP1_AD = {  # GSE222840, atopic dermatitis lesional biopsies
    "GSM6932986": "P74_AD1", "GSM6932987": "P75_AD2", "GSM6932988": "P77_AD3",
    "GSM6932989": "P81_AD4", "GSM6932990": "P96_AD5",
}
REP1_HC = {  # GSE173205, healthy control biopsies, same laboratory
    "GSM5534590": "P112_HC", "GSM5534591": "P115_HC",
    "GSM5534592": "P116_HC", "GSM5534593": "P121_HC",
}
REP2 = {  # GSE153760 biopsies only
    "GSM4653859": ("AD5", "AD"), "GSM4653860": ("AD6", "AD"),
    "GSM4653861": ("AD7", "AD"), "GSM4653862": ("AD8", "AD"),
    "GSM4653868": ("HC6", "Healthy"), "GSM4653869": ("HC7", "Healthy"),
}


def untar(series: str) -> Path:
    d = RAW / series
    d.mkdir(exist_ok=True)
    tar = RAW / f"{series}_RAW.tar"
    if tar.exists() and not any(d.glob("*matrix.mtx.gz")):
        with tarfile.open(tar) as t:
            t.extractall(d)
    return d


def read_triplet(d: Path, gsm: str):
    import scanpy as sc
    mtx = next(d.glob(f"{gsm}_*matrix.mtx.gz"))
    stem = mtx.name.replace("matrix.mtx.gz", "")
    tmp = OUT / "tmp" / gsm
    tmp.mkdir(parents=True, exist_ok=True)
    for src, dst in ((f"{stem}matrix.mtx.gz", "matrix.mtx.gz"),
                     (f"{stem}barcodes.tsv.gz", "barcodes.tsv.gz"),
                     (f"{stem}features.tsv.gz", "features.tsv.gz")):
        s = d / src
        if not s.exists():
            s = next(d.glob(f"{gsm}_*{dst.split('.')[0]}*"))
        shutil.copyfile(s, tmp / dst)
    a = sc.read_10x_mtx(tmp, var_names="gene_symbols", make_unique=True)
    shutil.rmtree(tmp)
    return a


def build(cohort: str, samples: list[tuple[str, Path, str, str, str]]):
    """samples: (gsm, dir, label, arm, donor)"""
    import scanpy as sc
    import anndata as ad
    parts = []
    for gsm, d, label, arm, donor in samples:
        a = read_triplet(d, gsm)
        a.var_names_make_unique()
        a.obs["sample_label"] = label
        a.obs["arm"] = arm
        a.obs["donor"] = donor
        a.obs["gsm"] = gsm
        a.obs_names = [f"{label}_{b}" for b in a.obs_names]
        sc.pp.filter_cells(a, min_genes=200)
        mt = [g for g in a.var_names if g.startswith("MT-")]
        a.obs["mt_frac"] = (np.asarray(a[:, mt].X.sum(1)).ravel()
                            / np.maximum(np.asarray(a.X.sum(1)).ravel(), 1)) if mt else 0.0
        a = a[a.obs.mt_frac < 0.20].copy()
        print(f"    {label} ({arm}): {a.n_obs} cells")
        parts.append(a)
    A = ad.concat(parts, join="inner", label="batch", keys=[p.obs.sample_label[0] for p in parts])
    print(f"  {cohort}: {A.n_obs} cells x {A.n_vars} genes")

    # Extract every per-cell quantity from the raw counts BEFORE any gene
    # subsetting: depth, panel counts and the §6 marker QC must come from the
    # full matrix. Scaling the full matrix would densify it (>14 GB here), so
    # clustering runs on the highly-variable subset only.
    C = sp.csr_matrix(A.X)
    idx = {g: i for i, g in enumerate(A.var_names)}
    counts_info = {"depth_retained": np.asarray(C.sum(1)).ravel().astype(np.int64)}
    for g in dict.fromkeys(PANEL + MAST_IDENTITY):
        counts_info[f"g_{g}"] = (np.asarray(C[:, idx[g]].todense()).ravel().astype(np.int32)
                                 if g in idx else np.zeros(A.n_obs, np.int32))
    for k, v in counts_info.items():
        A.obs[k] = v

    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key="batch")
    A = A[:, A.var.highly_variable].copy()
    sc.pp.scale(A, max_value=10)
    sc.tl.pca(A, n_comps=30)
    try:
        import scanpy.external as sce
        sce.pp.harmony_integrate(A, "batch")
        rep = "X_pca_harmony"
    except Exception:
        rep = "X_pca"
    sc.pp.neighbors(A, n_neighbors=15, use_rep=rep)
    sc.tl.leiden(A, resolution=1.0, key_added="leiden", flavor="igraph", n_iterations=2)
    sc.tl.umap(A)
    A.obs["umap1"] = A.obsm["X_umap"][:, 0]
    A.obs["umap2"] = A.obsm["X_umap"][:, 1]
    return A


def annotate(A) -> pd.DataFrame:
    """Score clusters on lineage markers; the mast cluster is the tryptase one.

    Counts-derived columns were attached to .obs in build() before the
    highly-variable-gene subset, so they reflect the full matrix.
    """
    obs = A.obs.copy()
    qc = np.vstack([obs[f"g_{g}"].values > 0 for g in MAST_QC]).sum(0)
    obs["mast_markers_n"] = qc
    obs["mast_qc_pass"] = qc >= 2

    sets = {"Mast": MAST_QC + ["MS4A2", "KIT", "HPGDS"], **CONTROL_MARKERS}
    score = {}
    for name, gs in sets.items():
        cols = [f"g_{g}" for g in gs if f"g_{g}" in obs.columns]
        r = 1e4 * obs[cols].sum(axis=1) / np.maximum(obs.depth_retained, 1)
        score[name] = np.log1p(r).groupby(obs["leiden"].values).mean()
    S = pd.DataFrame(score)
    obs["cluster_type"] = obs["leiden"].map(S.idxmax(axis=1)).astype(str)
    obs["mast_strict"] = (obs.cluster_type == "Mast") & obs.mast_qc_pass
    return obs, S


def analyse(obs: pd.DataFrame, cohort: str) -> pd.DataFrame:
    from scipy import stats
    rows = []
    per = []
    for (donor, arm), g in obs.groupby(["donor", "arm"], observed=True):
        m = g[g.mast_strict]
        per.append(dict(cohort=cohort, donor=donor, arm=arm, n_cells=len(g),
                        n_mast=len(m), mast_umi=int(m.depth_retained.sum()),
                        tissue_umi=int(g.depth_retained.sum()),
                        t9_mast=int(m[f"g_{TARGET}"].sum()) if len(m) else 0,
                        t9_tissue=int(g[f"g_{TARGET}"].sum()),
                        n_mast_t9pos=int((m[f"g_{TARGET}"] > 0).sum()) if len(m) else 0,
                        mast_median_depth=float(m.depth_retained.median()) if len(m) else np.nan,
                        stress_mast=float(m[[f"g_{x}" for x in HSP_GENES + IEG_GENES
                                             if f"g_{x}" in m.columns]].sum().sum()
                                          / max(m.depth_retained.sum(), 1) * 1e4) if len(m) else np.nan))
    P = pd.DataFrame(per)
    return P


def main() -> int:
    order = G.comparison_order([G.PRIMARY])
    G.check_model_terms(["intercept", "arm", "log10_depth(offset=log mast UMI)"],
                        grouping="donor")

    results = {}
    # ---------------- REP1 -------------------------------------------------
    print("REP1  GSE222840 (AD) + GSE173205 (healthy), same laboratory, 5' 10x")
    d222 = untar("GSE222840")
    d173 = RAW / "GSE173205"
    samples = [(g, d222, lab, "AD", lab.split("_")[0]) for g, lab in REP1_AD.items()]
    samples += [(g, d173, lab, "Healthy", lab.split("_")[0]) for g, lab in REP1_HC.items()]
    A1 = build("REP1", samples)
    obs1, S1 = annotate(A1)
    obs1.to_parquet(OUT / "rep1_obs.parquet")
    A1.obs[["umap1", "umap2"]].to_parquet(OUT / "rep1_umap.parquet")
    print("  cluster annotation (mean log1p rate per cluster, argmax):")
    print(obs1.groupby("cluster_type", observed=True).size().to_string())
    results["REP1"] = analyse(obs1, "REP1")

    # ---------------- REP2 -------------------------------------------------
    print("\nREP2  GSE153760 biopsies only (blisters excluded: no mast cells)")
    d153 = untar("GSE153760")
    samples = [(g, d153, lab, arm, lab) for g, (lab, arm) in REP2.items()]
    A2 = build("REP2", samples)
    obs2, S2 = annotate(A2)
    obs2.to_parquet(OUT / "rep2_obs.parquet")
    A2.obs[["umap1", "umap2"]].to_parquet(OUT / "rep2_umap.parquet")
    print(obs2.groupby("cluster_type", observed=True).size().to_string())
    results["REP2"] = analyse(obs2, "REP2")

    P = pd.concat(results.values())
    P.to_csv(TAB / "replication_sample_level.csv", index=False)
    print("\n" + "=" * 92)
    print(P.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
