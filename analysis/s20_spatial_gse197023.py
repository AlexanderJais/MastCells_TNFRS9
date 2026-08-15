"""Spatial arm — GSE197023: Visium whole-transcriptome spatial transcriptomics.

10x Visium of 7 AD patients (paired lesional and non-lesional) and 6 healthy
donors, all in one study and one pipeline (CLAUDE.md §5).

Why this arm matters for the mast-cell question: Visium never dissociates the
tissue. The single-cell estimate of per-mast-cell TNFRSF9 is conditional on
which mast cells survive enzymatic digestion (CLAUDE.md §4); dermal mast cells
are notoriously hard to release and are the shallowest-sequenced population in
the scRNA-seq cohort. Visium measures the intact dermis, so it is the check on
that conditioning.

Resolution caveat, stated up front: a 55 um Visium spot holds several cells, so
a spot is not a mast cell. Spot-level tryptase content is used as the in-situ
proxy for mast-cell abundance (CLAUDE.md §4: abundance and expression are one
quantity and in spot data they are not separable either).
"""
from __future__ import annotations

import subprocess
import sys
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

from genes import (CONTROL_MARKERS, HSP_GENES, IEG_GENES, MAST_IDENTITY,
                   MAST_QC, TARGET, TARGET_LIGAND, TYPE2_CONTEXT)  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
SPATIAL = RAW / "GSE197023"
OUT = ROOT / "data" / "processed" / "GSE197023"
OUT.mkdir(parents=True, exist_ok=True)

PANEL = list(dict.fromkeys(
    [TARGET, TARGET_LIGAND] + MAST_QC + MAST_IDENTITY + TYPE2_CONTEXT
    + HSP_GENES + IEG_GENES
    + [g for v in CONTROL_MARKERS.values() for g in v]))


def extract_all() -> list[Path]:
    SPATIAL.mkdir(parents=True, exist_ok=True)
    tar = RAW / "GSE197023_RAW.tar"
    with tarfile.open(tar) as t:
        names = [n for n in t.getnames() if n.endswith(".tar.gz")]
        for n in names:
            if not (SPATIAL / n).exists():
                t.extract(n, SPATIAL)
    out = []
    for gz in sorted(SPATIAL.glob("*.tar.gz")):
        label = gz.name.split("_", 1)[1].replace(".tar.gz", "")
        d = SPATIAL / label
        if not d.exists():
            with tarfile.open(gz) as t:
                t.extractall(SPATIAL)
        out.append(d)
    return out


def read_sample(d: Path) -> pd.DataFrame | None:
    import scanpy as sc
    h5 = d / "filtered_feature_bc_matrix.h5"
    if not h5.exists():
        cand = list(d.glob("*.h5"))
        if not cand:
            return None
        h5 = cand[0]
    a = sc.read_10x_h5(h5)
    a.var_names_make_unique()
    import scipy.sparse as sp
    X = sp.csr_matrix(a.X)
    X.data = np.rint(X.data).astype(np.int32)
    tot = np.asarray(X.sum(1)).ravel().astype(np.int64)
    ng = np.asarray((X > 0).sum(1)).ravel().astype(np.int32)
    gi = {g: i for i, g in enumerate(np.asarray(a.var_names))}
    df = pd.DataFrame(index=np.asarray(a.obs_names))
    df["sample"] = d.name
    df["total_counts"] = tot
    df["n_genes"] = ng
    for g in PANEL:
        df[f"g_{g}"] = np.asarray(X[:, gi[g]].todense()).ravel().astype(np.int32) if g in gi else 0
    # spatial coordinates, when present
    pos = d / "spatial" / "tissue_positions_list.csv"
    if pos.exists():
        p = pd.read_csv(pos, header=None,
                        names=["barcode", "in_tissue", "row", "col", "y", "x"]).set_index("barcode")
        df = df.join(p[["in_tissue", "row", "col", "y", "x"]], how="left")
    return df


def main() -> int:
    dirs = extract_all()
    frames = []
    for d in dirs:
        try:
            f = read_sample(d)
        except Exception as e:
            print(f"  {d.name}: ERROR {e}")
            continue
        if f is None:
            print(f"  {d.name}: no matrix")
            continue
        frames.append(f)
        print(f"  {d.name}: {len(f)} spots, median depth {np.median(f.total_counts):.0f}")
    spots = pd.concat(frames)

    lab = spots["sample"].astype(str)
    spots["donor"] = lab.str.replace(r"_(LS|NL)$", "", regex=True)
    spots["arm"] = np.where(lab.str.startswith("HE"), "Healthy",
                            np.where(lab.str.endswith("_LS"), "AD_LS", "AD_NL"))
    spots["arm_pooled"] = np.where(spots.arm == "Healthy", "Healthy", "AD")
    spots.to_parquet(OUT / "spots.parquet")
    print("\ntotal spots:", len(spots))
    print(spots.groupby(["arm"], observed=True).agg(
        samples=("sample", "nunique"), donors=("donor", "nunique"),
        spots=("total_counts", "size"), median_depth=("total_counts", "median")).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
