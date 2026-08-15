"""GSE204762 — primary discovery cohort: stream, extract, discard.

Multi-modal skin atlas (Nature Communications 2025 deposit): whole-skin 3' 10x
scRNA-seq of 11 adult atopic dermatitis patients (paired non-lesional and
lesional biopsies), 7 healthy donors and 2 scleroderma donors, all processed
through one pipeline. Raw counts are deposited per sample as h5ad.

The full archive is 15.7 GB, so each sample is downloaded, reduced and deleted
in turn. What is kept per cell:

  * clinical / technical metadata (donor, disease, lesional status, chemistry,
    batch, sex, author cell-type label)
  * retained library size and detected-gene count  (CLAUDE.md §6, depth control)
  * integer counts for the analysis panel (analysis/genes.py)

and, for every cell that is either labelled Mast by the authors or passes the
§6 marker QC independently, the complete count vector is retained so that
mast-cell analyses are not restricted to the panel.
"""
from __future__ import annotations

import gzip
import os
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

from genes import MAST_QC, PANEL  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE204762"
OUT = ROOT / "data" / "processed" / "GSE204762"
OUT.mkdir(parents=True, exist_ok=True)
RAW.mkdir(parents=True, exist_ok=True)

FILELIST = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE204nnn/GSE204762/suppl/filelist.txt"
SAMPLE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{stem}nnn/{gsm}/suppl/{fname}"

KEEP_OBS = [
    "sample_name", "subject_id", "ind_id", "disease_status", "lesional",
    "disease_lesional", "ind_disease_lesional", "sex", "sex_clinical",
    "chemistry_10x", "batch", "channel", "location", "side",
    "global_disease_assessment", "asthma", "ad_family_history",
    "topical_corticosteroids", "systemic_corticosteroids", "dupilumab",
    "n_umis", "n_counts", "n_genes", "mt_frac", "scrublet", "scrublet_score",
    "phase", "Cell type", "Cell type granular", "Compartment", "leiden",
]


def manifest() -> pd.DataFrame:
    import requests
    txt = requests.get(FILELIST, timeout=120).text
    rows = []
    for line in txt.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3 or parts[0] != "File":
            continue
        fname = parts[1]
        if not fname.endswith(".h5ad.gz"):
            continue
        gsm = fname.split("_", 1)[0]
        label = fname.split("_", 1)[1].replace(".h5ad.gz", "")
        rows.append(dict(gsm=gsm, fname=fname, label=label,
                         size=int(parts[3]) if len(parts) > 3 else 0))
    df = pd.DataFrame(rows)
    df = df[~df.label.str.startswith("mouse")].reset_index(drop=True)  # human only
    return df


def download(gsm: str, fname: str, dest: Path) -> Path:
    stem = gsm[:-3]
    url = SAMPLE_URL.format(stem=stem, gsm=gsm, fname=fname)
    gz = dest / fname
    if gz.exists() and gz.stat().st_size > 1_000_000:
        return gz
    for attempt in range(4):
        rc = subprocess.run(
            ["curl", "-sS", "-L", "--fail", "--max-time", "1800", "-o", str(gz), url]
        ).returncode
        if rc == 0 and gz.stat().st_size > 1_000_000:
            return gz
        import time
        time.sleep(2 ** (attempt + 1))
    raise RuntimeError(f"download failed: {url}")


def process(gsm: str, fname: str, label: str) -> str:
    obs_out = OUT / f"{gsm}_obs.parquet"
    mast_out = OUT / f"{gsm}_mast.npz"
    if obs_out.exists() and mast_out.exists():
        return f"{label}: cached"

    gz = download(gsm, fname, RAW)
    h5 = RAW / fname[:-3]
    if not h5.exists():
        with gzip.open(gz, "rb") as fi, open(h5, "wb") as fo:
            shutil.copyfileobj(fi, fo, length=1 << 24)

    import anndata as ad
    a = ad.read_h5ad(h5)
    counts = a.layers["counts"]
    counts = sp.csr_matrix(counts)
    counts.data = np.rint(counts.data).astype(np.int32)

    var = np.asarray(a.var_names)
    gene_pos = {g: i for i, g in enumerate(var)}

    # ---- per-cell depth (retained library size) and complexity --------------
    depth = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    ngene = np.asarray((counts > 0).sum(axis=1)).ravel().astype(np.int32)

    obs = pd.DataFrame(index=np.asarray(a.obs_names))
    for c in KEEP_OBS:
        if c in a.obs.columns:
            v = a.obs[c]
            obs[c] = v.astype(str).values if str(v.dtype) == "category" else v.values
    obs["gsm"] = gsm
    obs["sample_label"] = label
    obs["depth_retained"] = depth
    obs["n_genes_retained"] = ngene

    # ---- panel counts for every cell ---------------------------------------
    present = [g for g in PANEL if g in gene_pos]
    cols = [gene_pos[g] for g in present]
    sub = counts[:, cols].toarray().astype(np.int32)
    for j, g in enumerate(present):
        obs[f"g_{g}"] = sub[:, j]
    missing = [g for g in PANEL if g not in gene_pos]

    # ---- §6 marker QC, computed independently of the deposited label -------
    qc_cols = [gene_pos[g] for g in MAST_QC if g in gene_pos]
    qc = counts[:, qc_cols].toarray() > 0
    obs["mast_markers_n"] = qc.sum(axis=1).astype(np.int8)
    obs["mast_qc_pass"] = obs["mast_markers_n"] >= 2

    labelled = obs.get("Cell type", pd.Series(["?"] * len(obs))).astype(str).values == "Mast"
    obs["mast_labelled"] = labelled

    keep = labelled | obs["mast_qc_pass"].values
    obs.to_parquet(obs_out)

    mast = counts[keep, :]
    sp.save_npz(mast_out, mast.tocsr())
    np.savez_compressed(
        OUT / f"{gsm}_mastmeta.npz",
        barcodes=np.asarray(obs.index[keep], dtype=object),
        var=var.astype(object),
        missing_panel=np.asarray(missing, dtype=object),
    )

    n_mast = int(keep.sum())
    del a, counts, sub
    for p in (h5, gz):
        try:
            p.unlink()
        except FileNotFoundError:
            pass
    return (f"{label}: {len(obs)} cells, {int(labelled.sum())} labelled Mast, "
            f"{int(obs['mast_qc_pass'].sum())} QC-pass, {n_mast} stored")


def main() -> int:
    mf = manifest()
    mf.to_csv(OUT / "manifest.csv", index=False)
    print(f"{len(mf)} human samples")
    workers = int(os.environ.get("WORKERS", "3"))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(process, r.gsm, r.fname, r.label) for r in mf.itertuples()]
        for f in futs:
            try:
                print(" ", f.result(), flush=True)
            except Exception as e:  # keep going; report at the end
                print("  ERROR:", e, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
