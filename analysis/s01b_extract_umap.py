"""Second pass over GSE204762: retrieve the authors' integrated UMAP embedding.

The deposited per-sample h5ad files are slices of one Harmony-integrated object,
so their obsm['X_umap'] coordinates live in a common embedding and can be
concatenated into an atlas-wide map. Only the embedding and the cell index are
read here; the count matrices were already reduced in s01.
"""
from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE204762"
OUT = ROOT / "data" / "processed" / "GSE204762"
SAMPLE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{stem}nnn/{gsm}/suppl/{fname}"


def one(gsm: str, fname: str, label: str) -> str:
    dest = OUT / f"{gsm}_umap.parquet"
    if dest.exists():
        return f"{label}: cached"
    gz = RAW / fname
    h5 = RAW / fname[:-3]
    if not h5.exists():
        if not (gz.exists() and gz.stat().st_size > 1_000_000):
            url = SAMPLE_URL.format(stem=gsm[:-3], gsm=gsm, fname=fname)
            for attempt in range(4):
                rc = subprocess.run(["curl", "-sS", "-L", "--fail", "--max-time",
                                     "1800", "-o", str(gz), url]).returncode
                if rc == 0 and gz.stat().st_size > 1_000_000:
                    break
                import time
                time.sleep(2 ** (attempt + 1))
            else:
                raise RuntimeError(f"download failed {gsm}")
        with gzip.open(gz, "rb") as fi, open(h5, "wb") as fo:
            shutil.copyfileobj(fi, fo, length=1 << 24)

    import h5py
    with h5py.File(h5, "r") as f:
        um = np.asarray(f["/obsm/X_umap"][:], dtype=np.float32)
        idx = f["/obs/_index"][:]
        idx = np.array([x.decode() if isinstance(x, bytes) else x for x in idx])
    df = pd.DataFrame(um[:, :2], columns=["umap1", "umap2"], index=idx)
    df["gsm"] = gsm
    df.to_parquet(dest)
    for p in (h5, gz):
        try:
            p.unlink()
        except FileNotFoundError:
            pass
    return f"{label}: {len(df)} cells"


def main() -> int:
    mf = pd.read_csv(OUT / "manifest.csv")
    with ThreadPoolExecutor(max_workers=3) as ex:
        futs = [ex.submit(one, r.gsm, r.fname, r.label) for r in mf.itertuples()]
        for f in futs:
            try:
                print(" ", f.result(), flush=True)
            except Exception as e:
                print("  ERROR:", e, flush=True)
    parts = [pd.read_parquet(p) for p in sorted(OUT.glob("*_umap.parquet"))]
    allu = pd.concat(parts)
    allu.to_parquet(OUT / "umap_all.parquet")
    print(f"\nconcatenated UMAP: {len(allu)} cells")
    print(allu[["umap1", "umap2"]].describe().round(2).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
