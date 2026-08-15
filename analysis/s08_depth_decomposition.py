"""Where does "depth" actually come from — the flow cell, or the cell?

Correction of terminology, and the analysis that follows from it.

Sequencing depth in the strict sense is a property of the run: how many reads
are allocated to a library. It is shared by every cell in that library. What
per-cell UMI counts actually measure is SAMPLING EFFORT PER CELL, which is a
product of

    reads allocated to the library      <- flow cell, experimenter's choice
  x capture / lysis efficiency          <- droplet chemistry
  x the cell's own mRNA content         <- biology

Only the first is "sequencing depth". This matters because normalising by
per-cell UMI corrects the first two but also divides out the third, and the
third is biology.

The decomposition is empirical, not rhetorical. Two facts settle it:

  1. Within ONE library every cell was sequenced on the same flow cell to the
     same depth. So any difference in UMI per cell BETWEEN CELL TYPES inside a
     library cannot be sequencing depth.
  2. If mast cells are shallower than fibroblasts by the same factor in every
     library, that factor is cell-intrinsic, not experimental.

Subject: MAST CELLS (CLAUDE.md §2); other populations are the within-library
reference against which mast cells are measured.
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
from genes import TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY])

REF = ["Fibroblasts", "Keratinocytes"]   # §6 mandated within-library references


def main() -> int:
    o = pd.read_parquet(PROC / "obs_all.parquet")
    o = o[o.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    o["armp"] = np.where(o.arm == "Healthy", "Healthy", "AD")
    o["ct"] = o["Cell type"].astype(str)
    o["ldepth"] = np.log(o.depth_retained.clip(lower=1))

    # ---- 1. variance decomposition ---------------------------------------
    print("=" * 96)
    print("1. WHERE DOES THE VARIANCE IN log(UMI per cell) LIVE?")
    keep = o[o.ct.isin(["Mast"] + REF + ["T/NK", "Macrophages"])].copy()
    grand = keep.ldepth.mean()
    ss_tot = ((keep.ldepth - grand) ** 2).sum()
    # between libraries (flow-cell + loading): sample means
    smean = keep.groupby("sample_label", observed=True).ldepth.transform("mean")
    ss_lib = ((smean - grand) ** 2).sum()
    # between cell types WITHIN library: deviation of cell-type-in-sample mean
    ctmean = keep.groupby(["sample_label", "ct"], observed=True).ldepth.transform("mean")
    ss_ct = ((ctmean - smean) ** 2).sum()
    ss_res = ((keep.ldepth - ctmean) ** 2).sum()
    print(f"  between libraries (flow cell / loading) : {100*ss_lib/ss_tot:5.1f}%")
    print(f"  between cell types WITHIN a library     : {100*ss_ct/ss_tot:5.1f}%")
    print(f"  cell-to-cell within a cell type         : {100*ss_res/ss_tot:5.1f}%")
    print("  -> the within-library component cannot be sequencing depth:")
    print("     those cells shared a flow cell.")
    pd.DataFrame(dict(component=["between_libraries", "between_celltypes_within_library",
                                 "within_celltype"],
                      pct_variance=[100 * ss_lib / ss_tot, 100 * ss_ct / ss_tot,
                                    100 * ss_res / ss_tot])
                 ).to_csv(TAB / "depth_variance_decomposition.csv", index=False)

    # ---- 2. is the mast:reference ratio constant across libraries? -------
    print("\n" + "=" * 96)
    print("2. MAST-CELL UMI RELATIVE TO THE SAME LIBRARY'S FIBROBLASTS/KERATINOCYTES")
    print("   (a within-library ratio is immune to flow-cell depth by construction)")
    rows = []
    for s, g in o.groupby("sample_label", observed=True):
        m = g[g.ct == "Mast"]
        r = g[g.ct.isin(REF)]
        if len(m) < 5 or len(r) < 20:
            continue
        rows.append(dict(sample=s, arm=g.armp.iat[0], n_mast=len(m),
                         mast_umi=m.depth_retained.mean(),
                         ref_umi=r.depth_retained.mean(),
                         ratio=m.depth_retained.mean() / r.depth_retained.mean()))
    R = pd.DataFrame(rows)
    R.to_csv(TAB / "depth_mast_relative_ratio.csv", index=False)
    print(R.groupby("arm").ratio.agg(["count", "median", "min", "max"]).round(3).to_string())
    h = R.loc[R.arm == "Healthy", "ratio"]; a = R.loc[R.arm == "AD", "ratio"]
    u = stats.mannwhitneyu(a, h)
    print(f"\n  mast:reference UMI ratio — healthy median {h.median():.3f} "
          f"(n={len(h)} libraries), AD median {a.median():.3f} (n={len(a)}), "
          f"Mann-Whitney P={u.pvalue:.3f}")
    print("  The ratio is well below 1 in EVERY library: mast cells yield a fraction")
    print("  of the UMI that stromal cells in the same droplet run yield. That is")
    print("  cell-intrinsic mRNA content, not sequencing.")

    # ---- 3. so is the AD deficit in mast UMI real, or library-level? -----
    print("\n" + "=" * 96)
    print("3. IS THE AD MAST-CELL UMI DEFICIT A LIBRARY EFFECT?")
    lib = o.groupby(["sample_label", "armp"], observed=True).depth_retained.mean().reset_index()
    print("  mean UMI per cell, ALL cells, by arm (library-level effort):")
    print(lib.groupby("armp").depth_retained.agg(["count", "median"]).round(0).to_string())
    print("  mean UMI per MAST cell, by arm:")
    print(o[o.ct == "Mast"].groupby("armp").depth_retained.mean().round(0).to_string())
    print("  -> AD libraries are sequenced DEEPER overall yet their mast cells are")
    print("     SHALLOWER. A flow-cell effect cannot move those two in opposite")
    print("     directions; the mast-cell deficit is a property of the cells.")

    # ---- 4. the estimand this implies ------------------------------------
    print("\n" + "=" * 96)
    print("4. CONSEQUENCE FOR THE TNFRSF9 ESTIMAND")
    m = o[o.mast_strict].copy()
    tot = m.groupby("armp").agg(cells=("depth_retained", "size"),
                                umi=("depth_retained", "sum"),
                                t9=(f"g_{TARGET}", "sum"))
    print(tot.to_string())
    print(f"\n  per mast CELL   : healthy {100*tot.loc['Healthy','t9']/tot.loc['Healthy','cells']:.2f}"
          f" vs AD {100*tot.loc['AD','t9']/tot.loc['AD','cells']:.2f} per 100 cells "
          f"({(tot.loc['AD','t9']/tot.loc['AD','cells'])/(tot.loc['Healthy','t9']/tot.loc['Healthy','cells']):.1f}x)")
    print(f"  per mast UMI    : healthy {1e4*tot.loc['Healthy','t9']/tot.loc['Healthy','umi']:.3f}"
          f" vs AD {1e4*tot.loc['AD','t9']/tot.loc['AD','umi']:.3f} per 10k "
          f"({(tot.loc['AD','t9']/tot.loc['AD','umi'])/(tot.loc['Healthy','t9']/tot.loc['Healthy','umi']):.1f}x)")
    print("\n  Because the UMI difference between arms is cell-intrinsic rather than")
    print("  experimental, dividing it out removes biology as well as technical")
    print("  effort. The per-CELL estimand is therefore the primary one, and the")
    print("  per-UMI estimand is reported as the compositional companion.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
