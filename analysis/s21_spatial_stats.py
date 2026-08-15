"""Spatial statistics — GSE197023 Visium (in-situ, no dissociation).

Subject: MAST CELLS (CLAUDE.md §2). Fibroblast and keratinocyte content are the
mandated control regressions (§6): if TNFRSF9 tracks those as strongly as it
tracks mast content, the association is not mast-specific.

Two questions, in the contract's order (§3):
  A. Does whole-spot TNFRSF9 differ between healthy and AD skin?
  B. Within the tissue, does TNFRSF9 concentrate where mast cells are?

Both models carry log(spot depth) as an offset and cluster standard errors on
donor. The depth control is not cosmetic here: median spot depth is ~4x higher
in the lesional arm than in healthy skin, so an unadjusted comparison would be a
depth comparison (§6).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import MAST_QC, TARGET, TARGET_LIGAND  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE197023"
TAB = ROOT / "results" / "tables"
TAB.mkdir(parents=True, exist_ok=True)

SUBJECT = G.require_mast_subject("Mast")

# §6 mandated control populations for the co-localisation model.
CONTENT_SETS = {
    "Mast": MAST_QC,                                  # subject
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN"],       # control
    "Keratinocytes": ["KRT14", "KRT5", "KRT10"],      # control
}


def fit(d: pd.DataFrame, y: np.ndarray, X: np.ndarray, off: np.ndarray,
        groups: np.ndarray, names: list[str]) -> pd.DataFrame:
    res = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit(
        cov_type="cluster", cov_kwds={"groups": groups})
    rows = []
    for i, n in enumerate(names):
        b, se = float(res.params[i]), float(res.bse[i])
        rows.append(dict(term=n, log2FC=b / np.log(2), se_log2=se / np.log(2),
                         ci_lo=(b - 1.96 * se) / np.log(2),
                         ci_hi=(b + 1.96 * se) / np.log(2),
                         p=float(res.pvalues[i])))
    return pd.DataFrame(rows)


def main() -> int:
    s = pd.read_parquet(PROC / "spots.parquet")
    s = s[s.total_counts >= 200].copy()          # minimal spot QC
    s["arm_pooled"] = np.where(s.arm == "Healthy", "Healthy", "AD")

    print("§6 median spot depth by arm:")
    print(s.groupby("arm", observed=True)["total_counts"]
          .agg(["size", "median"]).to_string())

    # ---------------- A. group differences, contract order (§3) ------------
    order = G.comparison_order([G.PRIMARY, *G.SECONDARY, G.REFINEMENT])
    specs = {
        G.PRIMARY: ("arm_pooled", ("Healthy", "AD")),
        "Healthy_vs_AD_NL": ("arm", ("Healthy", "AD_NL")),
        "Healthy_vs_AD_LS": ("arm", ("Healthy", "AD_LS")),
        G.REFINEMENT: ("arm", ("AD_NL", "AD_LS")),
    }
    genes = [TARGET, TARGET_LIGAND] + MAST_QC + ["COL1A1", "KRT14"]
    rows = []
    for comp in order:
        key, (ref, alt) = specs[comp]
        d = s[s[key].isin([ref, alt])]
        x = (d[key] == alt).astype(float).values
        X = sm.add_constant(np.column_stack([x]), has_constant="add")
        off = np.log(d.total_counts.astype(float).values)
        G.check_model_terms(["intercept", "arm", "log10_depth(offset=log spot depth)"],
                            grouping="donor")
        for g in genes:
            col = f"g_{g}"
            if col not in d:
                continue
            r = fit(d, d[col].astype(float).values, X, off, d.donor.values,
                    ["intercept", "arm"]).iloc[1].to_dict()
            r.update(gene=g, comparison=comp,
                     n_spots_ref=int((~(d[key] == alt)).sum()),
                     n_spots_alt=int((d[key] == alt).sum()),
                     donors_ref=int(d.loc[d[key] == ref, "donor"].nunique()),
                     donors_alt=int(d.loc[d[key] == alt, "donor"].nunique()),
                     cp10k_ref=1e4 * d.loc[d[key] == ref, col].sum() / d.loc[d[key] == ref, "total_counts"].sum(),
                     cp10k_alt=1e4 * d.loc[d[key] == alt, col].sum() / d.loc[d[key] == alt, "total_counts"].sum())
            rows.append(r)
    A = pd.DataFrame(rows)
    A.to_csv(TAB / "spatial_group_glm.csv", index=False)
    print("\nA. TNFRSF9 by arm (spot-level, depth-offset, donor-clustered):")
    print(A[A.gene.isin([TARGET] + MAST_QC)][
        ["comparison", "gene", "log2FC", "ci_lo", "ci_hi", "p",
         "cp10k_ref", "cp10k_alt", "donors_ref", "donors_alt"]].to_string(index=False))

    # ---------------- B. in-situ co-localisation ---------------------------
    # Does TNFRSF9 concentrate in mast-rich spots, over and above depth and arm?
    rows = []
    for name, gset in CONTENT_SETS.items():
        cols = [f"g_{g}" for g in gset if f"g_{g}" in s.columns]
        content = s[cols].sum(axis=1)
        cpm = 1e4 * content / s.total_counts
        z = np.log1p(cpm)
        z = (z - z.mean()) / z.std()            # per-SD effect, comparable across sets
        for arm in ["Healthy", "AD_NL", "AD_LS"]:
            d = s[s.arm == arm]
            zz = z[s.arm == arm]
            X = sm.add_constant(np.column_stack([zz.values]), has_constant="add")
            off = np.log(d.total_counts.astype(float).values)
            G.check_model_terms(["intercept", f"{name}_content",
                                 "log10_depth(offset=log spot depth)"],
                                grouping="donor")
            r = fit(d, d[f"g_{TARGET}"].astype(float).values, X, off,
                    d.donor.values, ["intercept", "content"]).iloc[1].to_dict()
            r.update(content=name, arm=arm, n_spots=len(d),
                     donors=int(d.donor.nunique()))
            rows.append(r)
    B = pd.DataFrame(rows)
    B.to_csv(TAB / "spatial_colocalisation.csv", index=False)
    print("\nB. TNFRSF9 per SD of cell-type content (in situ, depth-offset):")
    print(B[["content", "arm", "log2FC", "ci_lo", "ci_hi", "p", "n_spots", "donors"]]
          .to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
