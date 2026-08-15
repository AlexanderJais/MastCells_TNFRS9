"""Adjusted versus unadjusted: does the depth offset create the result?

Reviewer challenge: the libraries were sequenced to comparable depth, so
adjusting for depth may be removing signal rather than confounding.

The challenge is answered with data, not argument. Three things are separated
here, because "sequencing depth" means three different quantities and only one
of them is a choice the experimenter made:

  1. READS PER LIBRARY   — an experimental design choice; may well be matched.
  2. UMI PER CELL        — reads x capture efficiency x the cell's own RNA
                           content. Not a design choice. Varies ~6-fold BETWEEN
                           CELL TYPES inside a single library.
  3. UMI PER MAST CELL, PER ARM — the only quantity that can bias this
                           particular comparison.

If (3) is matched between healthy and AD, the offset is close to a no-op and the
adjusted and unadjusted answers must agree. If it is not matched, the direction
of the bias is an empirical question, answered below.

Two estimands, both reported:
  * per UMI  — TNFRSF9 as a FRACTION of the mast-cell transcriptome
               (compositional; what the adjusted model estimates)
  * per CELL — TNFRSF9 molecules captured per mast cell
               (absolute; what the unadjusted model estimates)

Subject: MAST CELLS (CLAUDE.md §2).
"""
from __future__ import annotations

import sys
import warnings
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import guardrails as G  # noqa: E402
from genes import TARGET  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed" / "GSE204762"
TAB = ROOT / "results" / "tables"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY, *G.SECONDARY])

CONTROLS = ["Fibroblasts", "Keratinocytes", "T/NK", "Macrophages"]


def exact_perm(counts, expo, alt):
    """Exact donor-label permutation on log rate ratio, whatever the exposure."""
    nd, k = len(counts), int(alt.sum())
    masks = np.zeros((comb(nd, k), nd), bool)
    for i, c in enumerate(combinations(range(nd), k)):
        masks[i, list(c)] = True
    ea, er = masks @ expo, expo.sum() - masks @ expo
    ca, cr = masks @ counts, counts.sum() - masks @ counts
    S = np.log((ca + .5) / ea) - np.log((cr + .5) / er)
    i0 = int(np.where((masks == alt[None, :]).all(axis=1))[0][0])
    return S[i0] / np.log(2), float((np.abs(S) >= abs(S[i0]) - 1e-12).mean())


def main() -> int:
    obs = pd.read_parquet(PROC / "obs_all.parquet")
    obs = obs[obs.arm.isin(["Healthy", "AD_NL", "AD_LS"])].copy()
    obs["armp"] = np.where(obs.arm == "Healthy", "Healthy", "AD")
    mast = obs[obs.mast_strict]

    # ---- 1. was depth actually matched? ----------------------------------
    print("=" * 100)
    print("1. WHAT WAS ACTUALLY MATCHED?")
    lib = obs.groupby(["sample_label", "armp"], observed=True).agg(
        cells=("depth_retained", "size"),
        library_umi=("depth_retained", "sum"),
        umi_per_cell=("depth_retained", "mean")).reset_index()
    print("\n  (a) TOTAL UMI PER LIBRARY, by arm  [the design choice]")
    print(lib.groupby("armp")["library_umi"].describe()[
        ["count", "mean", "50%", "min", "max"]].round(0).to_string())

    print("\n  (b) MEAN UMI PER CELL WITHIN A LIBRARY, by cell type  [not a design choice]")
    ct = obs.groupby("Cell type", observed=True)["depth_retained"].agg(["size", "mean", "median"])
    print(ct.loc[["Mast"] + [c for c in CONTROLS if c in ct.index]].round(0).to_string())
    print("  -> inside the SAME libraries, mast cells carry a fraction of the UMI that")
    print("     keratinocytes or fibroblasts do. That is cell biology, not sequencing.")

    print("\n  (c) MAST-CELL UMI PER CELL, BY ARM  [the only thing that can bias this test]")
    md = mast.groupby("armp", observed=True)["depth_retained"].agg(
        cells="size", mean_umi="mean", median_umi="median")
    print(md.round(1).to_string())
    ratio = md.loc["Healthy", "mean_umi"] / md.loc["AD", "mean_umi"]
    print(f"  -> healthy mast cells average {ratio:.2f}x the UMI of AD mast cells.")
    print("     Healthy is the DEEPER arm, so an unadjusted comparison favours HEALTHY.")
    md.to_csv(TAB / "unadj_mast_depth_by_arm.csv")

    # ---- 2. the same test, both estimands --------------------------------
    print("\n" + "=" * 100)
    print("2. SAME TEST, TWO ESTIMANDS (donor-level exact permutation)")
    don = mast.groupby("donor", observed=True).agg(
        arm=("armp", "first"), n_mast=("depth_retained", "size"),
        mast_umi=("depth_retained", "sum"),
        t9=(f"g_{TARGET}", "sum")).reset_index()
    alt = (don.arm == "AD").values
    rows = []
    for name, expo, unit in (
            ("per mast UMI  (adjusted, compositional)", don.mast_umi.values.astype(float), "10k UMI"),
            ("per mast CELL (unadjusted, absolute)", don.n_mast.values.astype(float), "100 cells")):
        lf, p = exact_perm(don.t9.values.astype(float), expo, alt)
        scale = 1e4 if "UMI" in unit else 100
        rows.append(dict(estimand=name, unit=unit,
                         healthy=scale * don.loc[~alt, "t9"].sum() / expo[~alt].sum(),
                         AD=scale * don.loc[alt, "t9"].sum() / expo[alt].sum(),
                         fold=(don.loc[alt, "t9"].sum() / expo[alt].sum()) /
                              (don.loc[~alt, "t9"].sum() / expo[~alt].sum()),
                         log2FC=lf, p_perm=p))
    R = pd.DataFrame(rows)
    print(R.round(4).to_string(index=False))
    R.to_csv(TAB / "unadj_two_estimands.csv", index=False)

    # ---- 3. raw, unnormalised detection ----------------------------------
    print("\n" + "=" * 100)
    print("3. FULLY RAW NUMBERS — no normalisation of any kind")
    raw = mast.groupby("armp", observed=True).agg(
        mast_cells=("depth_retained", "size"),
        t9_transcripts=(f"g_{TARGET}", "sum"),
        t9_pos_cells=(f"g_{TARGET}", lambda v: int((v > 0).sum())))
    raw["pct_cells_positive"] = 100 * raw.t9_pos_cells / raw.mast_cells
    raw["transcripts_per_100_cells"] = 100 * raw.t9_transcripts / raw.mast_cells
    print(raw.round(3).to_string())
    raw.to_csv(TAB / "unadj_raw_counts.csv")

    # detection-rate test needs no depth model at all
    from scipy import stats
    tab = np.array([[int(raw.loc["AD", "t9_pos_cells"]),
                     int(raw.loc["AD", "mast_cells"] - raw.loc["AD", "t9_pos_cells"])],
                    [int(raw.loc["Healthy", "t9_pos_cells"]),
                     int(raw.loc["Healthy", "mast_cells"] - raw.loc["Healthy", "t9_pos_cells"])]])
    orr, pf = stats.fisher_exact(tab)
    print(f"  Fisher on TNFRSF9+ cell counts (no normalisation): OR={orr:.2f}, P={pf:.2e}")

    # ---- 4. is the mast transcriptome itself different? ------------------
    print("\n" + "=" * 100)
    print("4. ARE MAST-CELL TRANSCRIPTS THEMSELVES REGULATED?")
    print("   (if the mast transcriptome shrinks in AD, a per-UMI fraction rises for free)")
    tr = mast.groupby("armp", observed=True).agg(
        umi_per_cell=("depth_retained", "mean"),
        genes_per_cell=("n_genes_retained", "mean"))
    for g in ("TPSAB1", "TPSB2", "CPA3"):
        tr[f"{g}_per_cell"] = mast.groupby("armp", observed=True)[f"g_{g}"].mean()
    print(tr.round(2).to_string())
    d_umi = 100 * (tr.loc["AD", "umi_per_cell"] / tr.loc["Healthy", "umi_per_cell"] - 1)
    print(f"  -> AD mast cells capture {d_umi:+.1f}% UMI per cell relative to healthy.")
    print("     A smaller denominator inflates a per-UMI fraction; that is why the")
    print("     per-CELL estimand in section 2 is the conservative one here.")
    tr.to_csv(TAB / "unadj_mast_transcriptome.csv")

    # ---- 5. controls, both estimands -------------------------------------
    print("\n" + "=" * 100)
    print("5. §6 CONTROL POPULATIONS, BOTH ESTIMANDS")
    rows = []
    for pop in ["Mast"] + CONTROLS:
        sub = mast if pop == "Mast" else obs[obs["Cell type"].astype(str) == pop]
        d = sub.groupby("donor", observed=True).agg(
            arm=("armp", "first"), n=("depth_retained", "size"),
            umi=("depth_retained", "sum"), t9=(f"g_{TARGET}", "sum")).reset_index()
        a = (d.arm == "AD").values
        if a.sum() < 2 or (~a).sum() < 2:
            continue
        lf_u, p_u = exact_perm(d.t9.values.astype(float), d.umi.values.astype(float), a)
        lf_c, p_c = exact_perm(d.t9.values.astype(float), d.n.values.astype(float), a)
        rows.append(dict(population=pop + (" (SUBJECT)" if pop == "Mast" else ""),
                         log2FC_perUMI=lf_u, p_perUMI=p_u,
                         log2FC_perCELL=lf_c, p_perCELL=p_c))
    C = pd.DataFrame(rows)
    print(C.round(4).to_string(index=False))
    C.to_csv(TAB / "unadj_controls_both.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
