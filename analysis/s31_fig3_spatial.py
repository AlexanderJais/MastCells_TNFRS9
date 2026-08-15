"""Figure 3 — TNFRSF9 in situ (GSE197023, 10x Visium).

Subject: MAST CELLS (CLAUDE.md §2). Visium never dissociates the tissue, so this
arm tests whether the single-cell picture survives without enzymatic digestion.

Honesty note carried onto the figure itself: median UMI per spot differs ~4x
between arms (AD lesional 4,005 vs healthy 1,157 UMI), so every section is
annotated with its own median depth and all statistics carry a log-depth offset
(§6). A raw side-by-side image comparison would be a depth comparison.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402

import guardrails as G  # noqa: E402
from genes import MAST_QC, TARGET  # noqa: E402
from palette import DOUBLE_COL, MUTED, set_style  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE197023"
PROC = ROOT / "data" / "processed" / "GSE197023"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SUBJECT = G.require_mast_subject("Mast")
COMPARISONS = G.comparison_order([G.PRIMARY, *G.SECONDARY, G.REFINEMENT])
ARM_COL = {"Healthy": MUTED["healthy"], "AD_NL": MUTED["ad_nl"], "AD_LS": MUTED["ad_ls"]}
ARM_LAB = {"Healthy": "Healthy", "AD_NL": "AD non-lesional", "AD_LS": "AD lesional"}

# two sections per arm, chosen as the two deepest per arm so that the visual
# comparison is as depth-fair as the data allow
SECTIONS = ["HE_1", "HE_2", "AD_1_NL", "AD_7_NL", "AD_6_LS", "AD_7_LS"]


def section_image(name: str):
    d = RAW / name / "spatial"
    img = mpimg.imread(d / "tissue_lowres_image.png")
    sf = json.load(open(d / "scalefactors_json.json"))["tissue_lowres_scalef"]
    return img, sf


def spot_pitch(d: pd.DataFrame, sf: float) -> float:
    """Centre-to-centre spot spacing in display (lowres-image) units.

    Visium spot pitch is constant within a capture area but the fullres pixel
    scale differs between samples (109-227 px here), so a fixed marker size
    renders one section as a clean lattice and another as a merged blob.
    """
    from scipy.spatial import cKDTree
    xy = d[["x", "y"]].dropna().values.astype(float) * sf
    if len(xy) < 5:
        return 1.0
    dist, _ = cKDTree(xy).query(xy, k=2)
    return float(np.median(dist[:, 1]))


def panel(ax, s: pd.DataFrame, name: str, mode: str):
    """Draw one section. Marker sizes are set later, once the axes geometry is
    final, so that every section is drawn at true spot scale."""
    img, sf = section_image(name)
    ax.imshow(img)
    d = s[s["sample"] == name]
    x, y = d.x.values * sf, d.y.values * sf
    # crop to the tissue itself: the capture area is mostly empty slide
    cx, cy = (x.max() + x.min()) / 2, (y.max() + y.min()) / 2
    w, h = (x.max() - x.min()) * 1.10, (y.max() - y.min()) * 1.10
    side = max(w, h)                      # square crop -> no letterboxing
    ax.set_xlim(cx - side / 2, cx + side / 2)
    ax.set_ylim(cy + side / 2, cy - side / 2)
    pitch = spot_pitch(d, sf)
    if mode == "mast":
        tryp = d[[f"g_{g}" for g in MAST_QC]].sum(axis=1)
        v = np.log1p(1e4 * tryp / d.total_counts)
        sc = ax.scatter(x, y, c=v, s=1, lw=0, cmap="BuPu", vmin=0, vmax=4.2,
                        rasterized=True)
        return sc, pitch, 0.92
    sc = ax.scatter(x, y, c="#D2D2D2", s=1, lw=0, rasterized=True)
    pos = d[d[f"g_{TARGET}"] > 0]
    arm = d.arm.iat[0]
    hit = ax.scatter(pos.x.values * sf, pos.y.values * sf, s=1, lw=.4,
                     facecolor=ARM_COL[arm], edgecolor="#141414", zorder=5)
    return (sc, hit), pitch, 0.92


def size_for(ax, fig, pitch: float, frac: float) -> float:
    """Marker area in points^2 that makes a spot `frac` of the true spot pitch."""
    p0 = ax.transData.transform((0.0, 0.0))
    p1 = ax.transData.transform((pitch, 0.0))
    d_pt = abs(p1[0] - p0[0]) * 72.0 / fig.dpi
    return max((d_pt * frac) ** 2, 0.5)


def main() -> int:
    set_style()
    s = pd.read_parquet(PROC / "spots.parquet")
    s = s[s.total_counts >= 200]

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.72))
    gs = fig.add_gridspec(3, 6, height_ratios=[1, 1, 1.05], hspace=0.14, wspace=0.14)
    gsb = gs[2, :].subgridspec(1, 3, wspace=0.62)

    # ---- rows 1-2: in-situ maps ------------------------------------------
    pending = []           # (axes, scatter handles, pitch, size fraction)
    for j, name in enumerate(SECTIONS):
        d = s[s["sample"] == name]
        arm = d.arm.iat[0]
        ax = fig.add_subplot(gs[0, j])
        sc, pitch, frac = panel(ax, s, name, "mast")
        pending.append((ax, [sc], pitch, frac))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(f"{ARM_LAB[arm]}\n{name}", fontsize=5.6, pad=2,
                     color=ARM_COL[arm], fontweight="bold")
        if j == 0:
            ax.set_ylabel("tryptase/CPA3", fontsize=6)

        ax2 = fig.add_subplot(gs[1, j])
        (bg, hit), pitch2, frac2 = panel(ax2, s, name, "t9")
        pending.append((ax2, [bg], pitch2, frac2))
        pending.append((ax2, [hit], pitch2, 1.15))  # TNFRSF9+ spots: ~1 spot wide
        ax2.set_xticks([]); ax2.set_yticks([])
        for sp in ax2.spines.values():
            sp.set_visible(False)
        k = int((d[f"g_{TARGET}"] > 0).sum())
        ax2.set_title(f"{k} TNFRSF9$^+$ spots\nmedian {d.total_counts.median():.0f} UMI",
                      fontsize=5.2, pad=2)
        if j == 0:
            ax2.set_ylabel("TNFRSF9", fontsize=6)
    # geometry is final only after a draw; size every marker to its own grid
    fig.canvas.draw()
    for ax_, handles, pitch, frac in pending:
        s_pts = size_for(ax_, fig, pitch, frac)
        for h in handles:
            h.set_sizes([s_pts])

    cax = fig.add_axes([0.915, 0.68, 0.008, 0.15])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("log(1+tryptase per 10k)", fontsize=5.0)
    cb.ax.tick_params(labelsize=4.6)

    # ---- row 3a: group comparison (depth-adjusted) -----------------------
    ax = fig.add_subplot(gsb[0, 0])
    A = pd.read_csv(TAB / "spatial_group_glm.csv")
    A = A[(A.gene == TARGET)].set_index("comparison")
    comps = [G.PRIMARY, "Healthy_vs_AD_NL", "Healthy_vs_AD_LS"]
    labs = ["Healthy\nvs AD", "Healthy\nvs AD NL", "Healthy\nvs AD LS"]
    y = np.arange(len(comps))
    for i, c in enumerate(comps):
        r = A.loc[c]
        ax.plot([r.log2FC - r.se_log2, r.log2FC + r.se_log2], [i, i],
                color="#6B6B6B", lw=1.1)   # +/- 1 SEM
        ax.plot(r.log2FC, i, "o", ms=4, color=MUTED["ad_all"], zorder=5)
    ax.axvline(0, color="#B0B0B0", lw=.6, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.6)
    ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change", fontsize=6)
    ax.set_title("g   TNFRSF9 per spot  (± SEM)", loc="left", fontweight="bold")

    # ---- row 3b: in-situ co-localisation ---------------------------------
    ax = fig.add_subplot(gsb[0, 1])
    B = pd.read_csv(TAB / "spatial_colocalisation.csv")
    order = ["Healthy", "AD_NL", "AD_LS"]
    w = .26
    for k, content in enumerate(["Mast", "Fibroblasts", "Keratinocytes"]):
        sub = B[B.content == content].set_index("arm").loc[order]
        xs = np.arange(len(order)) + (k - 1) * w
        col = MUTED["mast"] if content == "Mast" else ("#C4C4C4" if content == "Fibroblasts" else "#E0E0E0")
        ax.bar(xs, sub.log2FC, w, color=col, edgecolor="none",
               label=content + (" (subject)" if content == "Mast" else " (control)"))
        ax.errorbar(xs, sub.log2FC, yerr=sub.se_log2, fmt="none",
                    ecolor="#6B6B6B", elinewidth=.7, capsize=1.4)   # +/- 1 SEM
    ax.axhline(0, color="#B0B0B0", lw=.6)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(["Healthy", "AD NL", "AD LS"], fontsize=6)
    ax.set_ylabel("log$_2$ TNFRSF9 per doubling\nof content  (± SEM)", fontsize=5.8, labelpad=1)
    ax.legend(fontsize=5, loc="upper left")
    ax.set_title("h   in-situ co-localisation", loc="left", fontweight="bold")

    # ---- row 3c: per-section rate ----------------------------------------
    ax = fig.add_subplot(gsb[0, 2])
    per = s.groupby("sample").apply(lambda d: pd.Series({
        "arm": d.arm.iat[0],
        "rate": 1e4 * d[f"g_{TARGET}"].sum() / d.total_counts.sum(),
        "depth": d.total_counts.median()}))
    per = per.sort_values(["arm", "rate"])
    ax.bar(np.arange(len(per)), per.rate.astype(float),
           color=[ARM_COL[a] for a in per.arm], edgecolor="none", width=.72)
    ax.set_xticks(np.arange(len(per)))
    ax.set_xticklabels(per.index, rotation=90, fontsize=4.6)
    ax.set_ylabel("TNFRSF9 per 10k spot UMI", fontsize=6, labelpad=1)
    ax.set_title("i   per section", loc="left", fontweight="bold")

    fig.savefig(FIG / "fig3_spatial.pdf")
    fig.savefig(FIG / "fig3_spatial.png")
    print("wrote fig3_spatial.pdf / .png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
