"""Muted manuscript palette and Matplotlib defaults (CLAUDE.md §8).

Figures: Nimbus Sans, muted palette, minimal on-plot text, legible at
manuscript scale (single column = 89 mm, double column = 183 mm).
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ---------------------------------------------------------------- palette ---
# Muted, colour-vision-safe. Reserved semantics:
#   HEALTHY  cool grey-blue      AD_NL  muted amber      AD_LS  muted red
MUTED = {
    "healthy":   "#6E8CA0",
    "ad_nl":     "#C89B4A",
    "ad_ls":     "#A65344",
    "ad_all":    "#B06A48",
    "mast":      "#7A5C7B",
    "tcell":     "#8FA36B",
    "macro":     "#7FA1A8",
    "kerat":     "#B8A08A",
    "fibro":     "#98867A",
    "endo":      "#8C8FB0",
    "other":     "#9E9E9E",
    "grey":      "#4D4D4D",
    "grey_mid":  "#8A8A8A",
    "grey_light": "#D6D6D6",
    "accent":    "#A65344",
    "null":      "#BFBFBF",
}

GROUP_COLORS = {
    "Healthy": MUTED["healthy"],
    "AD": MUTED["ad_all"],
    "AD non-lesional": MUTED["ad_nl"],
    "AD lesional": MUTED["ad_ls"],
    "AD_NL": MUTED["ad_nl"],
    "AD_LS": MUTED["ad_ls"],
}

CELLTYPE_COLORS = {
    "Mast": MUTED["mast"],
    "T/NK": MUTED["tcell"],
    "Macrophages": MUTED["macro"],
    "DC": MUTED["endo"],
    "Keratinocytes": MUTED["kerat"],
    "Fibroblasts": MUTED["fibro"],
    "VEC": MUTED["endo"],
    "Melanocytes": MUTED["other"],
}

# Ordered categorical sequence for generic use.
SEQUENCE = [
    MUTED["healthy"], MUTED["ad_ls"], MUTED["ad_nl"], MUTED["mast"],
    MUTED["tcell"], MUTED["macro"], MUTED["kerat"], MUTED["fibro"],
]


def _font() -> str:
    """Return Nimbus Sans if installed, else the closest metric-compatible fallback."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Nimbus Sans", "Nimbus Sans L", "Helvetica", "Arial",
                 "FreeSans", "DejaVu Sans"):
        if name in available:
            return name
    return "sans-serif"


FONT = _font()


def set_style() -> None:
    """Apply manuscript defaults. Call once at the top of every figure script."""
    plt.rcParams.update({
        "font.family": FONT,
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "lines.linewidth": 0.9,
        "figure.dpi": 200,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


MM = 1 / 25.4  # millimetres -> inches, for manuscript-width figures
SINGLE_COL = 89 * MM
DOUBLE_COL = 183 * MM


def panel_label(ax, letter: str, text: str = "", pad: float = 4,
                fontsize: float = 6.6) -> None:
    """Panel title: BOLD letter, PLAIN description — the single definition.

    Every figure in this project must call this rather than passing
    ``fontweight="bold"`` to ``set_title``, which bolds the description as well
    and produces the inconsistent labelling this helper exists to prevent.
    ``analysis/guardrails.py`` enforces it.
    """
    ax.set_title(rf"$\bf{{{letter}}}$   {text}" if text else rf"$\bf{{{letter}}}$",
                 loc="left", fontsize=fontsize, fontweight="normal", pad=pad)
