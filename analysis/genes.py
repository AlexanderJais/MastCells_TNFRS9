"""Gene panel definitions shared across the project.

Kept in one place so that the mast-cell marker QC (CLAUDE.md §6) and the
control-population definitions (§2) are identical in every script.
"""
from __future__ import annotations

# --- the question -----------------------------------------------------------
TARGET = "TNFRSF9"          # CD137 / 4-1BB
TARGET_LIGAND = "TNFSF9"    # 4-1BBL — receptor engagement needs a ligand source

# TNFRSF co-stimulatory relatives, for specificity of any TNFRSF9 result.
TNFRSF_FAMILY = ["TNFRSF9", "TNFSF9", "TNFRSF4", "TNFRSF18", "TNFRSF8",
                 "TNFRSF1A", "TNFRSF1B", "TNFRSF14", "CD40", "TNFRSF12A"]

# --- mast cells -------------------------------------------------------------
# §6 marker QC set: a cell must show >=2 of these to be counted as a mast cell.
MAST_QC = ["TPSAB1", "TPSB2", "CPA3"]

# Wider mast identity panel (lineage + granule + receptors).
MAST_IDENTITY = ["TPSAB1", "TPSB2", "TPSD1", "CPA3", "CMA1", "CTSG", "MS4A2",
                 "KIT", "HDC", "HPGDS", "GATA2", "SLC18A2", "VWA5A", "LTC4S",
                 "IL1RL1", "SIGLEC6", "RGS13", "CD9", "ADCYAP1", "SRGN",
                 "FCER1A", "FCER1G", "ENPP3", "CD34", "ITGB7"]

# Mast subtype axis (MC-TC = tryptase+chymase, MC-T = tryptase only).
MAST_SUBTYPE = ["CMA1", "CTSG", "CDH26", "TPSAB1", "TPSB2"]

# --- control populations (§2: controls only, never findings) -----------------
CONTROL_MARKERS = {
    "T/NK":          ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "IL7R", "NKG7", "CCL5"],
    "Macrophages":   ["CD68", "CD163", "MRC1", "AIF1", "LYZ", "C1QA", "C1QB"],
    "DC":            ["CD1C", "ITGAX", "LAMP3", "CLEC9A", "HLA-DRA"],
    "Keratinocytes": ["KRT14", "KRT5", "KRT1", "KRT10", "KRT6A", "KRT16",
                      "S100A8", "S100A9", "FLG", "LOR", "KRT2"],
    "Fibroblasts":   ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "APOD"],
    "VEC":           ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "LEC":           ["LYVE1", "PROX1", "CCL21"],
    "Melanocytes":   ["PMEL", "MLANA", "TYRP1", "DCT"],
}

# §6 requires the same test in fibroblasts + keratinocytes; T cells and
# macrophages are carried as immune-compartment controls.
REQUIRED_CONTROL_POPULATIONS = ["Fibroblasts", "Keratinocytes", "T/NK", "Macrophages"]

# --- §6 dissociation-stress burden -----------------------------------------
# Heat-shock and immediate-early genes: a stressed arm reads artificially low.
HSP_GENES = ["HSPA1A", "HSPA1B", "HSPA6", "HSPB1", "HSPH1", "HSP90AA1",
             "HSP90AB1", "DNAJA1", "DNAJB1", "DNAJB6", "HSPE1", "HSPD1"]
IEG_GENES = ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "IER2",
             "IER3", "ZFP36", "KLF6", "NR4A1", "DUSP1", "PPP1R15A", "SOCS3"]
STRESS_GENES = HSP_GENES + IEG_GENES

# --- mast-cell activation / type-2 context ---------------------------------
MAST_ACTIVATION = ["TNF", "IL6", "CXCL8", "VEGFA", "CCL2", "CCL3", "NR4A1",
                   "NR4A2", "NR4A3", "EGR2", "EGR3", "IL13", "IL4", "IL5",
                   "IL31", "OSM", "AREG", "TGFB1", "PTGS2", "CSF2"]

TYPE2_CONTEXT = ["IL13", "IL4", "IL5", "IL31", "TSLP", "IL33", "CCL17",
                 "CCL18", "CCL13", "CCL22", "IL9R", "IL17A", "IL22", "IFNG"]

HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "RPL13A", "RPS18", "TMSB4X", "UBC"]


def panel() -> list[str]:
    """Full deduplicated extraction panel, order-stable."""
    out: list[str] = []
    seen: set[str] = set()
    groups = [TNFRSF_FAMILY, MAST_IDENTITY, MAST_SUBTYPE, MAST_ACTIVATION,
              TYPE2_CONTEXT, STRESS_GENES, HOUSEKEEPING]
    groups += list(CONTROL_MARKERS.values())
    for grp in groups:
        for g in grp:
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out


PANEL = panel()
