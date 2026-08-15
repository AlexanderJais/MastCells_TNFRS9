# Datasets screened and used

Selection rule (CLAUDE.md §5): dedicated AD research first; healthy and AD arms
from the same study or the same laboratory pipeline; single-cell before
single-nucleus; spatial platforms included; bulk included.

## Used

| Cohort | Accession | Platform | Composition | Role |
|---|---|---|---|---|
| Discovery scRNA-seq | **GSE204762** | 10x 3′ scRNA-seq, whole-skin biopsy | 11 AD donors (paired non-lesional + lesional), 6 healthy donors, 280,518 cells passing filters | Primary. Raw counts, author cell-type labels, one pipeline for both arms |
| Spatial | **GSE197023** | 10x Visium | 7 AD lesional, 6 AD non-lesional, 6 healthy sections; 13,742 spots | In-situ arm: no dissociation |
| Bulk | **GSE121212** | Illumina RNA-seq, whole skin | 38 healthy, 27 AD lesional, 27 AD non-lesional; 28 psoriasis lesional + 27 non-lesional as disease control | Tissue-level, large n, disease-specificity control |
| Replication 1 | **GSE222840** + **GSE173205** | 10x 5′ scRNA-seq, whole-skin biopsy | 5 AD lesional (GSE222840) + 4 healthy (GSE173205) | Same laboratory, same protocol; healthy arm from the companion series — a cross-series comparator, disclosed per §5 |
| Replication 2 | **GSE153760** | 10x 3′ v3 scRNA-seq | 4 AD + 2 healthy **biopsies** | Single study, single chemistry |

## Screened and excluded, with reason

| Accession | Reason for exclusion |
|---|---|
| **GSE147424** (He et al. 2020, 5 AD / 7 healthy, cryopreserved biopsies) | **No mast cells recoverable.** TPSAB1, TPSB2, CPA3, MS4A2, CMA1, HDC and CTSG are absent from the deposited gene list of *all 17* samples, and no cell in any sample tested is triple-positive for the surrogate markers KIT/HPGDS/GATA2 (0/4,147, 0/2,578, 0/4,760). The deposit is also library-size-normalised rather than raw counts, so the §6 depth control cannot be applied. Reported here as a negative result about cryopreserved-biopsy dissociation, not silently dropped. |
| GSE153760 suction-blister samples (AD1–4, HC1–5) | The authors state that suction blistering does not recover mast cells; only the biopsy samples from this series are used. |
| GSE206391 (Schäbitz et al., Visium, 90 sections) | Retained as a secondary spatial resource; the primary spatial arm (GSE197023) already contains healthy, AD non-lesional and AD lesional skin from one study, which GSE206391's AD arm does not pair with an in-study healthy control. |
| GSE157194 (Möbus et al., bulk, 166 samples) | **No healthy arm** — AD lesional/non-lesional only. Usable for the §3 refinement comparison but not for the primary Healthy-vs-AD question, so it is not carried as a main cohort. |
| GSE204762 mouse samples, scleroderma samples | Out of scope: the question is human AD. |
| GSM5907087 (AD_4_NL, GSE197023) | Absent from the GEO supplementary archive although listed in the series matrix — a gap in the deposit. The spatial non-lesional arm therefore has 6, not 7, sections. |

## Data provenance notes affecting interpretation

* **GSE204762** deposits mitochondrial genes stripped from the count matrix, so
  `mt_frac` is carried from the authors' own QC rather than recomputed. Library
  size used throughout is the *retained* library size (sum over deposited genes),
  which is what governs detection in this matrix.
* **GSE204762** samples were deposited with two different gene axes (27,364 and
  29,991 genes); all cross-sample analyses are remapped onto the 27,030-gene
  intersection.
* **GSE121212** AD lesional comprises 21 samples labelled `lesional` and 6
  labelled `chronic_lesion`; both are treated as lesional AD.
* Chemistry is mixed within the discovery cohort (one 3′ v2 sample per arm,
  the rest v3) and between replication cohorts (5′ in REP1, 3′ v3 in REP2).

## In-vitro stimulation cohorts (added to test the mechanism directly)

| Cohort | Accession | Cells | Design | Role |
|---|---|---|---|---|
| Skin mast cells + AD milieu | **GSE196862** | Primary human **skin**-derived mast cells (HSMCs) | 24 libraries: unstimulated, IgE/Ag cross-linking, IL-33, TSLP, IL-25 and combinations; Cufflinks FPKM | Tests whether TNFRSF9 is inducible in the correct cell type and tissue by the AD-defining stimuli |
| Mast cells + activated T cells | **GSE235240** | Primary human mast cells, FACS-sorted after co-culture | 4 donors, paired: control, resting CD4⁺ Teff, anti-CD3/CD28-activated CD4⁺ Teff, IL-33, IgE/Ag | Independent confirmation with an independent stimulus and a paired donor design |

These are in-vitro systems, so they establish sufficiency of a stimulus, not the
in-vivo contribution. They are reported in §2B of the report as *measured* results
in primary cells, distinct from the *modelled* tissue estimates.
