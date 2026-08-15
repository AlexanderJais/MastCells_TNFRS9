# PROJECT CONTRACT — read before any analysis

This file is loaded automatically at the start of every session. It is binding. If an instruction
here conflicts with what seems locally convenient, **this file wins**. If the user overrides it
explicitly, update this file in the same commit. I had to set up this contract because there seems to be huge problems with Claude Code's rigor recently.

---

## 1. THE QUESTION

> **Is TNFRSF9 (CD137, 4-1BB) expressed in cutaneous MAST CELLS in atopic dermatitis, and does the expression of TNFRSF9 change
> between healthy skin and AD skin?**

Deliverable standard: Nature-portfolio manuscript quality.

## 2. THE SUBJECT IS MAST CELLS

**Every analysis must be about mast cells.** T cells, macrophages, keratinocytes and fibroblasts
appear in this project in exactly one role: **as control populations**, to demonstrate that a
mast-cell result is specific and not a batch or global effect.

**BANNED:** characterising T cells or macrophages as findings in their own right. Regressing tissue
signal onto T-cell content, describing T-cell activation states, or reporting "the signal comes from
T cells" as a result. If a control population is interesting, that is someone else's paper. OUR STUDY IS ABOUT MAST CELLS!


## 3. THE COMPARISON HIERARCHY — IN THIS ORDER

1. **PRIMARY — Healthy vs AD** (lesional and non-lesional pooled). This is the disease question.
2. **THEN — Healthy vs AD non-lesional**, and **Healthy vs AD lesional**, separately.
3. **THEN — AD lesional vs AD non-lesional** (paired within patient). This is a *refinement within
   disease*, not the headline.

**BANNED:** leading with lesional-vs-non-lesional because it happens to be the batch-clean
comparison. Methodological convenience must not choose the scientific framing. If the primary
comparison is confounded, **say so and still report it first**, then explain the confound.

## 4. ABUNDANCE AND EXPRESSION ARE ONE QUANTITY

Mast-cell-derived TNFRSF9 = **abundance × per-cell expression**. Never present these as separate
questions. Always report both factors and their product. In bulk data they are not separable at all;
in scRNA-seq the per-cell estimate is conditional on which mast cells survived dissociation.

## 5. DATA SOURCES — DEDICATED AD RESEARCH FIRST

Single-cell RNA-seq (scRNA-seq) datasets (cell-resolved transcriptomics; no single-nucleus datasets unless absolutely unavoidable)
Spatial transcriptomics datasets (Visium, CosMx, Xenium, MERSCOPE, Slide-seq, etc.). YES, VISIUM IS A SPATIAL PLATFORM! DON'T DISMISS VISIUM!! THERE ARE VISIUM DATASET WITH AD SKIN!
Bulk RNA-seq datasets

The final output should be suitable for a high-quality scientific manuscript submission (Nature portfolio style).

Healthy and AD arms must come from the **same study or the same laboratory pipeline**. Cross-study
healthy comparators are reportable only with the confound stated in the same sentence.

## 6. MANDATORY CONTROLS ON EVERY CLAIM

| Control | Why |
|---|---|
| Adjust for library size (`log10(depth)`) | TNFRSF9 detection is depth-dominated; unadjusted comparisons are depth comparisons |
| Cluster / aggregate by donor | cell-level tests pseudo-replicate |
| Report the same test in fibroblasts + keratinocytes | a shared shift means batch, not biology |
| Report per-group median sequencing depth | shows depth is not producing or hiding the effect |
| Marker-QC mast cells (≥2 of TPSAB1/TPSB2/CPA3) | deposited labels are unreliable |
| Dissociation-stress burden (HSP/IEG) per group | stressed arms read artificially low |

## 7. REPLICATION BEFORE CLAIMING

A finding from one cohort is a hypothesis. Test it in an independent dedicated AD cohort before it
is written as a result. **Report what fails to replicate in the same place as what succeeds** — the
FCER1A-low phenotype (§8 of the evidence doc) is the worked example.

## 8. REPORTING

- Lead with the finding. **Never bury a result in a subordinate clause** — if it is the answer, it is
  the headline, the abstract sentence, and a figure panel.
- State effect size, n, and the adjusted test. Bare P values are not results.
- Figures: Nimbus Sans; the muted palette in `analysis/palette.py`; minimal on-plot text; clean
  legends; legible at manuscript scale.
- Distinguish *measured*, *modelled*, and *inferred* explicitly.

## 9. RUN THE GUARDRAILS

```bash
python3 analysis/guardrails.py        # checks this contract against the repo, exits non-zero on violation
```

Every analysis script that makes a group comparison should import and use
`analysis/guardrails.py` (`require_mast_subject`, `comparison_order`, `check_model_terms`).

---

## Current standing answers (update when evidence changes)

*Last updated: 2026-08-15. Evidence: `results/REPORT.md`, tables in `results/tables/`.*

**Q0 — MECHANISM (in-vitro, primary cells).** TNFRSF9 is an ACTIVATION-INDUCED
gene in primary human mast cells. Resting human SKIN mast cells are essentially
negative (0.14 FPKM). IL-33 induces 29x, IgE/FcεRI cross-linking 7x, both
together 59x (8.66 FPKM); TSLP and IL-25 alone do nothing (GSE196862).
Reproduced in 4 paired donors: IL-33 log2FC +6.53 (P=0.001), IgE +5.51
(P=0.001), activated CD4+ Teff +2.97 (P=0.003), resting T cells no change
(GSE235240). Housekeeping and tryptase genes flat -> specific, not global.
TNFRSF9 has the lowest resting baseline and the largest induction of the
co-stimulatory TNFRSF members tested. **The two sufficient stimuli are the two
that define AD.** This explains the in-vivo heterogeneity (activation state at
biopsy) and why psoriasis matches AD at tissue level.

**Q1 — Is TNFRSF9 expressed in cutaneous mast cells in AD?**
**Yes, at low level.** 2.4–2.6% of AD mast cells carry a TNFRSF9 transcript
(0.36–0.47 per 10k mast UMI) in the two cohorts with adequate mast-cell
recovery. Reproducible. Mast cells nevertheless supply only 2.7–5.7% of the
TNFRSF9 in the tissue.

**Q2 — Does expression change between healthy and AD skin?**
**Higher in AD in 3/3 single-cell cohorts, but the magnitude does not
replicate.** Lead with the UNADJUSTED estimand (TNFRSF9 per mast cell): healthy
mast cells are the deeper arm (AD/H UMI ratio 0.74–1.04 across cohorts), so
unadjusted is the conservative choice and it also has lower heterogeneity.
Pooled **+1.24 (0.32–2.16) fixed, +1.32 (−0.07–2.71) random, I²=42%**. Adjusting
per mast UMI gives +1.44 / +1.53, I²=58%. Quote the random-effects estimate: a
~2.5-fold increase that does not reach conventional significance.

**On "depth".** Sequencing depth is a property of the run and is shared by every
cell in a library. Per-cell UMI is sampling effort = library reads x capture x
the cell's own mRNA content. Variance decomposition of log(UMI/cell): 27.4%
between libraries, 22.0% between cell types WITHIN a library (cannot be
sequencing), 50.6% cell-to-cell. Within-library mast:stromal UMI ratio is
0.08-0.30 in all 39 libraries and 34% lower in AD (0.134 vs 0.203, P<0.001).
AD libraries are sequenced deeper overall (4,292 vs 1,564 UMI/cell) yet AD mast
cells are shallower (693 vs 935) -- no flow-cell effect can do that. The deficit
is cell-intrinsic, so dividing by mast UMI removes biology. **Per-cell is the
primary estimand; per-UMI is reported beside it.**

**On depth adjustment.** The arms are NOT depth-matched: AD cells are 1.5–2.2x
deeper in fibroblasts, keratinocytes, macrophages and T/NK — but 0.74x in mast
cells. Consequences, both reported: (a) the mast result does not need the
adjustment (unadjusted 6.1x, P=0.030; raw TNFRSF9+ cell counts OR 5.35,
P=5.3e-4); (b) the §6 control populations DO need it — unadjusted, fibroblasts
(+2.05, P=0.011) and macrophages (+3.09, P=0.013) also "rise", which is the
shared-shift artefact §6 exists to catch.

**Decomposition (§4).** Abundance is unchanged (1.47% → 1.41% of cells,
log2OR +0.66, P=0.37); the per-cell factor carries the whole effect; the product
is directionally consistent but non-significant (log2FC +2.34, P=0.13).

**Three things NOT established.**
1. *De-novo induction.* Near-absence in healthy mast cells is a discovery-cohort
   property (3 transcripts, 5/6 donors zero) and **did not replicate** — REP2
   healthy mast cells are positive at 0.21 per 10k.
2. *Specificity to TNFRSF9.* 289 of 3,454 abundance-matched genes shift more
   (empirical P=0.084); the AD mast-cell transcriptome differs broadly.
3. *Specificity to AD.* Whole-skin TNFRSF9 rises as much in psoriasis
   (log2FC +2.20) as in AD (+2.48).

**Artefacts excluded.** Depth (stratified rate ratio 9.6×, i.e. depth suppresses
the effect); mast–T doublets (97/107 AD transcripts come from mast cells with no
T-cell transcript, log2FC +2.73, P=0.022); dissociation stress (healthy arm is
the *most* stressed, biasing against the finding); batch (fibroblasts and
keratinocytes show no shift).

**Binding limitation.** Mast-cell recovery ranges 0.11%–4.4% of cells across
cohorts and one dedicated AD atlas (GSE147424) yielded none at all. Every
number here is conditional on the dissociation protocol.
