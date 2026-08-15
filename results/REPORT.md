# TNFRSF9 (CD137/4-1BB) is an activation-induced gene in human cutaneous mast cells

### Absent at rest, switched on ~60-fold by IL-33 and IgE cross-linking, and modestly and heterogeneously elevated in atopic dermatitis skin

**Subject: cutaneous mast cells.** Fibroblasts, keratinocytes, T/NK cells and
macrophages appear in this report in exactly one role — as control populations,
to show that a mast-cell result is specific and not a batch or global effect
(CLAUDE.md §2).

Figures `results/figures/` · Tables `results/tables/` · Datasets and exclusions
`docs/datasets.md` · Contract compliance `python3 analysis/guardrails.py`

*Sections below are numbered 1–9. References of the form "CLAUDE.md §N" point to
the project contract, not to this document.*

---

## Abstract

Mast cells sit at the centre of type-2 skin inflammation, and TNFRSF9 (CD137,
4-1BB) is a co-stimulatory receptor of therapeutic interest, but whether
cutaneous mast cells express it — and whether that changes in atopic dermatitis
(AD) — has not been established. Using six public human cohorts spanning
single-cell, spatial, bulk and in-vitro transcriptomics (20 AD and 12 healthy
donors in vivo; primary human skin mast cells in vitro), we find that **TNFRSF9
is an activation-induced gene in mast cells rather than a disease marker.**
Resting primary human skin mast cells are essentially TNFRSF9-negative
(0.14 FPKM); IL-33 induces the gene 29-fold, IgE-receptor cross-linking 7-fold,
and the two together 59-fold, while TSLP and IL-25 alone do nothing. An
independent paired-donor cohort reproduces this (IL-33 log₂FC +6.53, P = 0.001,
n = 4 donors). In AD skin, mast-cell TNFRSF9 exceeds healthy skin in 3 of 3
single-cell cohorts, but the magnitude is heterogeneous and the pooled
random-effects estimate is not conventionally significant (log₂FC +1.32,
95% CI −0.07 to +2.71, P = 0.063, I² = 42%; 843 versus 4,315 mast cells).
Mast-cell abundance does not change, so the entire effect sits in per-cell
expression. Mast cells nevertheless contribute only 2.7–5.7% of the TNFRSF9 in
skin, and whole-tissue TNFRSF9 rises as much in psoriasis as in AD, so the
tissue-level signal marks inflamed skin rather than AD. The in-vitro mechanism
accounts for the in-vivo variance: TNFRSF9 reports how activated a mast cell was
at the moment of biopsy.

---

## 1. TNFRSF9 is expressed by cutaneous mast cells, at low level

In the discovery cohort (GSE204762: 11 AD donors with paired non-lesional and
lesional biopsies, 6 healthy donors, 280,518 cells, 39 libraries), 3,661 cells
passed the mast-cell marker QC of CLAUDE.md §6 (≥2 of TPSAB1/TPSB2/CPA3). Of
these, 94 (2.57%) carried at least one TNFRSF9 transcript. The replication cohort
with adequate mast-cell recovery (GSE153760) gave 33 of 1,379 (2.4%). Expression
is a reproducible property of the population, not an artefact of one dataset —
but it is a low-expression gene, and every analysis below is governed by counting
statistics rather than fold-change estimation (Fig. 1).

Mast cells are also the shallowest-sampled population in the atlas: a median of
546–622 UMI per cell against 3,157–4,061 for fibroblasts and keratinocytes
(Fig. 1e). Section 7.1 sets out what follows from that.

---

## 2. Healthy versus AD: the primary comparison (CLAUDE.md §3)

### 2.1 Abundance × per-cell expression, as one quantity (CLAUDE.md §4)

| Quantity | Healthy | AD | Effect | Test |
|---|---|---|---|---|
| **F1 · abundance** (mast % of all cells) | 1.47% | 1.41% | log₂OR +0.66 (−0.78 to 2.10) | binomial GLM, donor-clustered, P = 0.37 |
| **F2 · per-cell expression** (TNFRSF9 per mast cell) | 0.56 /100 cells | 3.43 /100 cells | log₂FC **+2.40 (6.1×)** | exact donor permutation, **P = 0.030** |
| **PRODUCT · F1 × F2** (mast TNFRSF9 per unit tissue) | 2.3×10⁻⁴ /10k | 1.2×10⁻³ /10k | log₂FC +2.34 (−0.72 to 5.39) | Poisson GLM, donor-clustered, P = 0.13 |

n = 6 healthy donors (538 mast cells), 11 AD donors (3,123 mast cells).

The two factors and their product are one quantity and are reported together.
**The per-cell factor carries the entire effect; abundance does not move.** The
product is directionally consistent but underpowered, inheriting the variance of
both factors.

Normalising per-cell expression by mast-cell UMI rather than cell count gives a
larger effect (0.057 → 0.469 per 10k, log₂FC +2.81, P = 0.012, n = 6 vs 11
donors). Section 7.1
explains why the per-cell number is the one quoted here.

### 2.2 Healthy versus AD non-lesional, and versus AD lesional

| Comparison | log₂FC (per mast UMI) | Exact permutation P | Transcripts | Donors |
|---|---|---|---|---|
| Healthy vs AD non-lesional | +3.05 | **0.0031** | 3 vs 58 | 6 vs 10 |
| Healthy vs AD lesional | +2.59 | 0.179 | 3 vs 49 | 6 vs 11 |

Both arms of AD skin sit above healthy skin. The lesional comparison misses
significance despite a similar effect size because lesional mast cells vary more
between donors — a power statement, not evidence of absence.

### 2.3 The refinement: AD lesional versus non-lesional

Within patients (10 donors with both samples), per-mast-cell TNFRSF9 did not
differ between lesional and non-lesional skin (log₂FC −0.46, 95% CI −1.97 to
1.06, P = 0.55). Mast-cell TNFRSF9 is a property of AD skin generally rather than
of the lesion — consistent with 2.2, where the non-lesional contrast is the
stronger of the two.

---

## 3. Replication, and what fails in it (CLAUDE.md §7)

Two independent cohorts, identical donor-level exact permutation test. Neither
deposit carries author cell-type labels, so mast cells were called by a rule
*calibrated on the discovery cohort's labels* — ≥2 of TPSAB1/TPSB2/CPA3 detected
**and** ≥50 tryptase+CPA3 transcripts per 10k UMI (precision 96.0%, sensitivity
79.3%). The magnitude requirement is not cosmetic: detection alone tags hundreds
of keratinocytes carrying ambient tryptase from lysed mast cells (Fig. 5a).

| Cohort | Donors H/AD | Mast cells H/AD | Transcripts H/AD | log₂FC per cell | log₂FC per UMI |
|---|---|---|---|---|---|
| Discovery GSE204762 (3′) | 6 / 11 | 538 / 3,123 | 3 / 107 | +2.40 | +2.81 |
| REP1 GSE222840+GSE173205 (5′) | 4 / 5 | 27 / 91 | 0 / 2 | +0.57 | +0.51 |
| REP2 GSE153760 biopsies (3′ v3) | 2 / 4 | 278 / 1,101 | 6 / 38 | +0.58 | +0.68 |
| **Pooled, fixed effects** | 12 / 20 | 843 / 4,315 | 9 / 147 | **+1.24** (0.32–2.16), P = 0.0085 | +1.44 (0.52–2.37) |
| **Pooled, random effects** | 12 / 20 | 843 / 4,315 | 9 / 147 | **+1.32** (−0.07–2.71), P = 0.063 | +1.53 (−0.16–3.22) |
| Heterogeneity | | | | **I² = 42%** | I² = 58% |

**What replicates:** the direction, 3 of 3 cohorts; and the fact of expression in
every cohort with usable mast-cell recovery.

**What does not:** the *magnitude* (the discovery effect is 4–6× either
replication), and statistical significance — no replication cohort reaches it
alone and neither does the random-effects pooled estimate.

**Why REP1 is uninformative rather than negative.** Mast-cell recovery failed in
that pipeline: 118 mast cells in 111,370 (0.11%), against 1.2–1.5% in discovery
and 4.4% in REP2, with no coherent mast-cell cluster in the embedding (Fig. 5a).
Two TNFRSF9 transcripts in the entire cohort cannot confirm or refute anything.

**Concentration warning for both.** REP2's AD signal is concentrated in one donor
(AD7: 31 of 38 transcripts, and an outlier at 13.3% mast cells), and its healthy
arm has 2 donors, one contributing 5 of 6 transcripts. The discovery cohort is
similarly concentrated (MGH108: 55 of 107). No in-vivo cohort examined here has
the donor numbers this question needs.

---

## 4. Mechanism: TNFRSF9 is induced by the two signals that define AD

No tissue analysis can separate "mast cells transcribe more TNFRSF9 in AD" from
"AD skin holds a different mast-cell population". A stimulation experiment can,
and two public datasets perform exactly that on **primary** human mast cells —
with no dissociation artefact, no ambient RNA, no doublets and no depth confound
(Fig. 6).

### 4.1 Primary human SKIN mast cells + the AD milieu (GSE196862)

Skin-derived mast cells, 24 libraries, stimulated with IgE-receptor cross-linking
and the epithelial alarmins IL-33, TSLP and IL-25, alone and combined.

| Condition | TNFRSF9 (FPKM) | log₂FC vs resting |
|---|---|---|
| **resting** | **0.14** | — |
| TSLP | 0.11 | −0.28 |
| IL-25 | 0.11 | −0.32 |
| IgE/Ag | 1.07 | +2.87 (7×) |
| IL-33 | 4.33 | +4.88 (29×) |
| IL-33 + IL-25 + TSLP | 3.76 | +4.67 |
| **IgE/Ag + IL-33** | **8.66** | **+5.88 (59×)** |
| IgE/Ag + all three alarmins | 7.96 | +5.75 |

Resting skin mast cells are effectively TNFRSF9-negative (replicates 0.011,
0.047, 0.006, 0.487 FPKM). IL-33 and FcεRI cross-linking each induce the gene and
synergise; TSLP and IL-25 alone do nothing, so this is not a generic cytokine
response.

**Specific, not global.** In the same libraries housekeeping genes are flat
(B2M +0.34, RPL13A −0.23, TMSB4X +0.19 under IgE+IL-33) and mast identity genes
barely move (TPSAB1 −0.29, TPSB2 +0.02). Among co-stimulatory TNF-receptor
relatives, TNFRSF9 has both the lowest resting level (0.14 FPKM versus 1.6–1.8
for TNFRSF18 and TNFRSF4) and the largest induction (+5.88 versus +3.20 and
+0.73): the most switch-like member of the family in this cell type.

### 4.2 Primary human mast cells + activated CD4⁺ T cells (GSE235240)

Four donors, paired, mast cells FACS-sorted after 24 h co-culture.

| Condition | TNFRSF9 (normalised counts) | paired log₂FC | P | n donors |
|---|---|---|---|---|
| resting | 103 | — | — | 4 |
| resting CD4⁺ T cells | 133 | +0.26 | 0.17 | 4 |
| **activated CD4⁺ T cells** | 768 | **+2.97** | **0.003** | 4 |
| IgE/Ag | 4,148 | +5.51 | 0.001 | 4 |
| **IL-33** | **9,079** | **+6.53 (92×)** | **0.001** | 4 |

The ordering reproduces with an independent stimulus and a paired design:
IL-33 > IgE/Ag > activated T cells > resting T cells ≈ unstimulated. Tryptase is
unchanged (TPSAB1 log₂FC −0.04, P = 0.61, n = 4 donors), so this is not a
differentiation or viability effect.

### 4.3 What this settles

**The mechanism.** TNFRSF9 in mast cells is switched on from a near-zero baseline
by exactly the two signals that define AD: IL-33 released by damaged
keratinocytes, and IgE cross-linking on FcεRI. This is a *measured* result in
primary cells, not an inference from tissue.

**Why the in-vivo data look as they do.** If TNFRSF9 reports activation state
rather than disease identity, then healthy mast cells should be near-zero (they
are); AD mast cells should be positive to a degree depending on how activated
they were at biopsy, predicting the observed between-cohort heterogeneity
(Section 3); and the tissue-level rise should not be AD-specific, because any
mast-activating inflammation would do the same (Section 6).

**What it does not settle.** These are in-vitro stimulations, not AD skin. They
establish that mast cells can express TNFRSF9 on demand and identify sufficient
stimuli. They do not establish how much of the in-vivo signal runs through this
pathway, nor that 4-1BB protein reaches the cell surface, nor any function.

---

## 5. In situ, without dissociation (GSE197023 Visium)

The per-cell estimate is conditional on which mast cells survive enzymatic
digestion (CLAUDE.md §4). Visium never dissociates the tissue and reproduces the
direction: whole-spot TNFRSF9 was higher in AD than healthy skin (log₂FC +2.06,
95% CI 0.35–3.77, P = 0.019; 13,742 spots, 6 healthy versus 13 AD sections,
donor-clustered with a log-depth offset), driven by lesional skin (+2.30,
95% CI 0.47–4.12, P = 0.014). Mast-cell content in situ was unchanged or slightly
lower (TPSAB1 log₂FC −0.80, P = 0.34; CPA3 −0.71, P = 0.021, 6 vs 13 sections),
matching the unchanged abundance in the single-cell arm.

Spot depth differs ~4-fold between arms (AD lesional 4,025 versus healthy 1,220
UMI), so every section in Fig. 3 is annotated with its own median depth.

**What does not hold spatially.** TNFRSF9 tracked mast-cell content in AD
*non-lesional* skin (log₂FC +1.13 per SD of tryptase content, 95% CI 0.40–1.87,
P = 0.0025, 3,727 spots, 6 donors) but **not** in AD lesional skin (+0.04,
95% CI −0.56 to 0.64, P = 0.90) or healthy skin (−0.05, P = 0.96). The mandated
control contents behaved as controls should in the same spots (fibroblast +0.73,
P = 0.13; keratinocyte −0.13, P = 0.50). In lesional skin, where whole-tissue
TNFRSF9 is highest, the signal is not spatially organised around mast cells —
consistent with Section 6.

---

## 6. Mast cells are a minor source of skin TNFRSF9, and the tissue signal is not AD-specific

| Arm | TNFRSF9 in mast cells | in whole tissue | mast share |
|---|---|---|---|
| Healthy | 3 | 104 | 2.9% |
| AD non-lesional | 58 | 1,018 | 5.7% |
| AD lesional | 49 | 1,844 | 2.7% |

Mast cells account for 0.22–0.41% of the tissue transcriptome and 2.7–5.7% of its
TNFRSF9. Induction within mast cells and their small share of the total are both
true; a receptor switched on in a rare cell can matter functionally without
dominating a bulk measurement, but neither can it explain a tissue-level signal.

Whole-skin RNA-seq (GSE121212, 147 samples) makes the specificity problem
explicit:

| Comparison | log₂FC | 95% CI | n | P |
|---|---|---|---|---|
| Healthy vs AD | +2.48 | 1.85 to 3.11 | 38 vs 54 | 1.0×10⁻¹⁴ |
| Healthy vs AD non-lesional | +1.38 | 0.75 to 2.01 | 38 vs 27 | 1.6×10⁻⁵ |
| Healthy vs AD lesional | +3.06 | 2.37 to 3.75 | 38 vs 27 | 2.9×10⁻¹⁸ |
| AD lesional vs non-lesional | +1.68 | 1.04 to 2.31 | 27 vs 27 | 2.4×10⁻⁷ |
| **Healthy vs psoriasis** | **+2.20** | 1.80 to 2.59 | 38 vs 55 | 8.7×10⁻²⁸ |

Poisson GLM on raw counts, log-library-size offset, errors clustered on patient.
Mast-cell content in the same samples was flat (TPSAB1 log₂FC −0.09, P = 0.74;
TPSB2 −0.10, P = 0.73; CPA3 −0.48, P = 0.037, n = 38 vs 54). In bulk tissue,
abundance and per-cell expression are not separable at all, and the mast-cell
contribution is not separable from any other source. What bulk establishes is
that skin TNFRSF9 rises steeply in AD, does so without any rise in mast-cell
content, and rises just as much in psoriasis — a marker of inflamed skin.

---

## 7. Rigour: estimand, controls, artefacts, calibration

### 7.1 What "depth" means, and which estimand follows

Sequencing depth is a property of the run and is shared by every cell in a
library. Per-cell UMI is not that: it is library reads × capture efficiency ×
**the cell's own mRNA content**. Only the first term is sequencing, and
normalising by per-cell UMI divides out the third, which is biology.

Variance decomposition of log(UMI per cell) across mast cells and the control
populations:

| Component | Share | Can it be sequencing depth? |
|---|---|---|
| Between libraries (flow cell, loading) | 27.4% | Yes |
| Between cell types **within** a library | 22.0% | **No** — those cells shared a flow cell |
| Cell to cell within a cell type | 50.6% | No |

Within-library ratios make this concrete: mast cells yield 0.08–0.30 of the UMI
that fibroblasts and keratinocytes yield **in the same run**, in all 39 libraries,
and that ratio is 34% lower in AD (median 0.134, n = 32 libraries) than healthy
(0.203, n = 7 libraries; Mann–Whitney P < 0.001). Decisively, the two move in
opposite directions: **AD libraries are sequenced more deeply overall (median
4,292 versus 1,564 UMI per cell) yet AD mast cells are shallower (693 versus
935).** No flow-cell effect can do that. AD mast cells are transcriptionally
smaller — 399 versus 477 genes detected, less tryptase per cell.

**Consequence.** Dividing by mast-cell UMI removes biology along with technical
effort, using a denominator 25% smaller in AD. The **per-cell estimand is
primary** here; the per-UMI estimand is reported beside it as the compositional
companion, as CLAUDE.md §6 requires. The result does not depend on the choice:

| Estimand | Healthy | AD | Fold | P |
|---|---|---|---|---|
| per mast UMI (adjusted) | 0.057 /10k | 0.469 /10k | 8.2× | 0.012 |
| **per mast CELL (unadjusted)** | 0.56 /100 | 3.43 /100 | **6.1×** | **0.030** |
| **no normalisation** (TNFRSF9⁺ cell counts) | 3/538 | 91/3,123 | OR 5.35 | **5.3×10⁻⁴** |

**But the adjustment cannot be dropped for the control populations**, which *are*
~2× deeper in AD. Unadjusted, TNFRSF9 appears to rise in almost every cell type —
the "shared shift means batch" pattern CLAUDE.md §6 exists to catch:

| Population (6 healthy vs 11 AD donors) | log₂FC per UMI (adjusted) | log₂FC per cell (unadjusted) |
|---|---|---|
| **Mast (subject)** | +2.81 (P = 0.012) | +2.40 (P = 0.030) |
| Fibroblasts | +0.91 (P = 0.23) | **+2.05 (P = 0.011)** |
| Keratinocytes | +0.26 (P = 0.72) | +1.37 (P = 0.15) |
| T/NK | +1.49 (P = 0.0001) | +2.06 (P = 0.0013) |
| Macrophages | +2.06 (P = 0.024) | +3.09 (P = 0.013) |

The adjustment is unnecessary for the mast-cell effect and essential for the
claim that the effect is specific. Both are reported.

### 7.2 The mandated controls (CLAUDE.md §6)

| Control | Result | Consequence |
|---|---|---|
| Library size | log-depth offset throughout; depth-stratified analysis below | Effect survives and strengthens |
| Donor aggregation | All inference donor-level; exact permutation over all 12,376 donor-label assignments | Cluster-robust SEs on 6 clusters not relied upon |
| Same test in fibroblasts + keratinocytes | +0.91 (P = 0.23) and +0.26 (P = 0.72), 6 vs 11 donors | **No shared shift → not batch** |
| Per-group median depth | Mast 546 / 541 / 622 UMI; tissue 1,999 / 2,833 / 3,306 (healthy / AD NL / AD LS) | Mast depth comparable across arms |
| Marker QC | 3,661 of 4,313 deposited "Mast" labels passed (84.9%); 433 QC-positive cells carried non-mast labels and were excluded | Deposited labels imperfect, as anticipated |
| Dissociation stress | HSP+IEG 986 per 10k mast UMI healthy versus 416 AD | **Healthy is the most stressed arm — biases against the finding**; Spearman(stress, rate) = −0.15, P = 0.57 |

**Depth-matched detection.** Healthy mast cells contributed zero TNFRSF9
transcripts in every depth stratum below 800 UMI, where AD mast cells contributed
35. The depth-stratified common rate ratio is **9.6×**, larger than the crude
estimate (Fig. 1f).

**Leave-one-donor-out.** Dropping any one of the 17 discovery donors (6 healthy,
11 AD) leaves log₂FC between +2.15 and +4.72, permutation P ≤ 0.043 throughout.
Removing the most influential donor (MGH108) gives +2.15 (P = 0.023).

### 7.3 The doublet artefact

T/NK cells carry TNFRSF9 at a rate comparable to AD mast cells, so a mast–T
doublet would present as exactly what is being counted. TNFRSF9⁺ mast cells do
carry T-cell transcripts more often than TNFRSF9⁻ ones (10.6% versus 3.4%;
depth-stratified OR 2.68 across six strata, n = 94 versus 3,567 cells), so some
contamination exists. No cell was flagged by the deposited scrublet call, and
TNFRSF9⁺ mast cells had *lower* doublet scores (median 0.02 versus 0.03).

The decisive test is deletion:

| Cell set | Mast cells | TNFRSF9 H / AD | log₂FC | Permutation P |
|---|---|---|---|---|
| All marker-QC mast cells | 3,661 | 3 / 107 | +2.81 | 0.012 |
| **Excluding any cell with a T-cell transcript** | 3,531 | 3 / 97 | **+2.73** | **0.022** |
| Excluding any T-cell *or* myeloid transcript | 3,187 | 3 / 78 | +2.40 | 0.050 |

97 of 107 AD transcripts come from mast cells with **no T-cell transcript at
all**, so doublets cannot be the source. The stricter purge removes 13% of mast
cells and pushes the test to the boundary: the effect is not eliminated, but this
cohort has no power to spare.

### 7.4 Calibration: is the effect exceptional?

* **Is the test calibrated?** Permuting arm labels *within* AD donors, where no
  disease contrast exists by construction, gives a median 4.0% of matched genes
  at P ≤ 0.05 (IQR 2.9–5.6%, 40 random splits). The exact permutation test is
  approximately calibrated.
* **Is TNFRSF9 exceptional?** No. The real contrast yields 17.8% of matched genes
  at P ≤ 0.05 — far outside the null range — so the mast-cell transcriptome
  differs widely between healthy and AD skin. Against that background TNFRSF9
  ranks 290 of 3,454 (empirical P = 0.084, Fig. 2e).

TNFRSF9 induction in AD mast cells is a real, calibrated nominal result in the
top 8% of a broad disease-associated shift — but it is *one of many* genes
moving, and the claim that TNFRSF9 is *specifically* dysregulated relative to the
rest of the mast-cell transcriptome is **not** supported.

---

## 8. Measured, modelled, inferred (CLAUDE.md §8)

* **Measured.** Transcript counts per cell, per spot, per bulk library and per
  in-vitro condition; marker expression; sequencing depth; the in-vitro induction
  of TNFRSF9 by IL-33 and IgE in primary skin mast cells.
* **Modelled.** All fold changes and intervals (Poisson/binomial GLMs with
  log-depth offsets and donor- or patient-clustered errors; exact donor-label
  permutation; inverse-variance meta-analysis).
* **Inferred.** That the in-vivo per-mast-cell difference reflects the same
  activation programme demonstrated in vitro; and that it is transcriptional
  induction rather than selective survival of a TNFRSF9⁺ subset through
  dissociation. The Visium arm supports the second indirectly but cannot resolve
  single mast cells at 55 µm.
* **Not established.** 4-1BB protein on the mast-cell surface (see
  `results/PROTEIN_EVIDENCE.md`); any functional consequence; AD-specificity of
  the tissue-level change; the in-vivo magnitude (I² = 42%); and selective
  dysregulation relative to the rest of the mast-cell transcriptome.

---

## 9. Limitations

1. **Small counts.** Nine TNFRSF9 transcripts across all healthy mast cells in
   all three in-vivo cohorts. Statements about healthy skin are upper bounds.
2. **Donor numbers.** 12 healthy and 20 AD donors in vivo, with each cohort's
   signal concentrated in one or two donors.
3. **Mast-cell recovery governs everything.** Recovery ranged 0.11%–4.4% of cells
   across cohorts under one rule. One dedicated AD atlas (GSE147424,
   cryopreserved biopsies) contained *no* recoverable mast cells and was excluded
   (`docs/datasets.md`); suction-blister sampling also fails to capture them.
4. **Heterogeneity.** I² = 42%; the random-effects pooled estimate does not reach
   significance.
5. **Not exceptional against background** (7.4), and **not AD-specific** at
   tissue level (Section 6).
6. **In vitro is sufficiency, not contribution.** Section 4 shows mast cells
   *can* be driven to express TNFRSF9; it does not quantify how much of the AD
   signal runs through that route.
7. **Transcript, not protein.** Everything above is RNA. Protein-level evidence
   is assessed separately in `results/PROTEIN_EVIDENCE.md`.
