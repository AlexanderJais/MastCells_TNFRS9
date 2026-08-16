# TNFRSF9 (CD137/4-1BB) is an activation-induced gene in human cutaneous mast cells

### Absent at rest, switched on ~60-fold by IL-33 and IgE cross-linking, and modestly and heterogeneously elevated in atopic dermatitis skin

**Subject: cutaneous mast cells.** Fibroblasts, keratinocytes, T/NK cells and
macrophages appear in this report in exactly one role — as control populations,
to show that a mast-cell result is specific and not a batch or global effect
(CLAUDE.md §2).

Figures `results/figures/` (Fig. 1–6, Fig. S1) · Tables `results/tables/` · Datasets and exclusions
`docs/datasets.md` · Contract compliance `python3 analysis/guardrails.py`

*Sections below are numbered 1–9. References of the form "CLAUDE.md §N" point to
the project contract, not to this document.*

*Figure conventions: error bars are **± 1 SEM** throughout. Intervals quoted in
the text are 95% confidence intervals and are labelled as such.*

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
n = 4 donors). In AD skin the discovery cohort shows a clear increase — 0.56% of
healthy mast cells carry a TNFRSF9 transcript versus 2.91% in AD (Fisher OR 5.35,
P = 0.001) — but **this does not reproduce in the one independent cohort with
adequate mast-cell recovery** (2.16% versus 2.45%, OR 1.14, P = 1.00), where the
healthy baseline is four-fold higher. Mast-cell abundance is unchanged in both,
so nothing in this system is driven by cell number. Mast cells nevertheless contribute only 2.7–5.7% of the TNFRSF9 in
skin, and whole-tissue TNFRSF9 rises as much in psoriasis as in AD, so the
tissue-level signal marks inflamed skin rather than AD. The in-vitro mechanism
supplies the likely explanation for the discordant in-vivo cohorts: TNFRSF9
reports how activated a mast cell was when the biopsy was taken, which is not a
fixed property of a diagnosis.

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
541–891 UMI per cell across the three tissue groups, against 1,006–3,373 for
fibroblasts and 1,268–4,061 for keratinocytes
(`results/tables/sc_depth_by_arm_celltype.csv`). Section 7.1 sets out what
follows from that, and why per-cell UMI is *not* "sequencing depth".

---

## 2. Healthy versus AD: the primary comparison (CLAUDE.md §3)

### 2.1 Mast-cell numbers and TNFRSF9⁺ mast-cell numbers

Each point in Fig. 2 is one sample (one donor × arm); statistics are donor-level
Mann–Whitney against the healthy arm, plus a pooled Fisher exact test on raw
TNFRSF9⁺ cell counts.

| Quantity | Healthy | AD non-lesional | AD lesional | MWU vs healthy (NL / LS) |
|---|---|---|---|---|
| Samples (donors) | 6 (6) | 10 (10) | 11 (11) | — |
| Cells profiled | 59,745 | 97,947 | 122,826 | — |
| **Mast cells recovered** | **538** | **1,564** | **1,559** | — |
| Mast cells per sample, median | 99.5 | 131 | 151 | P = 0.37 / 0.26 |
| **Mast cells, % of all cells** (median) | **1.47%** | **1.30%** | **1.21%** | **P = 0.79 / 1.00** |
| **TNFRSF9⁺ mast cells** | **3** | **52** | **39** | — |
| TNFRSF9⁺ per sample, median | 0 | 2.5 | 2.0 | P = 0.036 / 0.106 |
| **TNFRSF9⁺ as % of mast cells** | **0.56%** | **3.33%** | **2.50%** | **P = 0.020 / 0.090** |

**Pooled healthy versus AD: 3 of 538 versus 91 of 3,123 mast cells are TNFRSF9⁺
(Fisher OR 5.35, P = 5.3×10⁻⁴).**

Read together, as CLAUDE.md §4 requires: **mast-cell abundance does not differ
between healthy and AD skin (1.47% versus 1.30%/1.21% of cells, P = 0.79/1.00),
while the fraction of mast cells carrying TNFRSF9 rises about 5-fold (0.56% to
2.50–3.33%).** The change is in what mast cells express, not in how many there
are. Five of six healthy donors contributed no TNFRSF9⁺ mast cell at all.

A per-transcript version of the same contrast — TNFRSF9 molecules per mast cell,
which uses all 110 transcripts rather than only presence/absence — gives log₂FC
+2.40 (6.1×, exact donor permutation P = 0.030, n = 6 vs 11 donors); normalising
by mast-cell UMI instead gives +2.81 (P = 0.012). Section 7.1 explains why the
per-cell version is the one quoted.

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

## 3. Replication: the AD increase does not reproduce (CLAUDE.md §7)

**Which cohorts can test this.** A cohort contributes only if it recovered ≥100
marker-QC mast cells in each arm and ≥0.5% of its cells as mast cells — a
mast-cell QC threshold judged on tryptase/CPA3 marker data alone, independent of
TNFRSF9 and of the outcome.

| Cohort | Mast-cell recovery | Mast cells H / AD | Qualifies |
|---|---|---|---|
| GSE153760 | 4.41% | 278 / 1,101 | yes |
| GSE204762 (discovery) | 1.31% | 538 / 3,123 | yes |
| GSE222840 + GSE173205 | **0.106%** | **27 / 91** | **no** |

GSE222840+GSE173205 recovered mast cells 12–40× less efficiently than the other
two and yielded 2 TNFRSF9 transcripts in 111,370 cells (Fig. 5d). A cohort that
cannot recover the cell type cannot test a gene in it, so it is reported as an
attempted replication that failed on mast-cell recovery, not as evidence about
TNFRSF9, and it is not pooled.

**The result in the one qualifying independent cohort.**

| Cohort | TNFRSF9⁺ healthy | TNFRSF9⁺ AD | Fisher OR | P |
|---|---|---|---|---|
| Discovery GSE204762 | 3/538 (**0.56%**) | 91/3,123 (**2.91%**) | 5.35 | **0.001** |
| **GSE153760** | 6/278 (**2.16%**) | 27/1,101 (**2.45%**) | **1.14** | **1.00** |

**The AD increase does not reproduce.** In GSE153760 the TNFRSF9⁺ mast-cell
fraction is the same in healthy and AD skin. The discrepancy is not in the AD
arms, which agree closely (2.91% versus 2.45%); it is in the healthy arms, which
differ four-fold (0.56% versus 2.16%). GSE153760's healthy baseline leaves no
room for an increase.

Neither cohort is well placed to settle which healthy baseline is right:
GSE153760 has **two** healthy donors, one contributing 5 of its 6 TNFRSF9⁺ cells;
the discovery cohort has six, five of them contributing none. On the donor-level
test used throughout this report, GSE153760's 2-versus-4 donor design admits only
C(6,4) = 15 label assignments, so its exact permutation P cannot fall below
0.067 whatever the effect size: this cohort can fail to support an increase, but
it could not have confirmed one at the 5% level.

**The two cohorts do not select mast cells by the same rule, and that is not the
explanation.** The discovery deposit carries author cell-type labels and its mast
cells are those labelled Mast that also pass marker QC; the replication deposits
carry no labels, so mast cells there are defined by marker QC plus a tryptase
magnitude threshold calibrated on the discovery cohort (96.0% precision, 79.3%
sensitivity against its labels). Re-analysing the discovery cohort under the
replication's rule leaves the discovery increase intact and slightly larger — 2
of 559 (0.36%) healthy versus 82 of 3,002 (2.73%) AD mast cells, Fisher OR 7.82,
P = 1.1×10⁻⁴, against OR 5.35 under its own rule
(`results/tables/replication_rule_sensitivity.csv`). The disagreement between
cohorts is therefore about the cohorts, not about the cell-selection rule.

On the per-transcript rate the two cohorts agree in direction (discovery log₂FC
+2.40, GSE153760 +0.58), and pooled they give +1.27 (95% CI 0.33 to 2.22,
P = 0.0084) under fixed effects and +1.43 (95% CI −0.36 to 3.21, P = 0.12) under
random effects. But the primary quantification — what fraction of mast cells
carry the transcript — replicates only as *expression*, not as an *AD-associated
increase*.

**What replicates:** that cutaneous mast cells express TNFRSF9, at 2.2–2.9% of
mast cells in both qualifying cohorts.

**What does not:** the increase in AD.

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

The §6 control genes were carried through the identical group test: COL1A1
(fibroblast) and KRT14 (keratinocyte) showed no comparable rise, so the spatial
result is not a global shift (`results/tables/spatial_group_glm.csv`).

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
TPSB2 −0.10, P = 0.73; CPA3 −0.48, P = 0.037, n = 38 vs 54). Nor does whole-skin
TNFRSF9 track mast-cell content *within* any arm (Spearman ρ = +0.20 healthy,
+0.26 AD non-lesional, +0.12 AD lesional, −0.13 psoriasis lesional; all P > 0.18,
n = 27–38 samples per arm) — consistent with mast cells being a minor source
(Section 6 above). In bulk tissue,
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
(0.203, n = 7 libraries; Mann–Whitney P < 0.001). Because it is formed within a
library, this ratio cannot be produced by a flow-cell effect, and it is the one
statement here that does not depend on how libraries are weighted.

Absolute UMI per cell does depend on that choice, and both weightings are
reported (`results/tables/depth_arm_estimators.csv`, 6 healthy versus 11 AD
donors, 39 libraries):

| Estimator | All cells, AD / healthy | Mast cells, AD / healthy |
|---|---|---|
| Pooled over cells | 4,111 / 2,152 (1.9×) | 693 / 935 (**0.74×**) |
| Median of library means | 4,292 / 1,564 (2.7×) | 629 / 347 (1.8×) |

AD libraries are the deeper ones under either weighting. The mast-cell row is
not stable: pooled over cells AD mast cells are shallower, but weighting each
library equally they are deeper, because the healthy mast cells are concentrated
in a few deep libraries. We therefore do not argue from "the two move in
opposite directions"; the within-library ratio above is the estimator-free
evidence that the deficit is cell-intrinsic. AD mast cells are also
transcriptionally smaller — 399 versus 477 genes detected, less tryptase per
cell.

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
| Per-group median depth | Median UMI per cell, healthy / AD NL / AD LS: mast 891 / 541 / 622; all cells 1,155 / 2,467 / 2,992 | **Healthy mast cells are the deepest arm** while healthy tissue is the shallowest — the mast comparison is not depth-favoured in AD |
| Marker QC | 3,661 of 4,313 deposited "Mast" labels passed (84.9%); 433 QC-positive cells carried non-mast labels and were excluded | Deposited labels imperfect, as anticipated |
| Dissociation stress | HSP+IEG 986 per 10k mast UMI healthy versus 416 AD | **Healthy is the most stressed arm — biases against the finding**; Spearman(stress, rate) = −0.15, P = 0.57 |

**Sampling-matched detection.** Stratifying mast cells by UMI per cell, healthy
mast cells contributed zero TNFRSF9 transcripts in every stratum below 800 UMI,
where AD mast cells contributed 35. The stratified common rate ratio is **9.6×**,
larger than the crude estimate (`results/tables/sc_depth_matched.csv`).

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
  ranks 290 of 3,454 (empirical P = 0.084, Fig. S1b).

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
  the tissue-level change; the in-vivo magnitude (I² = 70% per cell, 78% per UMI
  across the two qualifying cohorts); and selective dysregulation relative to the
  rest of the mast-cell transcriptome.

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
4. **The AD increase is not replicated.** It is clear in the discovery cohort
   (OR 5.35, P = 0.001) and absent in the one qualifying independent cohort
   (OR 1.14, P = 1.00). The disagreement is in the healthy baseline (0.56% vs
   2.16% of mast cells), which neither cohort has the donors to settle.
5. **Not exceptional against background** (7.4), and **not AD-specific** at
   tissue level (Section 6).
6. **In vitro is sufficiency, not contribution.** Section 4 shows mast cells
   *can* be driven to express TNFRSF9; it does not quantify how much of the AD
   signal runs through that route.
7. **Transcript, not protein.** Everything above is RNA. Protein-level evidence
   is assessed separately in `results/PROTEIN_EVIDENCE.md`.
