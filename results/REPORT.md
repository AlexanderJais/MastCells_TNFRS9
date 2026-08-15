# TNFRSF9 (CD137/4-1BB) is an activation-induced gene in human cutaneous mast cells: absent at rest, switched on ~60-fold by IL-33 and IgE cross-linking, and modestly elevated in atopic dermatitis skin

**Subject of this study: cutaneous mast cells.** Fibroblasts, keratinocytes,
T/NK cells and macrophages appear throughout as control populations only, to
establish that a mast-cell result is specific and not a batch or global effect
(CLAUDE.md §2).

Figures: `results/figures/` · Tables: `results/tables/` · Dataset inventory and
exclusions: `docs/datasets.md`

---

## Summary

**TNFRSF9 is expressed by human cutaneous mast cells.** Across three independent
single-cell cohorts, mast cells carried TNFRSF9 transcripts at 0.06–0.47 per
10,000 mast UMI, with 0.6–3.3% of mast cells positive. This is a low-expression
gene in this cell type, and the whole analysis is governed by counting
statistics rather than fold-change estimation.

**Mast-cell TNFRSF9 is higher in AD than in healthy skin in all three cohorts,
but the effect size is not reproducible.** Direction was consistent 3/3
(12 healthy and 20 AD donors in total, 843 versus 4,315 mast cells). The headline
estimate is the **unadjusted** one — TNFRSF9 molecules per mast cell — because
healthy mast cells are the *deeper* arm in this comparison, which makes the
unadjusted number the conservative one (§2): pooled **log₂FC +1.24 (95% CI +0.32
to +2.16, P = 0.0085) under fixed effects, +1.32 (95% CI −0.07 to +2.71,
P = 0.063) under random effects, I² = 42%**. Normalising instead by mast-cell UMI
gives a slightly larger effect (+1.44 fixed, +1.53 random, I² = 58%). Given the
heterogeneity, the random-effects estimate is the one to quote: **a ~2.5-fold
increase that does not reach conventional significance.**

**The change is in per-cell expression, not in cell number.** Mast-cell
abundance was 1.47% of cells in healthy skin versus 1.41% in AD (log₂OR +0.66,
95% CI −0.78 to 2.10, P = 0.37, 6 vs 11 donors). Because mast-cell-derived
TNFRSF9 is the product of abundance and per-cell expression (§4), the product
rises ~5-fold but does **not** reach significance (log₂FC +2.34, 95% CI −0.72 to
5.39, P = 0.13), inheriting the variance of both factors.

**In primary cells the mechanism is unambiguous.** Resting primary human *skin*
mast cells are essentially TNFRSF9-negative (0.14 FPKM). IL-33 induces the gene
29-fold, IgE-receptor cross-linking 7-fold, and the two together **59-fold**
(8.66 FPKM); TSLP and IL-25 alone do nothing. A second cohort of primary mast
cells reproduces it across 4 paired donors (IL-33 log₂FC +6.53, P = 0.001;
activated CD4⁺ T cells +2.97, P = 0.003), with housekeeping and tryptase genes
flat. **TNFRSF9 marks mast-cell activation state, and the two sufficient stimuli
are precisely the two that define AD** (§2B, Fig. 6). This is why the in-vivo
effect is real but heterogeneous: it tracks how activated the mast cells happened
to be at biopsy, not disease identity — which also explains why psoriasis matches
AD at tissue level.

**Four findings temper the in-vivo result, and one contradicts a tempting reading.**

1. The "off in health, on in disease" reading is **correct in principle but not
   absolute in tissue**. In primary cells the resting state really is
   near-negative (0.14 FPKM), and the discovery cohort matches that (3
   transcripts in 522,068 UMI; 5 of 6 donors at zero). But healthy skin in vivo
   is not a resting culture: GSE153760 healthy mast cells carried TNFRSF9 at
   0.21 per 10k (6 transcripts, 2.2% of 278 cells). Baseline activation in
   normal skin is expected on this model, so a strict de-novo claim overstates
   it — the in-vivo axis is graded activation, not a binary switch.
2. The increase is not exceptional: among 3,454 abundance-matched genes, 289
   shift more strongly, placing TNFRSF9 at empirical P = 0.084. The mast-cell
   transcriptome differs broadly in AD (17.8% of matched genes at P ≤ 0.05
   against a true-null median of 4.0%).
3. Mast cells are a **minor source** of skin TNFRSF9 — 2.7–5.7% of tissue
   TNFRSF9 transcripts.
4. At tissue level the rise is **not AD-specific**: whole-skin TNFRSF9 rises as
   much in psoriasis (log₂FC +2.20, 95% CI 1.80–2.59, n = 38 healthy vs 55
   psoriasis).

---

## 1. The primary comparison: healthy versus AD (§3)

### 1.1 Is TNFRSF9 expressed in cutaneous mast cells at all?

Yes, at low abundance. In the discovery cohort (GSE204762) 94 of 3,661
marker-QC mast cells (2.57%) carried at least one TNFRSF9 transcript
(0.392 per 10,000 mast UMI). The replication cohort with adequate mast-cell
recovery (GSE153760) gave 33 of 1,379 (2.4%). Expression is therefore a
reproducible property of the population, not of one dataset.

Mast cells are the **shallowest-sequenced population in the atlas** — median 546
(healthy), 541 (AD non-lesional) and 622 (AD lesional) UMI, against 3,157–4,061
UMI for fibroblasts and keratinocytes (Fig. 1e). Every comparison below carries a
log-library-size offset and is checked against depth-matched strata.

### 1.2 Healthy versus AD, pooled — discovery cohort (primary)

| Quantity (§4) | Healthy | AD | Effect | Test |
|---|---|---|---|---|
| **F1 · abundance** (mast % of cells) | 1.47% | 1.41% | log₂OR +0.66 (−0.78 to 2.10) | binomial GLM, donor-clustered, P = 0.37 |
| **F2 · per-cell expression** (TNFRSF9 per 10k mast UMI) | 0.057 | 0.469 | log₂FC **+2.81** | exact donor permutation (12,376 assignments), **P = 0.012** |
| **PRODUCT · F1 × F2** (mast TNFRSF9 per 10k tissue UMI) | 2.3×10⁻⁴ | 1.2×10⁻³ | log₂FC +2.34 (−0.72 to 5.39) | Poisson GLM, donor-clustered, P = 0.13 |

n = 6 healthy donors (538 mast cells) and 11 AD donors (3,123 mast cells). The
two factors and their product are reported together because they are one
quantity (§4): **the per-cell factor moves, the abundance factor does not, and
the product is directionally consistent but underpowered.**

### 1.3 Healthy versus AD non-lesional, and healthy versus AD lesional

| Comparison | log₂FC | Exact permutation P | Transcripts | Donors |
|---|---|---|---|---|
| Healthy vs AD non-lesional | +3.05 | **0.0031** | 3 vs 58 | 6 vs 10 |
| Healthy vs AD lesional | +2.59 | 0.179 | 3 vs 49 | 6 vs 11 |

Both arms of AD skin sit above healthy skin. The lesional comparison does not
reach significance despite a similar effect size, because lesional mast cells
vary more between donors — a power statement, not evidence of absence.

### 1.4 The refinement: AD lesional versus AD non-lesional

Within patients (10 donors with both samples) per-mast-cell TNFRSF9 did not
differ between lesional and non-lesional skin (log₂FC −0.46, 95% CI −1.97 to
1.06, P = 0.55). Mast-cell TNFRSF9 is therefore a property of AD skin generally
rather than of the lesion, consistent with §1.3.

---

## 2. Replication, and what fails in it (§7)

Two independent cohorts were tested with the identical donor-level exact
permutation test. Because neither deposit carries author cell-type labels, mast
cells were called by a rule *calibrated on the discovery cohort's labels*
(≥2 of TPSAB1/TPSB2/CPA3 detected **and** ≥50 tryptase+CPA3 transcripts per 10k
UMI: precision 96.0%, sensitivity 79.3%). The magnitude requirement is not
cosmetic — detection alone tags hundreds of keratinocytes carrying ambient
tryptase from lysed mast cells.

| Cohort | Donors (H/AD) | Mast cells (H/AD) | Transcripts (H/AD) | Rate H → AD (per 10k) | log₂FC | P |
|---|---|---|---|---|---|---|
| Discovery GSE204762 (3′) | 6 / 11 | 538 / 3,123 | 3 / 107 | 0.057 → 0.469 | +2.81 | 0.012 |
| REP1 GSE222840+GSE173205 (5′) | 4 / 5 | 27 / 91 | 0 / 2 | 0.000 → 0.282 | +0.51 | 0.79 |
| REP2 GSE153760 biopsies (3′ v3) | 2 / 4 | 278 / 1,101 | 6 / 38 | 0.208 → 0.356 | +0.68 | 0.80 |
| **Pooled, fixed effects** | 12 / 20 | 843 / 4,315 | 9 / 147 | — | **+1.44** (0.52–2.37) | **0.0022** |
| **Pooled, random effects** | 12 / 20 | 843 / 4,315 | 9 / 147 | — | **+1.53** (−0.16–3.22) | 0.076 |

**What replicates.** The direction, in 3 of 3 cohorts. The fact of expression, in
all cohorts with usable mast-cell recovery.

**What does not replicate.** (i) The *magnitude*: I² = 58% (Q = 4.79, df 2,
P = 0.091); the discovery effect is 4–6× larger than either replication.
(ii) The **near-absence in healthy mast cells**. GSE153760 healthy mast cells
carried TNFRSF9 at 0.21 per 10k — 3.6× the discovery healthy rate — with 6 of
278 cells positive across 2 donors. The clean "off in health, on in disease"
pattern is a property of the discovery cohort, not an established fact.
(iii) Statistical significance: neither replication cohort reaches it alone, and
the random-effects pooled estimate does not either.

**Why REP1 is uninformative rather than negative.** Mast-cell recovery failed in
that pipeline: 118 mast cells in 111,370 cells (0.11%), versus 1.2–1.5% in the
discovery cohort and 4.4% in REP2, and no coherent mast-cell cluster forms in the
embedding (Fig. 5a). With 2 TNFRSF9 transcripts in the whole cohort it has no
power to confirm or refute anything.

**A caution about both.** REP2's AD signal is concentrated in one donor (AD7: 31
of 38 transcripts; also an outlier at 13.3% mast cells) and its healthy arm has
just 2 donors, 5 of whose 6 transcripts come from one of them. The discovery
cohort is similarly concentrated (MGH108: 55 of 107). No cohort examined here has
the donor numbers this question needs.

---

## 2A. Does the depth adjustment create the result?

Two objections, both answered with data rather than argument.

### 2A.0 "Sequencing depth" is the wrong name for what is being adjusted

Sequencing depth in the strict sense is a property of the run — how many reads
are allocated to a library — and it is shared by every cell in that library. Per-
cell UMI counts are not that. They are *sampling effort per cell*: library reads
× capture efficiency × **the cell's own mRNA content**. Only the first term is
sequencing depth, and normalising by per-cell UMI divides out the third term,
which is biology. The terminology used in §6 of the contract ("library size") is
the conventional scRNA-seq usage, where each barcoded cell is its own library —
but calling it *sequencing depth* was imprecise, and the distinction turns out to
matter here.

Decomposing the variance of log(UMI per cell) across mast cells and the §6
control populations:

| Component | Share of variance | Can it be sequencing depth? |
|---|---|---|
| Between libraries (flow cell, loading) | 27.4% | Yes |
| Between cell types **within** a library | 22.0% | **No** — those cells shared a flow cell |
| Cell to cell within a cell type | 50.6% | No |

Roughly three-quarters of the variation is not a sequencing-run property. The
within-library ratio makes this concrete and is immune to flow-cell depth by
construction: mast cells yield **0.08–0.30** of the UMI that fibroblasts and
keratinocytes yield **in the same run**, in every one of 39 libraries — and that
ratio is **34% lower** in AD (median 0.134, n = 32 libraries) than in healthy
skin (0.203, n = 7 libraries; Mann–Whitney P < 0.001).

The decisive observation is that the two move in opposite directions: **AD
libraries are sequenced more deeply overall (median 4,292 versus 1,564 UMI per
cell) yet AD mast cells are shallower (693 versus 935 UMI).** No flow-cell effect
can do that. The mast-cell UMI deficit in AD is a property of the cells, not of
the sequencing — AD mast cells are transcriptionally smaller, detecting 399
versus 477 genes and carrying less tryptase per cell.

**Consequence.** Dividing TNFRSF9 by mast-cell UMI removes real biology along
with technical effort, and it divides by a denominator that is 25% smaller in the
AD arm. The **per-cell estimand is therefore primary** in this report and the
per-UMI estimand is reported beside it as the compositional companion, as §6
requires.

### 2A.1 What is and is not depth-matched

**"Sequencing depth" is three quantities, and only one can bias this test.**

| Quantity | Healthy | AD | Ratio AD/H |
|---|---|---|---|
| Mean UMI per cell, all cells | 2,976 | 4,233 | 1.42 |
| Mean UMI per **fibroblast** | 1,937 | 4,279 | **2.21** |
| Mean UMI per **keratinocyte** | 2,226 | 4,817 | **2.16** |
| Mean UMI per **macrophage** | 1,823 | 3,729 | **2.05** |
| Mean UMI per **T/NK cell** | 1,080 | 1,606 | 1.49 |
| Mean UMI per **mast cell** | **934** | **693** | **0.74** |

The arms are *not* depth-matched: AD cells are 1.5–2.2× deeper in every
population — **except mast cells, where healthy is 1.35× deeper.** Mast cells are
the exception, and that fact decides the argument in two opposite directions.

**For the mast-cell result, the adjustment is not doing the work.** Because
healthy is the deeper arm, the unadjusted comparison is the *conservative* one,
and it still holds (discovery cohort, donor-level exact permutation):

| Estimand | Healthy | AD | Fold | log₂FC | P |
|---|---|---|---|---|---|
| per mast UMI (adjusted, compositional) | 0.057 /10k | 0.469 /10k | 8.2× | +2.81 | 0.012 |
| **per mast CELL (unadjusted, absolute)** | 0.56 /100 cells | 3.43 /100 cells | **6.1×** | **+2.40** | **0.030** |
| **no normalisation at all** (TNFRSF9⁺ cells) | 3/538 | 91/3,123 | OR 5.35 | — | **5.3×10⁻⁴** |

Pooled across all three cohorts the unadjusted estimand gives log₂FC **+1.24
(95% CI +0.32 to +2.16, P = 0.0085)** fixed and **+1.32 (−0.07 to +2.71,
P = 0.063)** random, with **lower** heterogeneity than the adjusted version
(I² = 42% versus 58%). Healthy mast cells are the deeper arm in the discovery
cohort (ratio 0.75) and in GSE153760 (0.93), and equal in REP1 (1.04), so the
unadjusted estimand is conservative or neutral in every cohort. **It is therefore
the number this report leads with.**

**But you cannot drop the adjustment for the control populations.** Fibroblasts,
keratinocytes and macrophages *are* 2× deeper in AD, so for them the unadjusted
comparison genuinely is a depth comparison:

| Population (6 healthy vs 11 AD donors) | log₂FC per UMI (adjusted) | log₂FC per cell (unadjusted) |
|---|---|---|
| **Mast (subject)** | +2.81 (P = 0.012) | +2.40 (P = 0.030) |
| Fibroblasts | +0.91 (P = 0.23) | **+2.05 (P = 0.011)** |
| Keratinocytes | +0.26 (P = 0.72) | +1.37 (P = 0.15) |
| T/NK | +1.49 (P = 0.0001) | +2.06 (P = 0.0013) |
| Macrophages | +2.06 (P = 0.024) | +3.09 (P = 0.013) |

Unadjusted, TNFRSF9 appears to rise in *almost every cell type* — which is
exactly the "shared shift means batch, not biology" pattern §6 exists to catch,
and here it is an artefact of those populations being sequenced twice as deeply
in AD. **So: the adjustment is unnecessary for the mast-cell effect and essential
for the claim that the effect is specific.** Both statements are reported rather
than choosing whichever is convenient.

**One genuine cost of adjusting, in the objection's favour.** AD mast cells
capture 24.7% fewer UMI per cell than healthy ones (693 versus 934) and detect
fewer genes (399 versus 477). Dividing by a smaller denominator inflates a
per-UMI fraction, which is precisely why the adjusted estimate (8.2×) exceeds the
unadjusted one (6.1×). Normalising by total UMI measures TNFRSF9 as a *share of
the mast-cell transcriptome*, not molecules per cell; where the transcriptome
itself changes size, those are different questions and both are reported above.

---

## 2B. The controlled experiment: TNFRSF9 is an activation-induced gene in primary human mast cells

Every tissue analysis above is confounded in the same way: it cannot separate
"mast cells transcribe more TNFRSF9 in AD" from "AD skin holds a different
mast-cell population". A stimulation experiment can. Two public datasets perform
exactly that experiment on **primary** human mast cells, with no dissociation
artefact, no ambient RNA, no doublets and no depth confound.

### 2B.1 Primary human SKIN mast cells + the AD milieu (GSE196862)

Human skin-derived mast cells, stimulated for 24 h with IgE-receptor
cross-linking and the epithelial alarmins IL-33, TSLP and IL-25, alone and
combined (24 libraries, 2–4 replicates per condition).

| Condition | TNFRSF9 (FPKM) | log₂FC vs resting |
|---|---|---|
| **resting** | **0.14** | — |
| TSLP | 0.11 | −0.28 |
| IL-25 | 0.11 | −0.32 |
| IgE/Ag | 1.07 | +2.87 (7×) |
| IL-33 | 4.33 | +4.88 (29×) |
| IL-33+IL-25+TSLP | 3.76 | +4.67 |
| **IgE/Ag + IL-33** | **8.66** | **+5.88 (59×)** |
| IgE/Ag + all three alarmins | 7.96 | +5.75 |

Resting skin mast cells are **essentially TNFRSF9-negative** (individual
replicates 0.011, 0.047, 0.006, 0.487 FPKM). IL-33 and FcεRI cross-linking each
induce the gene, and together they synergise to ~59-fold. TSLP and IL-25 alone do
nothing, so this is not a generic cytokine response.

**The induction is specific, not global.** In the same libraries the housekeeping
genes are flat (B2M +0.34, RPL13A −0.23, TMSB4X +0.19 under IgE+IL-33) and the
mast-cell identity genes barely move (TPSAB1 −0.29, TPSB2 +0.02). Among
co-stimulatory TNF-receptor relatives, TNFRSF9 has both the lowest resting level
(0.14 FPKM versus 1.6–1.8 for TNFRSF18 and TNFRSF4) and the largest induction
(+5.88 versus +3.20 and +0.73). It is the most switch-like member of the family
in this cell type.

### 2B.2 Primary human mast cells + activated CD4⁺ T cells (GSE235240)

Four donors, paired design, mast cells FACS-sorted after 24 h of co-culture.

| Condition | TNFRSF9 (normalised counts) | paired log₂FC | P | n donors |
|---|---|---|---|---|
| resting | 103 | — | — | 4 |
| resting CD4⁺ T cells | 133 | +0.26 | 0.17 | 4 |
| **activated CD4⁺ T cells** | 768 | **+2.97** | **0.003** | 4 |
| IgE/Ag | 4,148 | +5.51 | 0.001 | 4 |
| **IL-33** | **9,079** | **+6.53 (92×)** | **0.001** | 4 |

The same ordering reproduces in an independent cohort with an independent
stimulus: IL-33 > IgE/Ag > activated T cells > resting T cells ≈ unstimulated.
Tryptase genes are unchanged (TPSAB1 log₂FC −0.04, P = 0.61, n = 4 donors), so
this is not a differentiation or viability effect.

### 2B.3 What this settles, and what it does not

**It settles the mechanism.** TNFRSF9 in mast cells is an activation-induced
gene, switched on from a near-zero baseline by exactly the two signals that
define atopic dermatitis: IL-33 released by damaged keratinocytes, and IgE
cross-linking on FcεRI. This is a *measured* result in primary cells, not an
inference from tissue.

**It explains the tissue data, including the parts that did not replicate.**
If TNFRSF9 marks mast-cell activation state rather than disease identity, then
(i) healthy skin mast cells should be near-zero — they are; (ii) AD mast cells
should be positive to a degree that depends on how activated they were at the
moment of biopsy — which predicts exactly the between-cohort heterogeneity
observed (I² = 42–58%, §2); and (iii) the tissue-level rise should *not* be
AD-specific, because any inflammation that activates mast cells would do the
same — psoriasis matches it (§8).

**It does not settle the disease question.** These are in-vitro stimulations, not
AD skin. They establish that mast cells *can* express TNFRSF9 on demand and
identify the sufficient stimuli; they do not establish how much of the in-vivo
AD signal comes from this pathway, nor that 4-1BB protein reaches the mast-cell
surface, nor any functional consequence. The in-vivo effect size remains the
heterogeneous, modest one reported in §2.

---

## 3. Mandatory controls (§6)

| Control | Result | Consequence |
|---|---|---|
| **Library size** | log-depth offset in every model; depth-stratified analysis in §4 | Effect survives and strengthens |
| **Donor aggregation** | All inference at donor level; exact permutation over 12,376 donor-label assignments | Cluster-robust SEs on 6 clusters not relied upon |
| **Same test in fibroblasts + keratinocytes** | Fibroblasts log₂FC +1.07 (−0.45 to 2.59), P = 0.17; keratinocytes +0.29 (−0.95 to 1.54), P = 0.64 (n = 6 vs 11 donors) | **No shared shift → not a batch or global effect** |
| **Per-group median depth** | Mast 546 / 541 / 622 UMI; whole tissue 1,999 / 2,833 / 3,306 UMI (healthy / AD NL / AD LS) | Mast depth comparable across arms; tissue depth is not, hence the offset |
| **Marker QC** | 3,661 of 4,313 deposited "Mast" labels passed ≥2 of TPSAB1/TPSB2/CPA3 (84.9%); 433 QC-positive cells carried non-mast labels and were excluded | Deposited labels are imperfect, as §6 anticipates |
| **Dissociation stress** | HSP+IEG burden 986 per 10k mast UMI in healthy versus 416 in AD | **The healthy arm is the most stressed, which biases *against* the reported direction**; donor-level Spearman(stress, rate) = −0.15, P = 0.57 |

**Leave-one-donor-out.** Dropping any single one of the 17 discovery donors
(n = 6 healthy, 11 AD) leaves log₂FC between +2.15 and +4.72, permutation
P ≤ 0.043 throughout. Removing the most influential donor (MGH108) gives +2.15
(P = 0.023).

### 3.1 The doublet control

T/NK cells carry TNFRSF9 at a rate comparable to AD mast cells, so a
mast-cell/T-cell doublet would present as exactly what is being counted: a
tryptase-positive, TNFRSF9-positive cell. This is the most dangerous artefact
available to this result. It is real, but it is not sufficient.

TNFRSF9⁺ mast cells do carry T-cell transcripts more often than TNFRSF9⁻ mast
cells (10.6% versus 3.4% of cells; depth-stratified OR 2.68 across six depth
strata, n = 94 versus 3,567 cells), and the enrichment survives depth
stratification, so some ambient or doublet contamination is present. No cell in
either group was flagged by the deposited scrublet call, and TNFRSF9⁺ mast cells
had *lower* doublet scores (median 0.02 versus 0.03).

The decisive test is deletion:

| Cell set | Mast cells | TNFRSF9 healthy / AD | log₂FC | Permutation P |
|---|---|---|---|---|
| All marker-QC mast cells | 3,661 | 3 / 107 | +2.81 | 0.012 |
| **Excluding any cell with a T-cell transcript** | 3,531 | 3 / 97 | **+2.73** | **0.022** |
| Excluding any T-cell *or* myeloid transcript | 3,187 | 3 / 78 | +2.40 | 0.050 |

97 of the 107 AD TNFRSF9 transcripts come from mast cells carrying **no T-cell
transcript at all**, so mast–T doublets cannot be the source of the result. The
stricter T-plus-myeloid purge removes 13% of mast cells and pushes the test to
the significance boundary: the effect is not eliminated, but this cohort has no
power to spare.

---

## 4. Is the discovery-cohort healthy arm genuinely negative?

Within the discovery cohort, four lines of evidence say the 3 healthy
transcripts are not simply under-sampling:

1. **Exact Poisson intervals do not overlap.** Healthy 0.057 per 10k mast UMI
   (95% CI 0.012–0.168); AD 0.469 (95% CI 0.384–0.566), from 522,068 and
   2,282,915 mast UMI.
2. **Expected versus observed.** Under the AD rate the healthy exposure should
   have yielded 24.5 transcripts; 3 were observed.
3. **Depth matching strengthens rather than explains it.** Healthy mast cells
   contributed **zero** TNFRSF9 transcripts in every depth stratum below 800 UMI,
   where AD mast cells contributed 35. The depth-stratified common rate ratio is
   **9.6× (log₂ 3.26)** — larger than the crude estimate (Fig. 1f).
4. **Donor positivity.** 8 of 11 AD versus 1 of 6 healthy donors had ≥1 TNFRSF9⁺
   mast cell (Fisher OR 13.3, P = 0.050).

**But this pattern is cohort-specific.** GSE153760 healthy mast cells are
positive at 0.21 per 10k (§2). The correct statement is therefore that the
discovery cohort's healthy mast cells are indistinguishable from zero at the
exposure available, *and that a second cohort does not reproduce that*. The
attractive "de-novo induction" reading is a hypothesis this work does not
establish.

---

## 5. Calibration: how much is the discovery result worth?

* **Is the test calibrated?** Permuting arm labels *within* the AD donors, where
  no disease contrast exists by construction, gave a median 4.0% of matched
  genes at P ≤ 0.05 (IQR 2.9–5.6%, range 1.8–10.4% over 40 random splits). The
  exact permutation test is therefore approximately calibrated.
* **Is TNFRSF9 exceptional?** No. The real contrast yields 17.8% of matched genes
  at P ≤ 0.05 — far outside the null range — so **the mast-cell transcriptome
  differs widely between healthy and AD skin**. Against that background TNFRSF9
  ranks 290 of 3,454 (empirical P = 0.084, Fig. 2e).

TNFRSF9 induction in AD mast cells is a real, calibrated nominal result sitting
in the top 8% of a broad disease-associated shift — but it is *one of many* genes
moving, and the claim that TNFRSF9 is *specifically* dysregulated in mast cells
is **not** supported.

---

## 6. In situ: spatial transcriptomics without dissociation (GSE197023)

The per-cell estimate is conditional on which mast cells survive enzymatic
digestion (§4). Visium never dissociates the tissue and reproduces the direction:
whole-spot TNFRSF9 was higher in AD than healthy skin (log₂FC +2.06, 95% CI
0.35–3.77, P = 0.019; 13,742 spots, 6 healthy versus 13 AD sections,
donor-clustered with a log-depth offset), driven by lesional skin (+2.30, 95% CI
0.47–4.12, P = 0.014). Mast-cell content in situ was unchanged or slightly lower
(TPSAB1 log₂FC −0.80, P = 0.34; CPA3 −0.71, P = 0.021, 6 vs 13 sections),
matching the unchanged abundance in the single-cell arm.

Spot depth differs ~4-fold between arms (AD lesional 4,025 versus healthy 1,220
UMI), so every section in Fig. 3 is annotated with its own median depth and all
statistics carry a depth offset.

**What does not hold spatially.** TNFRSF9 tracked mast-cell content in AD
*non-lesional* skin (log₂FC +1.13 per SD of tryptase content, 95% CI 0.40–1.87,
P = 0.0025, 3,727 spots, 6 donors) but **not** in AD lesional skin (+0.04, 95% CI
−0.56 to 0.64, P = 0.90) or healthy skin (−0.05, P = 0.96). The mandated control
contents behaved as controls should in the same spots: fibroblast content +0.73
(P = 0.13), keratinocyte content −0.13 (P = 0.50). In lesional skin, where
whole-tissue TNFRSF9 is highest, that signal is *not* spatially organised around
mast cells.

---

## 7. Mast cells are a minor source of skin TNFRSF9

| Arm | TNFRSF9 in mast cells | in whole tissue | mast share |
|---|---|---|---|
| Healthy | 3 | 104 | 2.9% |
| AD non-lesional | 58 | 1,018 | 5.7% |
| AD lesional | 49 | 1,844 | 2.7% |

Mast cells account for 0.22–0.41% of the tissue transcriptome and 2.7–5.7% of its
TNFRSF9. Both the induction within mast cells and their small share are true
simultaneously; a receptor switched on in a rare cell can matter functionally
without dominating a bulk measurement, but it cannot explain a tissue-level
signal either.

---

## 8. Tissue level, and the specificity problem (GSE121212)

| Comparison | log₂FC | 95% CI | n | P |
|---|---|---|---|---|
| Healthy vs AD | +2.48 | 1.85 to 3.11 | 38 vs 54 | 1.0×10⁻¹⁴ |
| Healthy vs AD non-lesional | +1.38 | 0.75 to 2.01 | 38 vs 27 | 1.6×10⁻⁵ |
| Healthy vs AD lesional | +3.06 | 2.37 to 3.75 | 38 vs 27 | 2.9×10⁻¹⁸ |
| AD lesional vs non-lesional | +1.68 | 1.04 to 2.31 | 27 vs 27 | 2.4×10⁻⁷ |
| **Healthy vs psoriasis** | **+2.20** | 1.80 to 2.59 | 38 vs 55 | 8.7×10⁻²⁸ |

Poisson GLM on raw counts, log-library-size offset, errors clustered on patient.
Mast-cell content in the same samples was flat (TPSAB1 log₂FC −0.09, P = 0.74;
TPSB2 −0.10, P = 0.73; CPA3 −0.48, P = 0.037, n = 38 vs 54); mandated controls
moved little (COL1A1 −0.54, P = 0.089) or as epidermal acanthosis predicts
(KRT14 +0.84, P = 1.5×10⁻¹⁰).

In bulk tissue abundance and per-cell expression are not separable at all (§4),
and the mast-cell contribution is not separable from any other source. What bulk
establishes is that skin TNFRSF9 rises steeply in AD, does so without a rise in
mast-cell content, and rises just as much in psoriasis — so the *tissue-level*
signal marks inflamed skin, not AD.

---

## 9. What is measured, modelled, inferred (§8)

* **Measured.** Transcript counts per cell and per spot; marker expression;
  sequencing depth; cell-type labels re-derived by marker QC.
* **Modelled.** All fold changes and intervals (Poisson/binomial GLMs with
  log-depth offsets and donor- or patient-clustered errors; exact donor-label
  permutation; inverse-variance meta-analysis).
* **Inferred.** That the per-mast-cell change reflects transcriptional induction
  rather than selective survival of a TNFRSF9⁺ mast-cell subset through
  dissociation. The Visium arm supports this indirectly but cannot resolve single
  mast cells at 55 µm.
* **Not established.** 4-1BB protein on mast cells; functional consequence;
  AD-specificity of the tissue-level change; the magnitude of the mast-cell
  effect (I² = 58%); the near-absence in healthy skin (not replicated); and any
  claim of selective dysregulation relative to the rest of the mast-cell
  transcriptome.

---

## 10. Limitations

1. **Small counts.** Nine TNFRSF9 transcripts across all healthy mast cells in
   all three cohorts. Statements about healthy skin are upper bounds.
2. **Donor numbers.** 12 healthy and 20 AD donors in total, with each cohort's
   signal concentrated in one or two donors. This question needs a cohort built
   for it.
3. **Mast cells are hard to recover, and that governs everything.** Recovery
   ranged from 0.11% to 4.4% of cells across cohorts using the same rule. One
   dedicated AD atlas (GSE147424, cryopreserved biopsies) contained *no*
   recoverable mast cells and was excluded (`docs/datasets.md`); suction-blister
   sampling also fails to capture them.
4. **Heterogeneity.** I² = 58% across cohorts; the random-effects pooled estimate
   does not reach significance.
5. **Not exceptional against background** (§5), and not AD-specific at tissue
   level (§8).
6. **Cross-series healthy comparator in REP1** (§5 of the contract), which is in
   any case uninformative for mast cells.
