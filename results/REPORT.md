# TNFRSF9 (CD137/4-1BB) is absent from healthy cutaneous mast cells and induced in atopic dermatitis, while mast-cell abundance is unchanged

**Subject of this study: cutaneous mast cells.** Fibroblasts, keratinocytes,
T/NK cells and macrophages appear throughout as control populations only, to
establish that a mast-cell result is specific and not a batch or global effect
(CLAUDE.md §2).

---

## Summary

**TNFRSF9 is expressed by human cutaneous mast cells, and the expression is a
disease-associated event: it is undetectable in healthy mast cells and present
in atopic dermatitis (AD).** Across 3,661 marker-QC mast cells from 17 donors in
the discovery cohort (GSE204762), TNFRSF9 was carried by 3 of 538 healthy mast
cells (0.56%; 3 transcripts in 522,068 mast UMI from 6 donors, 5 of whom
contributed none) versus 91 of 3,123 AD mast cells (2.91%; 107 transcripts in
2,282,915 mast UMI from 11 donors, 8 of whom contributed at least one). On the
donor-aggregated rate this is a log₂ fold change of **+2.81 (8-fold), exact
permutation P = 0.012** over all 12,376 donor-label assignments.

**The change is in per-cell expression, not in cell number.** Mast-cell
abundance was 1.47% of cells in healthy skin and 1.41% in AD (log₂OR 0.66, 95%
CI −0.78 to 2.10, P = 0.37, 6 vs 11 donors). Because mast-cell-derived TNFRSF9
is the product of abundance and per-cell expression (§4), the product — mast
TNFRSF9 per unit tissue transcriptome — rises about 5-fold but does **not**
reach significance (log₂FC +2.34, 95% CI −0.72 to 5.39, P = 0.13), since it
inherits the variance of both factors.

**Three findings temper this.** (i) The increase is *not exceptional*: among
3,454 genes matched to TNFRSF9 for abundance in mast cells, 289 shift more
strongly between healthy and AD, placing TNFRSF9 at empirical P = 0.084 against
the genome-wide background. The mast-cell transcriptome differs broadly in AD
(17.8% of matched genes at P ≤ 0.05, against a true-null median of 4.0%).
(ii) Mast cells are a **minor source** of skin TNFRSF9, contributing 2.7–5.7% of
all TNFRSF9 transcripts in the tissue. (iii) At tissue level the TNFRSF9 rise is
**not AD-specific**: whole-skin TNFRSF9 increases just as much in psoriasis
(log₂FC +2.20, 95% CI 1.80–2.59, n = 38 healthy vs 55 psoriasis).

---

## 1. The primary comparison: healthy versus AD (§3)

### 1.1 Is TNFRSF9 expressed in cutaneous mast cells at all?

Yes, but at low abundance. Pooling all arms of the discovery cohort, 94 of 3,661
mast cells (2.57%) carried at least one TNFRSF9 transcript, a rate of 0.392
transcripts per 10,000 mast UMI. This is a genuinely low-expression gene in this
cell type, and the entire analysis is therefore governed by counting statistics
rather than by fold-change estimation.

Mast cells are also the **shallowest-sequenced population in the atlas** — median
622 UMI per cell in AD lesional, 541 in AD non-lesional and 546 in healthy skin,
against 3,157–4,061 UMI for fibroblasts and keratinocytes (Fig. 1e). Since
TNFRSF9 detection is depth-dominated, every comparison below carries a
log-library-size offset and is additionally checked against depth-matched strata
(§6).

### 1.2 Healthy versus AD, pooled (primary)

| Quantity (§4) | Healthy | AD | Effect | Test |
|---|---|---|---|---|
| **F1 · abundance** (mast % of cells) | 1.47% | 1.41% | log₂OR +0.66 (−0.78 to 2.10) | binomial GLM, donor-clustered, P = 0.37 |
| **F2 · per-cell expression** (TNFRSF9 per 10k mast UMI) | 0.057 | 0.469 | log₂FC **+2.81** | exact donor permutation, **P = 0.012** |
| **PRODUCT · F1 × F2** (mast TNFRSF9 per 10k tissue UMI) | 2.3×10⁻⁴ | 1.2×10⁻³ | log₂FC +2.34 (−0.72 to 5.39) | Poisson GLM, donor-clustered, P = 0.13 |

n = 6 healthy donors (538 mast cells) and 11 AD donors (3,123 mast cells).

The two factors and their product are reported together because they are one
quantity (§4). The honest reading is that **the per-cell factor moves and the
abundance factor does not**, and that the product is directionally consistent but
underpowered.

### 1.3 Healthy versus AD non-lesional, and healthy versus AD lesional

| Comparison | log₂FC | Exact permutation P | Counts | Donors |
|---|---|---|---|---|
| Healthy vs AD non-lesional | +3.05 | **0.0031** | 3 vs 58 | 6 vs 10 |
| Healthy vs AD lesional | +2.59 | 0.179 | 3 vs 49 | 6 vs 11 |

Both arms of AD skin sit above healthy skin. The lesional comparison does not
reach significance despite a similar effect size because lesional mast cells are
more variable between donors; this is a power statement, not evidence of absence.

### 1.4 The refinement: AD lesional versus AD non-lesional

Within patients (10 donors with both samples), per-mast-cell TNFRSF9 did not
differ between lesional and non-lesional skin (log₂FC −0.46, 95% CI −1.97 to
1.06, P = 0.55). TNFRSF9 induction in mast cells is therefore a property of AD
skin generally rather than of the lesion, which is consistent with the
non-lesional comparison above being the stronger of the two.

---

## 2. Is healthy skin genuinely negative? (de-novo induction)

The healthy arm contributes 3 TNFRSF9 transcripts in total, so the question of
whether healthy mast cells are simply *under-sampled* has to be settled directly.
Four lines of evidence say they are not:

1. **Exact Poisson intervals do not overlap.** Healthy: 0.057 per 10k mast UMI
   (95% CI 0.012–0.168). AD: 0.469 (95% CI 0.384–0.566), from 522,068 and
   2,282,915 mast UMI respectively.
2. **Expected-versus-observed.** Under the AD rate, the healthy exposure should
   have yielded 24.5 transcripts; 3 were observed.
3. **Depth matching removes the concern rather than creating it.** Stratifying
   mast cells by sequencing depth, healthy mast cells contributed **zero**
   TNFRSF9 transcripts in every stratum below 800 UMI, where AD mast cells
   contributed 35. The depth-stratified common rate ratio is **9.6×
   (log₂ 3.26)** — *larger* than the crude estimate, so depth suppresses rather
   than manufactures the difference (Fig. 1f).
4. **Donor-level positivity.** 8 of 11 AD donors versus 1 of 6 healthy donors had
   at least one TNFRSF9⁺ mast cell (Fisher OR 13.3, P = 0.050).

The pattern is therefore better described as **induction from a near-zero
baseline** than as a fold-change between two expressing states.

---

## 3. Mandatory controls (§6)

| Control | Result | Consequence |
|---|---|---|
| **Library size** | log-depth offset in every model; depth-stratified analysis in §2 | Effect survives and strengthens |
| **Donor aggregation** | All inference at donor level; exact permutation over 12,376 donor-label assignments | Cluster-robust SEs on 6 clusters not relied upon |
| **Same test in fibroblasts + keratinocytes** | Fibroblasts log₂FC +1.07 (−0.45 to 2.59), P = 0.17; keratinocytes +0.29 (−0.95 to 1.54), P = 0.64 (n = 6 vs 11 donors) | **No shared shift → not a batch or global effect** |
| **Per-group median depth** | Mast: 546 (healthy), 541 (AD NL), 622 (AD LS) UMI; whole tissue 1,999 / 2,833 / 3,306 UMI | Mast depth comparable across arms; tissue depth is not, hence the offset |
| **Marker QC** | 3,661 of 4,313 deposited "Mast" labels passed ≥2 of TPSAB1/TPSB2/CPA3 (84.9%); a further 433 QC-positive cells carried non-mast labels and were excluded from the strict set | Deposited labels are imperfect, as §6 anticipates |
| **Dissociation stress** | HSP+IEG burden 986 per 10k mast UMI in healthy versus 416 in AD | **The healthy arm is the most stressed, which biases *against* the reported direction**; donor-level Spearman(stress, TNFRSF9 rate) = −0.15, P = 0.57 |

**Leave-one-donor-out.** Dropping any single one of the 17 donors (n = 6 healthy,
11 AD) leaves log₂FC between +2.15 and +4.72 with permutation P ≤ 0.043 throughout. The most influential donor
(MGH108, 55 of 107 AD transcripts) reduces the effect to +2.15 (P = 0.023) when
removed; the result does not depend on it.

---

## 4. Calibration: how much is this result worth?

A nominal P value is not enough when 3,454 comparably expressed genes are
available to be tested. Two questions were separated:

* **Is the test calibrated?** Permuting arm labels *within* the AD donors, where
  no disease contrast exists by construction, gave a median of 4.0% of matched
  genes at P ≤ 0.05 (IQR 2.9–5.6%, range 1.8–10.4% over 40 random splits). The
  exact permutation test is therefore approximately calibrated.
* **Is TNFRSF9 exceptional?** No. The real healthy-versus-AD contrast yields
  17.8% of matched genes at P ≤ 0.05 — far outside the null range — so **the
  mast-cell transcriptome differs widely between healthy and AD skin**. Against
  that background TNFRSF9 ranks 290 of 3,454 (empirical P = 0.084, Fig. 2e).

**Interpretation.** TNFRSF9 induction in AD mast cells is a real, calibrated
nominal result (P = 0.012) that sits in the top 8% of the disease-associated
shift, but it is *one of many* genes moving, and the claim that TNFRSF9 is
*specifically* or *uniquely* dysregulated in mast cells is **not** supported.

---

## 5. In situ: spatial transcriptomics without dissociation (GSE197023)

The single-cell per-cell estimate is conditional on which mast cells survive
enzymatic digestion (§4). Visium never dissociates the tissue, and reproduces the
direction: whole-spot TNFRSF9 was higher in AD than healthy skin (log₂FC +2.06,
95% CI 0.35–3.77, P = 0.019; 13,742 spots, 6 healthy vs 13 AD sections,
donor-clustered with a log-depth offset), driven by lesional skin (+2.30, 95% CI
0.47–4.12, P = 0.014). Mast-cell content in situ was unchanged or slightly lower
(TPSAB1 log₂FC −0.80, P = 0.34; CPA3 −0.71, P = 0.021), matching the unchanged
mast abundance in the single-cell arm.

Spot depth differs ~4-fold between arms (AD lesional 4,025 versus healthy 1,220
UMI), so every section in Fig. 3 is annotated with its own median depth and all
statistics carry a depth offset.

**What does not replicate spatially.** TNFRSF9 tracked mast-cell content in AD
*non-lesional* skin (log₂FC +1.13 per SD of tryptase content, 95% CI 0.40–1.87,
P = 0.0025, 3,727 spots, 6 donors) but **not** in AD lesional skin (+0.04, 95% CI
−0.56 to 0.64, P = 0.90) or healthy skin (−0.05, P = 0.96). The mandated control
contents behaved as controls should: fibroblast content +0.73 (P = 0.13) and
keratinocyte content −0.13 (P = 0.50) in the same non-lesional spots. So in
lesional skin, where whole-tissue TNFRSF9 is highest, that signal is *not*
spatially organised around mast cells — consistent with mast cells being a minor
contributor (§6 below).

---

## 6. Mast cells are a minor source of skin TNFRSF9

| Arm | TNFRSF9 transcripts in mast cells | in whole tissue | mast share |
|---|---|---|---|
| Healthy | 3 | 104 | 2.9% |
| AD non-lesional | 58 | 1,018 | 5.7% |
| AD lesional | 49 | 1,844 | 2.7% |

Mast cells account for 0.22–0.41% of the tissue transcriptome and 2.7–5.7% of its
TNFRSF9. Both the induction within mast cells and their small share of the total
are true simultaneously, and the second does not cancel the first: a receptor
switched on in a rare cell can matter functionally without dominating a bulk
measurement.

---

## 7. Tissue level, and the specificity problem (GSE121212)

Whole-skin RNA-seq in an independent 147-sample cohort confirms a large
tissue-level increase and, importantly, shows it is **not specific to AD**.

| Comparison | log₂FC | 95% CI | n | P |
|---|---|---|---|---|
| Healthy vs AD | +2.48 | 1.85 to 3.11 | 38 vs 54 | 1.0×10⁻¹⁴ |
| Healthy vs AD non-lesional | +1.38 | 0.75 to 2.01 | 38 vs 27 | 1.6×10⁻⁵ |
| Healthy vs AD lesional | +3.06 | 2.37 to 3.75 | 38 vs 27 | 2.9×10⁻¹⁸ |
| AD lesional vs non-lesional | +1.68 | 1.04 to 2.31 | 27 vs 27 | 2.4×10⁻⁷ |
| **Healthy vs psoriasis** | **+2.20** | 1.80 to 2.59 | 38 vs 55 | 8.7×10⁻²⁸ |

Poisson GLM on raw counts, log-library-size offset, standard errors clustered on
patient. Mast-cell content in the same samples was flat (TPSAB1 log₂FC −0.09,
P = 0.74; TPSB2 −0.10, P = 0.73; CPA3 −0.48, P = 0.037, n = 38 vs 54), and the
mandated controls moved little (COL1A1 −0.54, P = 0.089) or in a direction
attributable to epidermal acanthosis (KRT14 +0.84, P = 1.5×10⁻¹⁰).

In bulk tissue, abundance and per-cell expression are not separable at all (§4),
and neither is the mast-cell contribution separable from any other source. What
bulk establishes is that (a) skin TNFRSF9 rises steeply in AD, (b) it does so
without a rise in mast-cell content, and (c) psoriasis does the same — so the
*tissue-level* signal is a marker of inflamed skin rather than of AD.

---

## 8. What is measured, what is modelled, what is inferred (§8)

* **Measured.** Transcript counts per cell and per spot; mast-cell marker
  expression; sequencing depth; cell-type labels re-derived by marker QC.
* **Modelled.** All fold changes and intervals (Poisson/binomial GLMs with
  log-depth offsets and donor- or patient-clustered errors; exact donor-label
  permutation for the discovery contrast).
* **Inferred.** That the per-mast-cell change reflects transcriptional induction
  rather than selective survival of a TNFRSF9⁺ mast-cell subset through
  dissociation. The Visium arm supports this indirectly (no dissociation, same
  direction) but cannot resolve single mast cells at 55 µm.
* **Not established.** Protein-level 4-1BB on mast cells; functional consequence;
  AD-specificity of the tissue-level change (psoriasis matches it); and any
  claim that TNFRSF9 is *selectively* dysregulated relative to the rest of the
  mast-cell transcriptome.

---

## 9. Limitations

1. **Small counts.** The healthy arm rests on 3 transcripts. Every conclusion
   about healthy skin is a statement about an upper bound (95% CI ≤ 0.168 per
   10k mast UMI), not about a measured non-zero level.
2. **Unbalanced design.** 6 healthy versus 11 AD donors, and 538 versus 3,123
   mast cells. Healthy skin is the limiting arm in every cohort examined.
3. **Mast cells are hard to recover.** One dedicated AD atlas (GSE147424,
   cryopreserved biopsies) contained *no* recoverable mast cells at all and had
   to be excluded (see `docs/datasets.md`); suction-blister sampling likewise
   fails to capture them. Any mast-cell result in skin is conditional on the
   dissociation protocol.
4. **Not exceptional against background.** §4 above.
5. **Cross-series healthy comparator in replication 1.** The healthy arm comes
   from a companion series of the same laboratory and protocol, not the same
   study (§5).
