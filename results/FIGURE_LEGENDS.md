# Figure legends

Conventions used throughout. Error bars show the standard error of the mean.
Intervals given in the report text are 95% confidence intervals and are labelled
as such. Mast cells were identified by marker quality control, requiring at least
two of TPSAB1, TPSB2 and CPA3 to be detected in a cell. All group comparisons
were made at the level of donors or samples, never individual cells. Statistical
tests are two-sided.

---

## Figure 1 | Cutaneous mast cells and TNFRSF9 in the discovery cohort

**a**, Uniform manifold approximation and projection (UMAP) of 280,518 cells from
39 whole-skin biopsies taken from 17 adults (6 healthy donors and 11 patients
with atopic dermatitis, each sampled at lesional and non-lesional sites;
GSE204762). Cells are coloured by the cell-type annotation supplied with the
dataset. Mast cells are shown in purple.

**b**, The same embedding coloured by the combined counts of the mast-cell
granule genes TPSAB1, TPSB2 and CPA3. The mast-cell island is the only region of
the atlas in which these transcripts are abundant.

**c**, The 3,661 cells that passed mast-cell marker quality control, coloured by
tissue group.

**d**, Mast cells from each tissue group, with TNFRSF9-positive cells shown as
filled circles. Values beneath each panel give the number of TNFRSF9-positive
mast cells, the total number of mast cells and the corresponding percentage.
Three of 538 mast cells were positive in healthy skin (0.6%), compared with 52 of
1,564 in non-lesional atopic dermatitis (3.3%) and 39 of 1,559 in lesional skin
(2.5%).

---

## Figure 2 | Mast-cell numbers and TNFRSF9-positive mast-cell numbers

Quantification of the cells shown in Fig. 1d. Each point is one sample, defined
as one donor at one tissue site (6 healthy, 10 non-lesional and 11 lesional
samples). Horizontal bars show group medians. P values above the brackets are
from Mann-Whitney U tests against the healthy group.

**a**, Mast cells recovered per sample. Recovery was similar across groups
(medians of 99.5 cells in healthy skin, 131 in non-lesional and 151 in lesional
atopic dermatitis; n = 6, 10 and 11 samples; P = 0.37 and P = 0.26).

**b**, Mast cells as a percentage of all cells recovered from the same sample.
Abundance did not differ between groups (median 1.47% in healthy skin, 1.30% in
non-lesional and 1.21% in lesional atopic dermatitis; P = 0.79 and P = 1.00).

**c**, TNFRSF9-positive mast cells per sample.

**d**, TNFRSF9-positive mast cells as a percentage of the mast cells in the same
sample. Pooling cells across samples, 3 of 538 mast cells were positive in
healthy skin and 91 of 3,123 in atopic dermatitis (odds ratio 5.35, Fisher exact
P = 5.3 x 10^-4). Five of the six healthy donors contributed no
TNFRSF9-positive mast cell.

---

## Figure 3 | TNFRSF9 in intact skin measured by spatial transcriptomics

Visium spatial transcriptomics of 19 skin sections from 6 healthy donors and 7
patients with atopic dermatitis (GSE197023), comprising 13,742 spots that passed
a minimum depth filter of 200 unique molecular identifiers (UMIs). Sections were
not dissociated, so the measurement does not depend on enzymatic release of mast
cells. Marker size is scaled to the true spot pitch of each section, which
differs between samples because the underlying scans differ in resolution.

**Top row**, Six representative sections, two per group, with spots coloured by
the combined TPSAB1, TPSB2 and CPA3 content per 10,000 UMIs and overlaid on the
haematoxylin and eosin image.

**Second row**, The same sections with TNFRSF9-positive spots highlighted. The
number of positive spots and the median sequencing depth are given above each
section. Depth per spot differs about fourfold between groups, so it is stated
for every section and adjusted for in all models.

**g**, TNFRSF9 per spot in each comparison, expressed as a log2 fold change
against healthy skin. Estimates come from a Poisson generalised linear model with
log UMIs per spot as an offset and standard errors clustered on donor. The
fibroblast and keratinocyte genes required as controls were carried through the
same test and showed no comparable rise.

**h**, TNFRSF9 per 10,000 spot UMIs for every section in the cohort.

---

## Figure 4 | Whole-skin TNFRSF9 in a large tissue cohort with a disease control

Bulk RNA sequencing of 147 whole-skin biopsies (GSE121212): 38 healthy, 27
non-lesional and 27 lesional atopic dermatitis, and 27 non-lesional and 28
lesional psoriasis. Psoriasis is included as a disease control, to ask whether
any change is specific to atopic dermatitis. Each point is one biopsy and
horizontal bars show medians.

**a**, Whole-skin TNFRSF9 in counts per million. Median values rise from 0.36 in
healthy skin to 0.56 in non-lesional and 1.52 in lesional atopic dermatitis, and
reach 2.46 in lesional psoriasis.

**b**, Combined TPSAB1, TPSB2 and CPA3 content in the same samples, used here as
the tissue-level proxy for mast-cell abundance. Content did not rise with
TNFRSF9. Within groups, whole-skin TNFRSF9 showed no correlation with mast-cell
content (Spearman rho between -0.13 and +0.26, all P > 0.18).

**c**, Log2 fold changes against healthy skin for TNFRSF9 and for TPSAB1, shown
for atopic dermatitis and psoriasis. Estimates come from a Poisson generalised
linear model on raw counts with log library size as an offset and standard errors
clustered on patient. TNFRSF9 rose to a similar degree in both diseases while
mast-cell content did not change, indicating that the tissue-level signal marks
inflamed skin rather than atopic dermatitis specifically.

---

## Figure 5 | Replication in an independent cohort

**a**, UMAP of 31,266 cells from six whole-skin biopsies in GSE153760 (4 patients
with atopic dermatitis, 2 healthy donors). The 1,379 cells meeting mast-cell
criteria are shown in purple. Suction-blister samples from the same study were
excluded because they do not recover mast cells.

**b**, Mast cells from this cohort with TNFRSF9-positive cells highlighted and
coloured by group.

**c**, TNFRSF9-positive mast cells as a percentage of mast cells, shown for the
discovery cohort and for GSE153760 on the same axis. Each point is one sample and
horizontal bars show medians. Odds ratios and Fisher exact P values were computed
on pooled cell counts within each cohort. The increase seen in the discovery
cohort (0.56% to 2.91%, odds ratio 5.35, P = 0.001) was not reproduced in
GSE153760 (2.16% to 2.45%, odds ratio 1.14, P = 1.00). The two cohorts agree on
the atopic dermatitis arms and differ fourfold on the healthy arms.

**d**, Mast-cell recovery for the three single-cell cohorts screened, expressed as
the percentage of all recovered cells. The dashed line marks the inclusion
threshold of 0.5%. GSE222840 combined with GSE173205 recovered 0.106% of cells as
mast cells, giving 27 mast cells in its healthy arm and 2 TNFRSF9 transcripts
across 111,370 cells, and was therefore excluded before any TNFRSF9 comparison
was made.

---

## Figure 6 | TNFRSF9 is induced in primary human mast cells by IL-33 and IgE receptor cross-linking

**a**, TNFRSF9 expression in primary human skin-derived mast cells stimulated for
24 h (GSE196862, 24 libraries, 2 to 4 replicates per condition). Resting cells
express almost no TNFRSF9 (0.14 FPKM). IL-33 raised expression 29-fold and IgE
receptor cross-linking 7-fold; together they raised it 59-fold, to 8.66 FPKM.
TSLP and IL-25 given alone had no effect.

**b**, Log2 fold changes against resting cells for TNFRSF9, related TNF receptor
family members, activation genes, mast-cell identity genes and housekeeping
genes. Housekeeping and granule genes were unchanged, so the induction is
specific rather than a global shift in the transcriptome.

**c**, TNFRSF9 in primary human mast cells co-cultured for 24 h with autologous
resting or activated CD4+ T cells, or stimulated with IL-33 or IgE and antigen
(GSE235240). Mast cells were sorted after culture. Four donors were tested and
each condition is paired within donor. P values are from paired t-tests on log2
normalised counts. Note the logarithmic axis.

---

## Supplementary Figure 1 | Controls and calibration for the discovery cohort

**a**, TNFRSF9 in atopic dermatitis compared with healthy skin, tested separately
in mast cells and in the control populations required by the project protocol.
Estimates are log2 fold changes from Poisson generalised linear models using each
population's own total UMIs as the offset, with standard errors clustered on
donor. Fibroblasts and keratinocytes showed no shift, indicating that the
mast-cell result is not a batch or global effect.

**b**, Distribution of log2 fold changes between healthy skin and atopic
dermatitis for 3,454 genes matched to TNFRSF9 for abundance in mast cells.
TNFRSF9 is marked. It ranks 290th of 3,454, an empirical P value of 0.084, so it
sits in the upper tail of a broad disease-associated shift rather than standing
apart from it.

**c**, Dissociation-stress burden, measured as combined heat-shock and
immediate-early gene counts per 10,000 mast-cell UMIs, for each sample. The
healthy group carried the highest stress burden. Since stressed cells return
lower counts for other transcripts, this works against rather than towards the
reported increase in atopic dermatitis.
