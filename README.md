# TNFRSF9 (CD137/4-1BB) in cutaneous mast cells in atopic dermatitis

Re-analysis of six public human transcriptomic cohorts — single-cell, spatial,
bulk and in-vitro — to answer one question:

> **Is TNFRSF9 expressed in cutaneous MAST CELLS in atopic dermatitis, and does
> that expression change between healthy and AD skin?**

The project runs under a binding contract, [`CLAUDE.md`](CLAUDE.md), which fixes
the subject (mast cells), the comparison order, the mandatory controls and the
reporting standard. `analysis/guardrails.py` enforces it and exits non-zero on
violation.

---

## The answer

**TNFRSF9 is an activation-induced gene in mast cells, not a disease marker.**

* **Expressed?** Yes, at low level. 2.4–2.6% of AD mast cells carry a transcript;
  reproducible across cohorts.
* **Changed in AD?** Higher in 3 of 3 single-cell cohorts, but heterogeneous:
  pooled random-effects log₂FC **+1.32 (95% CI −0.07 to +2.71), P = 0.063,
  I² = 42%**. Abundance does not change — the whole effect is per-cell
  expression.
* **Why heterogeneous?** Because it tracks activation state, not disease. In
  **primary human skin mast cells**, resting cells are essentially negative
  (0.14 FPKM) and the gene is induced **29× by IL-33, 7× by IgE cross-linking and
  59× by both** — the two signals that define AD. TSLP and IL-25 alone do
  nothing.
* **Caveats that matter.** Mast cells supply only 2.7–5.7% of skin TNFRSF9;
  whole-tissue TNFRSF9 rises as much in psoriasis as in AD; and among
  abundance-matched genes TNFRSF9 sits at empirical P = 0.084, so it is not
  *selectively* dysregulated.

Full write-up: **[`results/REPORT.md`](results/REPORT.md)**.
Protein-level evidence: **[`results/PROTEIN_EVIDENCE.md`](results/PROTEIN_EVIDENCE.md)**.

---

## Data

Selection rule (CLAUDE.md §5): dedicated AD research first, healthy and AD arms
from the same study or laboratory pipeline. Full inventory including exclusions:
[`docs/datasets.md`](docs/datasets.md).

| Role | Accession | Design |
|---|---|---|
| Discovery scRNA-seq | GSE204762 | 11 AD (paired NL/LS) + 6 healthy donors, 280,518 cells |
| Replication 1 | GSE222840 + GSE173205 | 5 AD + 4 healthy, same lab, 5′ chemistry |
| Replication 2 | GSE153760 | 4 AD + 2 healthy biopsies, 3′ v3 |
| Spatial | GSE197023 | Visium: 7 AD LS, 6 AD NL, 6 healthy sections |
| Bulk | GSE121212 | 38 healthy, 54 AD, 55 psoriasis (disease control) |
| In vitro (mechanism) | GSE196862 | Primary human **skin** mast cells ± IgE/IL-33/TSLP/IL-25 |
| In vitro (mechanism) | GSE235240 | Primary mast cells ± activated CD4⁺ T cells, 4 paired donors |

**Excluded and documented:** GSE147424 contains *no recoverable mast cells* —
tryptase/CPA3 genes are absent from all 17 deposited matrices. Reported as a
negative result about cryopreserved-biopsy dissociation, not dropped silently.

---

## Reproducing

Raw data are downloaded from GEO by the scripts; nothing is vendored.
Requires ~25 GB of free disk (the discovery archive is streamed sample by sample
and deleted as it goes).

```bash
python3 -m venv .venv
.venv/bin/pip install numpy pandas scipy statsmodels matplotlib h5py \
                     anndata scanpy pyarrow requests igraph leidenalg harmonypy
apt-get install -y fonts-urw-base35        # Nimbus Sans, for the figure standard
```

Run in order:

```bash
.venv/bin/python analysis/s01_extract_gse204762.py    # stream, reduce, discard
.venv/bin/python analysis/s01b_extract_umap.py        # authors' integrated embedding
.venv/bin/python analysis/s02_mast_sc_core.py         # cohort assembly + depth audit
.venv/bin/python analysis/s03_mast_tnfrsf9_discovery.py   # §4 factors and product
.venv/bin/python analysis/s04_robustness.py           # exact permutation, LOO, depth-matched
.venv/bin/python analysis/s05_calibration.py          # true-null calibration
.venv/bin/python analysis/s06_doublet_control.py      # mast-T doublet test
.venv/bin/python analysis/s07_unadjusted.py           # adjusted vs unadjusted estimands
.venv/bin/python analysis/s08_depth_decomposition.py  # flow cell vs cell-intrinsic
.venv/bin/python analysis/s10_bulk_gse121212.py       # bulk + psoriasis control
.venv/bin/python analysis/s20_spatial_gse197023.py    # Visium extraction
.venv/bin/python analysis/s21_spatial_stats.py        # in-situ co-localisation
.venv/bin/python analysis/s40_replication.py          # replication cohorts
.venv/bin/python analysis/s41_replication_test.py     # calibrated mast rule + test
.venv/bin/python analysis/s42_meta.py                 # meta-analysis, both estimands
.venv/bin/python analysis/s50_invitro_stimulation.py  # the mechanism
.venv/bin/python analysis/s3*.py                      # figures 1-6, fig S1
python3 analysis/guardrails.py                        # contract audit
```

---

## Repository layout

```
CLAUDE.md                     binding project contract + current standing answers
analysis/
  guardrails.py               executable contract enforcement (import API + repo audit)
  palette.py                  Nimbus Sans, muted manuscript palette
  genes.py                    marker panels: mast QC, controls, stress, target family
  s0*.py                      discovery cohort: extraction, core analysis, rigour checks
  s1*.py                      bulk arm
  s2*.py                      spatial arm
  s3*.py                      figures 1-6 and supplementary figure 1
  s4*.py                      replication + meta-analysis
  s5*.py                      in-vitro mechanism
docs/datasets.md              every dataset screened, used or excluded, with reasons
results/
  REPORT.md                   the manuscript
  PROTEIN_EVIDENCE.md         protein-level evidence assessment
  figures/                    fig1-6, PDF + PNG
  tables/                     every number in the report, as CSV
data/                         downloaded and processed data (git-ignored)
```

---

## Methodological notes worth knowing before reading the numbers

* **Inference is donor-level and exact.** With 6 healthy donors, cluster-robust
  standard errors are unreliable, so the primary test enumerates all 12,376
  donor-label assignments rather than assuming a distribution.
* **"Depth" is decomposed, not assumed.** Only 27.4% of the variance in
  log(UMI per cell) is between libraries; 22.0% is between cell types *within* a
  library and therefore cannot be sequencing. Both the per-cell (unadjusted) and
  per-UMI (adjusted) estimands are reported; per-cell is primary because AD mast
  cells are transcriptionally smaller, which makes it the conservative choice.
* **Mast cells are re-derived, never trusted.** Deposited labels are checked
  against ≥2 of TPSAB1/TPSB2/CPA3; in cohorts without labels a magnitude
  threshold calibrated on the discovery cohort (96% precision) is added, because
  ambient tryptase otherwise tags hundreds of keratinocytes.
* **Negative results are reported in the same place as positive ones**, per
  CLAUDE.md §7 — including the failure of the discovery effect size to replicate
  and the failure of spatial co-localisation in lesional skin.
