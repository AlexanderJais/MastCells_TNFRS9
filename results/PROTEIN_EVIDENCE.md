# Does TNFRSF9 (4-1BB/CD137) protein data exist for skin mast cells?

**Short answer: no — not for cutaneous mast cells, in health or in atopic
dermatitis.** 4-1BB protein *has* been demonstrated on mast cells, and CD137
protein *has* been demonstrated on granulocytes in AD skin, but never on the
intersection of the two. Everything in `results/REPORT.md` is transcript-level,
and this gap is the single most consequential limitation of that report.

Searches run: PubMed (mast cell × 4-1BB/CD137/TNFRSF9; skin mast cell ×
co-stimulatory; 4-1BB × AD/skin), Human Protein Atlas API (ENSG00000049249),
GEO (CITE-seq skin, mass cytometry skin mast cell, AD surface proteomics).

---

## 1. What protein evidence exists

### 1.1 4-1BB protein on mast cells — established, but not human skin

**Nishimoto et al., *Blood* 2005;106(13):4241–8** — "Costimulation of mast cells
by 4-1BB, a member of the TNF receptor superfamily, with the high-affinity IgE
receptor."

* 4-1BB is induced **at mRNA *and protein* levels** on stimulation through FcεRI.
* Agonistic anti-4-1BB antibodies **enhance FcεRI-induced cytokine production**.
* **4-1BB-deficient mast cells show reduced degranulation and cytokine
  production**; 4-1BBL-deficient cells support the same conclusion.
* Mechanism: defective Ca²⁺ flux, reduced Lyn/Btk/PLCγ2 activity, and a
  constitutive 4-1BB–Lyn interaction.

This is the paper that matters most. It establishes protein expression *and* a
function, and it independently anticipates — twenty years earlier — the
FcεRI-driven induction recovered transcriptomically in §4 of the report.

**Caveats.** The loss-of-function work uses 4-1BB- and 4-1BBL-deficient mast
cells, which implies mouse bone-marrow-derived mast cells; the abstract does not
establish whether human mast cells were also assayed at protein level, and it is
not skin-derived mast cells. So this demonstrates *mast cells can express 4-1BB
protein on FcεRI stimulation*, not that human cutaneous mast cells do so in AD.

**Wensman et al., *Mol Immunol* 2012;50(4):210–9** — Tnfrsf9 was among the most
strongly upregulated genes in mast cells exposed to Lewis lung carcinoma
conditioned medium (array + qRT-PCR). Mouse, tumour context, transcript-level;
supporting but peripheral.

### 1.2 CD137 protein in AD skin — but on eosinophils, not mast cells

**Heinisch et al., *J Allergy Clin Immunol* 2001;108(1):21–8** — "Functional
CD137 receptors are expressed by eosinophils from patients with IgE-mediated
allergic responses but not by eosinophils from patients with non-IgE-mediated
eosinophilic disorders."

* CD137 protein measured by **flow cytometry and by immunohistochemistry in skin
  tissue**.
* Blood and tissue eosinophils from **atopic dermatitis** and extrinsic asthma
  express CD137; eosinophils from **healthy controls express neither mRNA nor
  protein**.
* Expression is induced in vitro by supernatants from activated T cells.

This is a striking parallel: in a different granulocyte, in the right disease and
the right tissue, CD137 shows exactly the pattern this project found for mast
cells — **absent at baseline, present in IgE-mediated disease, induced by
activation signals.** It also demonstrates that CD137 immunohistochemistry works
in human skin sections, which matters for the experiment proposed below.

### 1.3 Human Protein Atlas — no usable protein data, but notable RNA support

| Field | Value |
|---|---|
| Protein tissue specificity | **Not detected** |
| Protein tissue distribution | **Not detected** |
| Antibody | HPA071425 (single) |
| Reliability (IHC) | **None** — no validated immunohistochemistry |
| Reliability (IF) | Approved |
| Subcellular main location | Plasma membrane |

**HPA has no validated tissue protein data for TNFRSF9 anywhere in the body**, so
it cannot be used for skin. Any IHC attempt would need antibody validation first.

The same resource does, however, provide independent RNA support from its
pan-tissue single-cell atlas — where **mast cells are the highest-expressing
cell type for TNFRSF9**:

| Cell type | TNFRSF9 (nCPM) |
|---|---|
| **Mast cells** | **89.6** |
| T cells | 44.5 |
| NK cells | 26.3 |
| Innate lymphoid cells | 19.1 |
| cDC | 14.0 |
| Neutrophils | 12.1 |

This is orthogonal to every dataset analysed in the report and ranks mast cells
above T cells, which is worth stating plainly given how much of §7.3 is spent
excluding mast–T doublets. It is still RNA, and HPA aggregates across tissues and
disease states, so it speaks to capacity rather than to resting skin.

### 1.4 What does *not* exist

* **No CITE-seq or mass-cytometry dataset of human skin with a CD137 antibody.**
  GEO searches for skin CITE-seq, AD surface proteomics and skin mast-cell CyTOF
  returned nothing with CD137/4-1BB in the panel.
* **No flow-cytometry study of CD137 on human skin-derived mast cells**, in
  health or AD.
* **No 4-1BB immunohistochemistry co-stained with tryptase** in any skin
  condition.
* Systemic/soluble 4-1BB in AD serum is a separate quantity: measurable (TNFRSF9
  sits on the Olink Target 96 Inflammation panel, and AD systemic-biomarker
  proteomics exists), but soluble, systemic and not attributable to mast cells.
  The companion interstitial-fluid proteomics of GSE153760 is a plausible place
  to look for a skin-compartment soluble measurement.

---

## 2. What this does to the report's conclusions

**Strengthens.** The mechanism in §4 is not a transcriptomic curiosity: 4-1BB
protein is induced on mast cells by FcεRI cross-linking and is functionally
costimulatory (Nishimoto 2005). The AD-restricted, activation-induced pattern is
independently documented for CD137 protein in AD skin in a neighbouring
granulocyte (Heinisch 2001). And HPA's single-cell atlas independently ranks mast
cells as the top TNFRSF9-expressing cell type.

**Does not change.** The in-vivo effect size remains transcript-level, modest and
heterogeneous (log₂FC +1.32, 95% CI −0.07 to +2.71, I² = 42%, 843 vs 4,315 mast
cells). No protein measurement exists to corroborate or contradict it in skin.

**Remains open.** Whether 4-1BB reaches the surface of human cutaneous mast cells
in vivo; how many mast cells are positive in AD versus healthy skin; and whether
surface density is sufficient for the costimulatory function Nishimoto describes.

---

## 3. The experiment that would close the gap

Ordered by cost, with the specific reagents:

1. **Multiplex IF/IHC on archival FFPE skin** — tryptase (or CPA3) + CD137
   co-stain on AD lesional, AD non-lesional and healthy sections. Cheapest, uses
   existing blocks, and Heinisch 2001 establishes that CD137 IHC works in skin.
   Requires antibody validation first: the sole HPA antibody (HPA071425) has no
   IHC reliability score, so a clinically used clone (e.g. 4B4-1) with an
   isotype control and a CD137⁺ positive-control tissue is the safer route.
2. **Flow cytometry of enzymatically dispersed skin** — gate
   CD45⁺CD117⁺FcεRIα⁺ mast cells, stain CD137, compare healthy versus AD, with
   and without ex-vivo IL-33 or anti-IgE. Directly tests the report's central
   prediction: near-zero at rest, induced on stimulation. Note the report's own
   finding that mast-cell recovery from skin ranges 0.11–4.4% depending on
   protocol — the dissociation method will dominate the yield.
3. **CITE-seq with TotalSeq anti-CD137 on AD skin biopsies** — gives protein and
   transcript in the same cell and would settle the doublet question definitively,
   but is the most expensive option and inherits the same dissociation problem.

The prediction to test, stated so it can fail: **4-1BB protein should be
undetectable or near-undetectable on mast cells in healthy skin, detectable on a
minority of mast cells in AD lesional skin, and strongly inducible on
skin-derived mast cells within hours of IL-33 or anti-IgE ex vivo.**
