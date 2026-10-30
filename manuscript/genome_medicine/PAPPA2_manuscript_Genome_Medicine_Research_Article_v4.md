# Unit-aware analysis reveals opposing donor and within-donor disease signals in human kidney single-cell data

**Article type:** Research Article

**Authors:** Peng Zhi, Xiaohe Yan, Chaofan Zhang

**Affiliations:** Henan Provincial Chest Hospital, Zhengzhou, China

**Corresponding author:** Peng Zhi, Department of General Practice, Henan
Provincial Chest Hospital, Zhengzhou, Henan, China. Email:
zhiyoyo0409@163.com. ORCID: 0009-0003-2201-5131

**Co-author contacts:** Xiaohe Yan. Email: yanxiaohe123456@163.com. ORCID:
0009-0006-2259-9271; Chaofan Zhang. Email: zcf1765081785@163.com. ORCID:
0009-0009-9088-0016

## Abstract

Single-cell disease studies frequently describe a cell state as
disease-associated without specifying whether the signal reflects donor
composition, local neighborhood structure, or coupling within a donor. We
developed a unit-aware framework that reports these quantities separately and
applied it to human kidney thick ascending limb (TAL) states. We analyzed
106,851 KPMP TAL cells from 223 donors using donor pseudobulk, Milo
neighborhood abundance, neighborhood-expression association, and within-donor
cell-rank association. In aTAL2, donor pseudobulk and neighborhood-expression
direction were positive, whereas within-donor association was negative in
combined disease and AKI (FDR = 0.0086 for each). In a prespecified 400-gene
panel, donor-within reversal was 19.0%, compared with 12.7% in
expression-matched random sets. A conditional effect-pairing null gave 20.8%,
showing that the initial excess reflected the selected effect-size
distribution rather than negative coupling between units. Simulations
distinguished opposing from concordant relationships and quantified a
screening-threshold null. Two non-KPMP cohorts supported
donor-level PAPPA2 disease elevation but not the within-donor effect.
Donor-stratified signature transfer achieved 93.9% held-out accuracy. PAPPA2
was a strong kidney tubulointerstitial eGene (FDR = 1.64e-14) without
glomerular eQTL support. Decomposition attributed the overall PAPPA2
difference to subtype-expression rather than subtype-composition effects. The
analysis unit is therefore an effect modifier in single-cell disease studies.
Donor-aware and within-donor findings should be reported separately, and
matched-set comparisons should be accompanied by an effect-aware null.

**Keywords:** single-cell genomics, analysis unit, kidney disease, thick
ascending limb, donor-aware inference, software framework

## Introduction

Chronic kidney disease and acute kidney injury are major causes of morbidity
and mortality worldwide [1,2]. Kidney injury is not distributed uniformly
across the nephron. Proximal tubules, thick ascending limb (TAL) cells, distal
tubules, endothelial cells, and immune populations can adopt distinct injury,
adaptation, and repair states [3,4]. The biological and clinical meaning of
these states depends on both tissue context and the way in which the disease
signal is measured. The central question is not only whether a gene is
associated with disease, but whether that association survives a change in
the unit of inference.

Single-cell and spatial technologies now make it possible to map kidney cell
states at high resolution [5]. However, the conclusions of such studies depend
on the unit used for inference. Donor pseudobulk asks whether mean expression
differs between disease and Reference donors; Milo neighborhood abundance asks
whether neighborhood abundance differs between disease and Reference;
neighborhood-expression association asks whether neighborhoods with higher
expression also have higher disease logFC; and within-donor analysis asks
whether cells with greater neighborhood exposure have higher expression inside
each donor. Donor-level robustness is used as a stability check rather than as
a separate unit of inference. These approaches answer different questions and
may not give the same answer. This ambiguity is important because a disease
signal attributed to one cell state could instead reflect donor composition,
local neighborhood structure, or a difference between cells within the same
donor. Direct comparison of established analysis units within the same disease
atlas is therefore needed.

TAL cells perform active sodium, potassium, and chloride transport, are
metabolically demanding, and are sensitive to ischemic and inflammatory injury
[6,7]. Several transcriptionally distinct TAL subtypes have been reported in
human kidney atlases, but whether these subtypes share the same disease-related
neighborhood states remains unclear. This question is clinically relevant
because TAL injury and adaptation may affect tubular workload, inflammatory
signaling, and repair programs. PAPPA2 is not the central hypothesis; it is a
case gene with a less well-defined renal role, used to test whether a
gene-state association remains consistent when the analysis unit changes
[8,9,10].

In this study, we analyzed the complete KPMP single-nucleus resource, including
106,851 TAL cells, and compared four established units of inference. The study
had three prespecified aims: first, to determine whether donor, neighborhood,
neighborhood-expression, and within-donor analyses produce concordant
subtype-level disease conclusions; second, to test whether the principal
directions transfer to non-KPMP kidney cohorts; and third, to determine
whether the PAPPA2 case signal has a human kidney regulatory context. The study
was not designed as a new algorithm, causal test, or mechanistic experiment.
It was designed to determine which conclusions survive a change in the unit of
inference.

The central finding is that changing the analysis unit changes subtype-level
conclusions. aTAL2 provided the strongest example of cross-unit discordance,
frTAL showed a weaker related pattern, and PAPPA2 marked a disease-associated
TAL state without establishing causality. We then tested the key directions in
two non-KPMP cohorts and evaluated kidney tubulointerstitial eQTL and snATAC
context. These analyses separate a stable disease-level observation from claims
that require exact subtype replication, chromatin validation, or causal
inference.

## Methods

### Study design and data sources

We designed a framework-development and empirical validation study using
public or controlled-access human kidney data. No new human samples were
collected, and no new ethical approval was required for secondary analysis of
de-identified data. The study implemented a reusable unit-aware workflow and
tested whether disease conclusions remained stable when the unit of inference
changed, rather than proposing a new molecular mechanism (Fig. 1).

[[FIGURE:Fig1]]

We selected data sources to provide one discovery resource and separate
directional support. The discovery resource was the KPMP human kidney
single-nucleus RNA sequencing dataset [5], which contained 1,388,643 cells and
35,455 features. We extracted all cells annotated as C-TAL-A, C-TAL-B, aTAL2,
or frTAL, yielding 106,851 TAL cells. A previously generated 4,446-cell
derivative matrix was retained only as a coverage and directional comparator.
Spatial cross-modality analysis used the KPMP GSE183456 Visium resource, part
of the KPMP GEO superseries GSE183279, which contained 23 sections from 22
donors [5,11]. Because GSE183456 shares KPMP consortium provenance with the
discovery data, it was treated as same-consortium support rather than an
external cohort. Separate directional eGFR testing used GSE137570 [12], and
multi-dataset GFR association testing used Nephroseq v5 [13]. Independent
validation used the non-KPMP GSE211785 multi-omic kidney atlas and the
independent GSE254185 kidney multiome resource.

### TAL subtypes and the four analysis units

We restricted the analysis to four TAL subtypes and defined the analysis units
before comparing their results. The subtypes were C-TAL-A, C-TAL-B, aTAL2, and
frTAL. All primary analyses used the complete 106,851-cell TAL set and the
atlas-provided UMAP coordinates.

We compared four primary analysis units because each answers a different
question about disease-related TAL state. Donor pseudobulk asks whether mean
expression differs between disease and Reference donors. Milo neighborhood
abundance asks whether neighborhood abundance differs between disease and
Reference. Neighborhood-expression association asks whether neighborhoods
with higher expression also have higher disease logFC. Within-donor
cell-rank association asks whether cells with greater neighborhood exposure
have higher expression inside each donor.

We treated donor-level robustness as a stability check rather than as a fifth
analysis unit. This distinction avoids treating resampling stability as an
independent biological measurement. The primary Milo threshold was SpatialFDR
< 0.05. Milo used logMS normalization and a fixed random seed of 20260927. One
AKI donor with zero neighborhood counts was excluded from that comparison.
Augur was retained only as a secondary subset-level analysis of disease versus
Reference separability [14], using 20 subsamples, 30 cells per subsample, and
3-fold cross-validation.

### Donor-level pseudobulk expression

Donor pseudobulk provided the reference analysis unit, because donor-level mean
expression is the most common way to relate cell-state genes to disease.
Raw counts for PAPPA2, CXCL12, NOTCH3, HES1, and EPHB2 were extracted directly
from the sparse source matrix. Counts were normalized using the full per-cell
RNA count and transformed with log1p. Expression was averaged within each
donor and TAL subtype. Disease groups were compared with Reference donors using
two-sided Mann-Whitney tests, and donor-level linear models were adjusted for
age midpoint and sex. Benjamini-Hochberg correction was applied across the 60
gene, subtype, and comparison tests. The 4,446-cell subset result was retained
only as a comparator, and its original 120-test correction was not mixed with
the complete-TAL analysis.

### Composition and expression decomposition

We decomposed the overall four-subtype TAL PAPPA2 difference into a subtype
composition component and a subtype-expression component. For each donor, the
overall TAL mean was the cell-weighted sum of subtype proportions and
subtype-specific mean expression. Group differences were separated into
composition at reference means, expression at reference proportions, and an
exact remainder containing interaction and within-group covariance terms.
Missing subtype means were set to zero contribution when the corresponding
proportion was zero. Confidence intervals were obtained from 2,000 donor-level
bootstrap resamples within disease and Reference groups.

### Milo neighborhood abundance

Milo neighborhood abundance tested whether disease was associated with shifts
in local TAL neighborhoods rather than only with donor-level average
expression. Neighborhoods were rebuilt from the complete TAL set using the
atlas UMAP coordinates, k = 30, a neighbor ratio of 0.1, and refined sampling.
We tested disease combined, AKI, and CKD against Reference [15]. Results were
interpreted as differences in TAL-associated neighborhood state rather than as
evidence of direct cell-cell contact.

### Neighborhood-expression association

Neighborhood-expression association tested whether neighborhoods with higher
gene expression were also those with higher disease logFC. Each Milo
neighborhood was represented by its mean PAPPA2 expression and its Milo logFC,
and Spearman correlation was calculated across neighborhoods for each
comparison and subtype. Stability was assessed with 500 donor-stratified
resamples, leave-one-donor-out analysis, and a partial Spearman correlation controlling for neighborhood purity (Additional file 1). These resampling intervals describe
donor-level stability and are not independent cell-level p values.

### Within-donor cell-rank association

Within-donor cell-rank association tested whether cells with greater
neighborhood exposure had higher expression inside the same donor. For each
cell, neighborhood exposure was defined as the mean Milo logFC of the
neighborhoods assigned to that cell's TAL subtype. Within each donor, Spearman
correlation was calculated between PAPPA2 expression and neighborhood exposure,
so each donor contributed one estimate. Statistical support was assessed with
1,000 within-donor expression permutations. The primary analysis required at
least 5 cells per donor, and sensitivity analyses required at least 3 and 10
cells. BH correction was applied across 60 gene, subtype, and comparison tests.
Because neighborhoods overlap, conventional correlation p values were treated
as descriptive rather than as independent-sample inference.

### Analysis-unit comparison and internal robustness

Analysis-unit comparison tested whether the same gene, subtype, and disease
comparison led to the same direction across the four units. For each gene,
subtype, and comparison, we recorded the direction in donor pseudobulk, Milo
neighborhood abundance, neighborhood-expression association, and within-donor
cell-rank association and classified the result as concordant or discordant. A
1,000-replicate donor bootstrap was used for complete-TAL pseudobulk effects.
The complete-TAL pseudobulk analysis, the previous subset pseudobulk analysis,
the neighborhood-level resampling result, and the within-donor permutation
result were joined into a single analysis-unit table. This table was used to
evaluate stability across units rather than to treat the most significant unit
as definitive.

Internal robustness tested whether the donor-level directions were dependent
on individual donors, age or sex, or tissue region. We repeated the
pseudobulk analysis after removing one donor at a time, resampling donors
within condition strata, adjusting for age and sex, and stratifying by cortex
or medulla. These analyses were treated as stability checks and not as
additional discovery units.

### Prespecified 400-gene analysis-unit panel

A 400-gene panel was frozen before the panel outcome analysis. The panel
contained 100 TAL subtype markers, 100 MSigDB Hallmark injury, inflammatory,
hypoxic, or metabolic genes, 100 genes coexpressed with PAPPA2 in the 156
training donors, and 100 internal expression-matched negative controls. A
separate 2,000-gene pool was matched to the target expression and detection
distribution and used to generate 1,000 random 400-gene sets.

Only the three gene-specific units were compared across this panel: donor
pseudobulk, neighborhood-expression association, and within-donor cell-rank
association. A reversal required opposite signs and absolute standardized
effects of at least 0.5 in both units. The primary endpoint was the
donor-versus-within reversal fraction. One-sided empirical P values compared
the target panel with the random sets. An expression-stratified random-set
analysis was performed as a sensitivity check. Because expression matching
does not match effect-size distributions, a post hoc conditional pairing
permutation was then performed with 10,000 permutations. Within each subtype
and comparison, right-unit effects were reassigned across genes without
replacement while the left-unit effect and both effect-size marginals were
preserved. Valid strata were averaged within each gene before genes were
averaged, matching the original panel endpoint. One-sided empirical P values
used the 10,000 conditional null values. Threshold sensitivity was evaluated
from 0.5 to 2.0. The target-panel hash and frozen definitions are provided in
Additional file 21.

### Unit-aware software framework

We implemented the analytical workflow as `unitaware` v0.2.0, a
platform-independent Python package. The package provides functions for donor pseudobulk,
neighborhood-expression association, within-donor cell-rank association,
standardized reversal scoring, composition-expression decomposition, and
simulation. It also provides a conditional effect-pairing permutation for
calibrating reversal estimates against the observed effect-size structure.
Analysis-unit definitions and thresholds are explicit inputs rather than
hidden defaults. The package includes a command-line interface, six unit tests,
API documentation, and an MIT software license, and is distributed with the
archived analysis code.

### Simulation benchmark

We simulated six prespecified scenarios: null, donor-positive only,
within-positive only, concordant positive, donor positive with within
negative, and donor negative with within positive. Each scenario used 25
Reference and 25 disease donors, 30-80 cells per donor, gene-specific effects,
donor random intercepts, 50 exposure-defined neighborhoods, and 100 replicates.
The same reversal threshold used in the empirical panel was applied. The
simulation tested whether the framework distinguished genuine opposing effects
from null and concordant settings.

### KPMP eGFR analysis

The KPMP eGFR analysis tested whether the PAPPA2-associated TAL state was
linked to kidney function at the donor level. Donor-level PAPPA2-positive
fraction was calculated within TAL cells. Baseline eGFR was treated as an
ordinal outcome with four categories: 0-29, 30-59, 60-89, and 90 or more
mL/min/1.73 m². The primary ordinal logistic model included
PAPPA2-positive fraction per 10 percentage points, age, sex, and disease
category. Additional models added C-TAL-A, C-TAL-B, aTAL2, and frTAL fractions,
followed by mean detected genes and assay platform. Odds ratios below 1
indicate lower eGFR with higher PAPPA2-positive fraction.

### Candidate-gene specificity

Candidate-gene specificity tested whether the PAPPA2 direction was distinctive
relative to pre-specified control genes. CXCL12, NOTCH3, HES1, and EPHB2 were
analyzed with the same donor-aware and eGFR models used for PAPPA2. The
comparison assessed directional consistency, BH-corrected support, and tissue
specificity. Because multiple genes can correlate with GFR, the aim was not to
identify a unique biomarker but to determine whether PAPPA2 showed a distinct
negative direction among the tested candidates.

### KPMP GSE183456 spatial cross-modality analysis

Spatial cross-modality analysis tested whether the PAPPA2-associated TAL state
co-occurred with a repair-failure spatial program within the same KPMP
consortium resource. The frozen spatial analysis used nine main-branch genes:
THSD4, TEX41, EGF, CA10, CACNB2, CABP1, CXCL12, PRKG1, and NOS1. The
repair-failure signature used TMEM207, KLRD1, RNF144B, JAG1, TMEM52B, DACH1,
OSBPL10, LRMDA, and KCNJ10. For each section, genes were normalized by counts
per 10,000, transformed with log1p, and converted to within-section z scores.
Branch scores were the mean z scores of the signature genes, and high-state
thresholds were fixed at the 75th percentile. The primary metric was the odds
ratio for co-occurrence of high main-branch and high repair-failure states.
The same parameters were applied to all 23 sections without re-tuning.

### GSE137570 public eGFR testing

GSE137570 testing provided independent but underpowered directional evidence
for the PAPPA2-eGFR relationship. GSE137570 Cohort 1 contained 24 donors with
public eGFR, age, sex, and fibrosis. PAPPA2 expression was extracted from the
normalized read-count matrix. Spearman correlation was calculated between
PAPPA2 and eGFR, and a linear model adjusted for age, sex, and fibrosis. The
same model was applied to the four control genes. This analysis was regarded
as directional support rather than definitive external validation.

### Nephroseq GFR analysis

Nephroseq analysis tested whether the PAPPA2-GFR direction was observed across
multiple precomputed kidney datasets. Nephroseq v5 was queried through its
public API using the GFR filter and thresholds pValue < 0.05,
abs(rValue) ≥ 0.5, and abs(foldChange) ≥ 1.5. PAPPA2 and the four control
genes were queried with the same thresholds. Each returned analysis was
summarized by dataset, tissue type, GFR method, sample count, correlation
coefficient, and p value. Analyses were not treated as independent cohorts
because several analyses can originate from the same dataset.

### Independent non-KPMP cohort analyses

GSE211785 was generated by the Susztak laboratory and contains control kidney,
diabetic kidney disease (DKD), and hypertensive kidney disease samples.
We used the non-KPMP merged object and restricted analyses to cells annotated
as C_TAL or M_TAL in the SC_RNA and SN_RNA libraries. After requiring at least
five TAL cells per donor, 29,847 TAL cells from 35 donors with TAL coverage
were available (16 control, 7 DKD, and 12 HKD donors).

GSE254185 contains 46,957 nuclei from 7 normal kidneys and 5 kidneys with
obstructive nephropathy. Its public TAL annotations included cTAL1,
cTAL2, mTAL, TAL Injured, and TAL Inflammatory states, providing 9,912 annotated
TAL cells. For both cohorts, counts were normalized to counts per 10,000 and
log1p transformed. Donor pseudobulk values were calculated within TAL state.
Donor-label permutations used 20,000 random permutations unless exact
enumeration was feasible. Neighborhood features were calculated from the
public UMAP coordinates using k = 30 nearest neighbors. These analyses
provide independent-cohort validation, not experimental or causal validation.

### Donor-stratified TAL signature transfer

To test whether the four KPMP TAL signatures could identify an aTAL2-like
state outside the discovery annotations, donors were divided into 156 training
and 67 held-out donors. The top 35 markers per subtype were selected only in
the training donors. Each cell received the mean within-subtype marker z
score, and the difference between the top and second-ranked subtype score was
defined as the marker margin. The prespecified high-confidence threshold was
a marker margin of at least 0.2.

Held-out KPMP cells were used to estimate coverage, accuracy, macro F1, and
aTAL2 precision and recall. The same frozen marker panel was then projected to
GSE211785 and GSE254185. In each external cohort, donor pseudobulk and
within-donor analyses were repeated in the predicted high-confidence
aTAL2-like cells. This analysis tests transfer of a transcriptional signature.
It does not reproduce the original KPMP subtype labels and is not described as
exact aTAL2 validation.

### Kidney eQTL and snATAC context

We used NephQTL2 gene-level eQTL results from 240 glomerular and 311
tubulointerstitial human kidney samples [16]. PAPPA2 and the fixed
top-35 aTAL2 markers were evaluated without re-selecting genes. We reported
tubulointerstitial and glomerular gene-level FDR values and checked whether
PAPPA2 appeared among FDR-significant PTWAS associations for eGFR or urinary
albumin-to-creatinine ratio. Absence from the significant PTWAS list was
recorded as non-significance, not as a negative association.

GSE211785 snATAC nuclei annotated as C_TAL or M_TAL were analyzed using the
gene-activity representation as a chromatin-accessibility proxy. Donor means
were calculated after requiring at least five snATAC nuclei per donor and
state. Disease and control donors were compared with exact donor-label
permutation when feasible. Seven control and four disease donors contributed
to both states. Benjamini-Hochberg correction was applied across the fixed gene
and state tests. GSE211785 RNA and snATAC are from the same cohort, so the
snATAC result is an orthogonal-modality analysis rather than an independent
cohort.

### Reproducibility and software

Reproducibility was addressed by fixing analysis parameters and linking each
reported result to versioned code and summary tables. Analysis scripts, frozen
analysis protocols, validation tables, and result hashes were prepared for
release in the public repository. The fixed-seed Milo analysis, GSE183456
spatial analysis, and Nephroseq queries were rerun with fixed parameters.
Random seeds and versioned outputs are listed in the Availability of data and
materials section.

The authors used OpenAI Codex CLI to assist with organization of
methodological text, reference checking, and English-language editing. The
authors reviewed and edited all content and retain responsibility for the
final manuscript.

## Results

### Complete TAL coverage provided a common substrate for analysis-unit comparison

Complete TAL coverage allowed all four analysis units to be evaluated in the
same cells. The complete set contained 106,851 cells from the four target TAL
subtypes, whereas the derivative comparator contained 4,446 cells. The
coverage difference was largest for subtypes that were sparse in the
derivative matrix: aTAL2 contained 6,624 Reference cells from 83 donors in the
complete data but only 102 Reference cells from 38 donors in the subset; frTAL
contained 6,656 Reference cells from 86 donors in the complete data but only
142 cells from 48 donors in the subset. The complete set therefore provided a
different donor and state composition from the derivative matrix and served as
the primary substrate for all subsequent comparisons.

### TAL subtypes did not share a single disease-related neighborhood state

Complete-data Milo showed that disease-related neighborhood changes were
strongly subtype dependent. In the combined disease comparison, C-TAL-A had
1,797 significant neighborhoods among 4,012 annotated neighborhoods, C-TAL-B
had 715 among 2,805, aTAL2 had 632 among 1,465, and frTAL had 1,023 among
1,332. The direction also differed: frTAL had 996 positive and 27 negative
significant neighborhoods, whereas C-TAL-A had 333 positive and 1,464 negative
neighborhoods.

Disease-specific analyses preserved this subtype dependence. In AKI, aTAL2 had
842 significant neighborhoods, 822 positive and 20 negative, with a median
logFC of 1.364. In CKD, frTAL had 982 significant neighborhoods, 956 positive
and 26 negative, with a median logFC of 1.138. By contrast, CKD aTAL2 had a
median logFC of -0.406. These results show that TAL neighborhood remodeling
cannot be summarized as one shared disease response (Fig. 2;
Additional files 2 and 3).

[[FIGURE:Fig2]]

The secondary subset-level Augur analysis provided complementary evidence that
subtype states were not equally separable. AKI C-TAL-A had the highest mean
AUC at 0.753, whereas CKD aTAL2 had a mean AUC below 0.5 at 0.479. The moderate
AUC values and broad subsample intervals indicate that disease and Reference
states were only partially separable. Together with Milo, these results
establish subtype-specific neighborhood structure before the analysis units are compared directly (Fig. 3).

[[FIGURE:Fig3]]

### Subset and complete pseudobulk results changed direction for selected genes

The derivative subset preserved broad pseudobulk trends but changed selected
gene and subtype conclusions. In the complete analysis, 20 of 60 gene, subtype,
and comparison tests had FDR < 0.05 using unadjusted Mann-Whitney tests, and 25
of 60 remained significant after adjustment for age and sex. Directions agreed
between the subset and complete pseudobulk analyses in 48 of 60 comparisons,
whereas 12 comparisons changed direction. Those discordant tests were
concentrated in CXCL12, HES1, NOTCH3, and EPHB2, and none involved PAPPA2.
Thus, the subset result was useful as a comparator but was not sufficient to represent the complete TAL population (Fig. 4; Additional file 4).

[[FIGURE:Fig4]]

### The four analysis units were not interchangeable

The four analysis units produced related but non-equivalent views of TAL state.
Across 12 comparison and subtype groups, neighborhood PAPPA2 expression was
positively correlated with Milo logFC in all 12 groups, with Spearman
correlations ranging from 0.029 in AKI aTAL2 to 0.729 in combined-disease
C-TAL-A. Donor-stratified bootstrap intervals excluded zero in 11 of 12 groups;
the exception was AKI aTAL2 (95% interval -0.086 to 0.255). Partial
correlations controlling for neighborhood purity were similar to the primary
correlations (Fig. 5; Additional files 5 and 6).

[[FIGURE:Fig5]]

The units differed in what they measured and did not always support one shared direction (Additional file 7). Donor pseudobulk, neighborhood-expression association, and
within-donor cell-rank association were therefore analyzed as distinct units
rather than as interchangeable estimates of the same quantity. This difference
was clearest in aTAL2 and was also present in a weaker form in frTAL.

### Analysis-unit reversal was frequent but depended on effect structure

The donor-versus-within reversal fraction was 19.0% in the 400-gene target
panel, compared with 12.7% in expression-matched random sets (empirical
P = 0.000999, the minimum attainable value with 1,000 sets). Restricting to the
300 biological genes increased the fraction to 20.8% versus 12.7%. An
expression-stratified sensitivity analysis gave the same result (19.0% versus
12.8%; empirical P = 0.000999).

The panel result was not driven by the internal negative controls. Their
donor-versus-within reversal fraction was 13.6%, statistically
indistinguishable from random 100-gene sets from the external control pool
(12.6%; two-sided empirical P = 0.485). Reversal fractions differed across
biological groups: 24.0% for PAPPA2 coexpression genes, 21.6% for TAL subtype
markers, and 16.8% for injury, inflammatory, and metabolic genes.

Donor-versus-neighborhood reversal was also higher (21.9% versus 14.5%;
empirical P = 0.000999). Neighborhood-versus-within reversal showed a smaller
full-panel excess (21.0% versus 19.6%; empirical P = 0.005), but the biological
300-gene subset did not (P = 0.189).

The conditional pairing null changed the interpretation. For
donor-versus-within reversal, the observed fraction was 19.0% and the
effect-aware null was 20.8% (P = 0.999). For
donor-versus-neighborhood reversal, the observed fraction was 21.9% and the
null was 21.2% (P = 0.111). Thus, the apparent excess over expression-matched
genes was explained by selection of gene-state pairs with larger effects in
both units. The selected panel demonstrated frequent sign disagreement, but
did not establish broad negative coupling beyond the effect-size distribution.

Neighborhood-versus-within reversal retained a small conditional excess
(21.0% versus 19.4%; P = 0.0017). This result was threshold sensitive and was
not significant above the 0.5 threshold. Individual gene metrics, random-set
results, conditional null distributions, and threshold sensitivity are
provided in Additional files 22-24 and 26-29. The analysis therefore supports
unit dependence as a measurement property, not a claim that selected TAL genes
are uniquely prone to reversal.

### Simulation benchmark distinguished opposing from concordant unit effects

In the simulation benchmark, donor-versus-within reversal was 95.5% when the
true donor effect was positive and the within-donor effect negative, and 97.3%
in the reciprocal scenario. Concordant positive effects produced a 2.3%
reversal fraction, donor-only and within-only scenarios produced intermediate
fractions, and the null scenario produced 20.1%. The null fraction reflects the
prespecified effect-size threshold and finite sampling, not a true opposing
effect. The value is therefore a screening-threshold calibration rather than
a biological false-positive rate. These results show that the framework can
recover strong opposing relationships while retaining explicit null calibration
(Additional files 24, 25, and 27).

### aTAL2 showed the strongest cross-unit discordance

aTAL2 was the clearest example because donor-level and within-donor analyses
pointed in opposite directions. Complete donor pseudobulk showed higher PAPPA2
expression in aTAL2 for all three comparisons, with positive mean differences
and FDR < 0.05 in each case. The same positive direction was retained in the
age- and sex-adjusted models.

Neighborhood-expression association was positive but less stable in aTAL2.
The full-data Spearman correlation was 0.145 for combined disease and 0.029
for AKI, and the AKI donor-stratified bootstrap interval included zero. This
result indicates that neighborhood expression and donor pseudobulk were
aligned in direction but not equally strong.

Within-donor cell-rank association moved in the opposite direction. The
disease-minus-Reference rank difference was -0.2857 in AKI and -0.2063 in
combined disease, with FDR = 0.0086 for each comparison. In CKD, the difference
was also negative (-0.0886) but did not pass FDR (0.224). Therefore the aTAL2
conflict was strongest in aTAL2 and significant in two of three disease
comparisons, rather than uniform across all comparisons (Fig. 6; Additional file 8).

[[FIGURE:Fig6]]

### frTAL showed a weaker related cross-unit pattern

frTAL showed the same qualitative direction conflict as aTAL2, but with weaker
statistical support. Donor pseudobulk and neighborhood-expression association
were positive for frTAL. In contrast, the within-donor disease-minus-Reference
rank difference was negative in all three comparisons: -0.0897 in combined
disease (FDR = 0.018), -0.1175 in AKI (FDR = 0.146), and -0.0557 in CKD
(FDR = 0.224). Thus, frTAL supports the general analysis-unit sensitivity
finding but does not reach the same strength as aTAL2.

### PAPPA2 was a directionally consistent case gene across TAL subtypes

Complete donor-level pseudobulk analysis showed higher PAPPA2 expression in all
four TAL subtypes for AKI, CKD, and combined disease comparisons. All 12
PAPPA2 gene, subtype, and comparison combinations had positive effect
direction, and 11 of 12 passed FDR < 0.05; AKI frTAL was directionally positive
but did not pass FDR (0.080). For example, the mean increase in the combined
disease comparison was 0.894 for C-TAL-A, 0.477 for C-TAL-B, 0.411 for aTAL2,
and 0.370 for frTAL. The same positive pattern was observed in the age- and sex-adjusted models (Additional file 9).

Leave-one-donor analysis preserved the PAPPA2 direction in all 12 groups (Fig. 7; Additional file 10).
Bootstrap intervals were used to describe donor-level stability, not
independent-cell significance. Cortex and medulla strata agreed in direction
in 11 of 12 groups; the only exception was AKI frTAL in medulla, a group that
was not significant in the complete pseudobulk analysis. These results
identify PAPPA2 as a consistent disease-associated TAL case gene while
preserving the distinction between average disease-group elevation and
within-donor state association.

[[FIGURE:Fig7]]

### PAPPA2 enrichment reflected subtype expression rather than composition

The overall four-subtype TAL PAPPA2 difference was not driven by subtype
proportion changes. In the combined disease comparison, the observed difference
was +0.563 (95% bootstrap interval +0.483 to +0.639), the subtype-expression
component was +0.641 (+0.560 to +0.721), the composition component was +0.010
(-0.003 to +0.023), and the interaction/covariance remainder was -0.088
(-0.133 to -0.044). The expression component accounted for 86.7% of the
absolute component magnitude, compared with 1.3% for composition.

The same pattern was present in disease-specific comparisons. In AKI, the
observed difference was +0.379, expression was +0.522, and composition was
-0.009. In CKD, the observed difference was +0.650, expression was +0.698, and
composition was +0.019 (Additional file 20). Thus, higher PAPPA2 in disease
donors was not explained by a shift toward TAL subtypes with higher baseline
PAPPA2. The between-donor expression difference and the negative within-donor
neighborhood association therefore describe separable levels of variation.
This decomposition does not establish cell-intrinsic induction; donor-level
exposures or unmeasured tissue context may still contribute.

### PAPPA2-positive fraction was associated with lower eGFR in KPMP

In KPMP donors, higher TAL PAPPA2-positive fraction was associated with lower
eGFR in ordinal models. The primary model included 202 donors. For every 10
percentage-point increase in PAPPA2 fraction, the odds ratio was 0.819
(95% CI 0.682-0.985, p = 0.0335) after adjustment for age, sex, and disease
category. The direction remained negative after adding TAL subtype fractions
(OR 0.696, 95% CI 0.530-0.915) and technical covariates
(OR 0.661, 95% CI 0.500-0.875).

The signal was concentrated mainly in CKD. The primary model in disease donors
included 130 donors and produced an OR of 0.789 (95% CI 0.648-0.960). CKD
included 94 donors and produced an OR of 0.733 (95% CI 0.582-0.924). AKI
included 36 donors and did not show a clear association. These results support
an internal adjusted association between PAPPA2 and kidney function, but they
do not constitute external clinical validation (Additional file 11).

### The PAPPA2 signal was directionally distinct among tested genes

In the within-donor analysis, PAPPA2 had seven FDR-supported comparisons at the
primary threshold, compared with two for CXCL12 and none for NOTCH3, HES1, and
EPHB2. The advantage of PAPPA2 was primarily directional consistency rather
than the largest effect size.

The same pattern was observed in Nephroseq. PAPPA2 had seven significant
tissue-specific negative GFR analyses, compared with six for HES1, three for
EPHB2, two for NOTCH3, and one for CXCL12. HES1 and EPHB2 had multiple
significant results, but their directions were mixed or less restricted to
tubulointerstitium. CXCL12 did not reproduce the negative PAPPA2 pattern.

### PAPPA2-associated spatial signatures showed directional consistency in the KPMP Visium resource

GSE183456 is the KPMP Visium spatial resource and belongs to the same
consortium as the discovery snRNA-seq data. We therefore used it for spatial
cross-modality analysis, not as an independent external cohort. Across 23
sections, 20 had a co-occurrence odds ratio above 1 for the two
PAPPA2-associated branch signatures. All 11 DKD sections had an odds ratio
above 1, and all 10 DKD donors remained above 1 after donor-level
deduplication. The result was directionally consistent with the discovery
analysis.

The spatial result was not disease-specific. The DKD versus Reference
comparison had p = 0.078, and 5 of 6 Reference sections also had a
co-occurrence odds ratio above 1. The appropriate interpretation is
cross-modality directional consistency of a tissue-state relationship within
the KPMP resource, not a disease-specific or external spatial validation.

### Public external datasets provided directional PAPPA2-GFR support

GSE137570 Cohort 1 provided public eGFR for 24 donors. PAPPA2 was negatively
associated with eGFR (Spearman rho = -0.365, p = 0.079), but the adjusted model
was not significant (p = 0.351). EPHB2 had a stronger unadjusted correlation
(rho = -0.646, p = 0.0006), indicating that this cohort cannot support a claim of PAPPA2-specific significance (Additional file 12).

Nephroseq provided stronger external association support. In 50 GFR analyses,
PAPPA2 had 10 significant results, of which 9 were negative. Seven negative
analyses came from tubulointerstitial datasets. The strongest results included
Sampson Nephrotic Syndrome TubInt (r = -0.558, n = 49), Ju CKD TubInt 2
(r = -0.530, n = 42), ERCB Nephrotic Syndrome TubInt (r = -0.712, n = 17; and
r = -0.686, n = 16), and Mariani Nephrotic Syndrome TubInt 2 (r = -0.837,
n = 10). One positive analysis came from a cadaveric transplant kidney dataset
(r = 0.589, n = 12).

Together, GSE137570 and Nephroseq support the direction of the PAPPA2-GFR
relationship. The Nephroseq analyses are precomputed and not independently
adjusted in a common model, so they should be reported as multi-dataset
external association support rather than as independent cohort counts (Additional file 13).

### Non-KPMP cohorts independently supported PAPPA2 disease elevation

In GSE211785, donor-level PAPPA2 expression was higher in disease than control
donors in both TAL compartments. For C_TAL, the difference was +0.283
(permutation p = 0.0315; FDR = 0.0630 across six TAL-state and comparison
tests), with a bootstrap 95% CI of +0.053 to +0.514. For M_TAL, the difference
was +0.200 (p = 0.105; FDR = 0.133), with a bootstrap 95% CI of +0.026 to
+0.400. DKD-specific contrasts were larger (+0.443 in C_TAL, p = 0.0174;
+0.414 in M_TAL, p = 0.0227), whereas HKD-specific contrasts were smaller
(+0.189 and +0.083). Leave-one-donor-out estimates remained positive for both
states.

Within-donor injury-score associations were negative in C_TAL (-0.033;
p = 0.611) and M_TAL (-0.052; p = 0.392), but were not statistically
supported. Using neighborhood disease logFC as the exposure also produced a
negative C_TAL estimate (-0.074; p = 0.175) and a near-zero M_TAL estimate
(+0.014; p = 0.862). Thus, GSE211785 independently supported the donor-level
PAPPA2 association but did not establish the within-donor negative effect.

### Injured TAL states in GSE254185 reproduced analysis-unit divergence

GSE254185 provided a second independent test using cTAL1, cTAL2, mTAL,
TAL Injured, and TAL Inflammatory annotations. PAPPA2 donor pseudobulk was
higher in mTAL (+0.341; exact p = 0.0121; FDR = 0.0604 across the five TAL
states) and in TAL Injured (+0.168; p = 0.711). Neighborhood abundance was also
positive in TAL Injured. In contrast, the neighborhood-expression association
in TAL Injured was negative (rho = -0.096; FDR = 7.1e-05), and the
within-donor rank difference was negative using injury-marker exposure
(-0.290; exact p = 0.0702) and neighborhood disease logFC exposure
(-0.062; p = 0.825). The same qualitative positive-donor/negative-within-donor
pattern was therefore observed in a different TAL state, although the exact
aTAL2 and frTAL labels were not available (Additional file 14).

### Donor-stratified signature transfer recovered an aTAL2-like state

The frozen four-subtype marker panel assigned 86.66% of held-out KPMP cells at
the prespecified margin threshold. Among assigned cells, accuracy was 93.90%,
macro F1 was 93.17%, aTAL2 precision was 83.28%, and aTAL2 recall was 93.93%
(Additional files 15 and 16). The high held-out performance justified a
cautious projection to non-KPMP cohorts; it did not establish cross-cohort
label equivalence.

In GSE254185, the transfer identified 1,176 high-confidence aTAL2-like cells,
predominantly from TAL Inflammatory and TAL Injured states. Their donor-level
PAPPA2 difference was +0.554 (permutation p = 0.0043), whereas the
within-donor injury-score direction was negative (-0.581; p = 0.125). In
GSE211785, 3,751 high-confidence aTAL2-like cells were mostly derived from
M_TAL, with a smaller C_TAL contribution. The donor-level PAPPA2 difference
was +0.570 (p = 0.0235), and the within-donor neighborhood-logFC direction was
negative but not statistically supported (-0.095; p = 0.405)
(Additional file 17).

The external predictions therefore reproduced the donor-level direction and
the qualitative donor/within-donor conflict. They did not provide an exact
aTAL2 label, and the negative within-donor estimates were underpowered. The
transferred cells in GSE254185 were enriched for inflammatory and injured
states rather than a clean cortical aTAL2 population.

### Kidney eQTL supported a tubulointerstitial context but not disease chromatin

PAPPA2 was a strong tubulointerstitial eGene in NephQTL2
(gene-level FDR = 1.64e-14; model size 1.035), with no glomerular eQTL support
(FDR = 0.526) [16] (Additional file 18). This compartment
pattern is consistent with a tubule-centered signal. It is a regulatory
association and does not establish a disease effect.

PAPPA2 was not among NephQTL2's FDR-significant PTWAS associations for eGFR or
urinary albumin-to-creatinine ratio. This absence is reported as
non-significance rather than as a negative association, because instrument
availability and power can affect PTWAS inclusion. The fixed aTAL2 marker set
also showed no overall enrichment for tubulointerstitial eGenes: 11 of 35
markers were eGenes, with an exploratory odds ratio of 0.58 relative to the
NephQTL2 gene background. Individual markers had strong eQTLs, but the marker
program as a whole should not be described as enriched.

GSE211785 snATAC gene activity did not show higher disease accessibility for
PAPPA2. Using seven control and four disease donors, the disease-minus-control
difference was -0.029 in C_TAL (exact p = 0.843) and -0.016 in M_TAL
(exact p = 0.937; FDR = 1.0 for both) (Additional file 19). The same
cohort therefore supports the tubulointerstitial eQTL context but not a
disease-associated PAPPA2 chromatin increase. This is an orthogonal-modality
boundary result, not an independent-cohort replication.

### Candidate pathway results did not establish downstream mechanism

We did not find consistent support for a unified JAG1-NOTCH3-HES1 TAL program.
NOTCH3 was significant only in CKD frTAL in the donor pseudobulk analysis, and
HES1 showed mixed directions across subtypes. HES1 also had mixed positive and
negative GFR associations in Nephroseq.

Similarly, EPHB2 and STAT6 did not form a consistent TAL axis. STAT6 did not
show significant donor-level changes, and EPHB2 had subtype-dependent
directions. These results are retained as negative or exploratory findings and
do not support a causal pathway model.

## Discussion

### The main contribution is a unit-aware framework with empirical validation

The main contribution is a reusable unit-aware framework that makes the
analysis unit an explicit part of single-cell disease inference. The framework
does not introduce a new clustering algorithm. It combines donor pseudobulk,
neighborhood-expression association, within-donor cell rank, direction
standardization, composition-expression decomposition, and matched-null
calibration in one workflow. Applying it to the complete TAL resource showed
that the units were not interchangeable and that their non-equivalence was
observable at both the gene and state level. The result changes how a
disease-associated cell state should be
interpreted: donor, neighborhood, and within-donor analyses describe different
aspects of disease, so agreement between them cannot be assumed.

This interpretation matters because single-cell studies often report a
disease-associated state at one level of analysis and then generalize it to a
cell type or disease process. Our results show that such generalization is not
automatically justified. The central advance is therefore not a new molecular
mechanism, but a sharper framework for deciding what a single-cell disease
association does and does not represent.

### Matched random sets and effect-aware nulls answer different questions

The prespecified 400-gene panel addresses the principal generalizability
limitation of a single-gene study. Reversal between donor pseudobulk and
within-donor coupling was nearly 50% more frequent in biologically selected
TAL genes than in expression-matched controls (19.0% versus 12.7%). The result
persisted after expression stratification and was strongest for
PAPPA2-coexpressed genes and TAL subtype markers. Internal negative controls
were statistically indistinguishable from the external control pool.

That comparison nevertheless does not isolate gene-wise cross-unit coupling.
When the same target effects were re-paired under a conditional null, observed
donor-versus-within reversal was lower than expected from independent pairing
(19.0% versus 20.8%). Donor-versus-neighborhood reversal showed a small,
non-significant excess (21.9% versus 21.2%). The initial random-set enrichment
therefore demonstrated that selected genes had larger effects across units; it
did not demonstrate a panel-specific tendency toward opposite directions.

The neighborhood-versus-within pair retained a small conditional excess
(21.0% versus 19.4%), but this finding disappeared at higher effect thresholds
and was not significant in the biological 300-gene subset. The defensible
conclusion is graded: analysis-unit disagreement is common and can change
gene-level interpretation, whereas broad negative coupling between units is
not established by this panel.

### aTAL2 and frTAL show why analysis units cannot be merged

aTAL2 provides the strongest evidence that average disease elevation and
within-donor state coupling are distinct quantities. Donor pseudobulk and
neighborhood-expression association were positive, whereas within-donor
association was negative in combined disease and AKI and did not reach FDR in
CKD. This does not indicate a contradiction that must be forced into one
direction. Instead, it shows that the two units answer different questions:
one asks whether disease groups have higher average expression, and the other
asks whether neighborhood exposure tracks expression within an individual
donor.

The decomposition shows that the overall TAL PAPPA2 difference is not a
consequence of shifting subtype proportions. Across the four TAL subtypes, the
observed disease-reference difference was +0.563, the subtype-expression
component was +0.641, the composition component was +0.010, and the
interaction/covariance term was -0.088. AKI and CKD showed the same pattern.
Thus, higher subtype-specific expression between donors can coexist with a
negative within-donor neighborhood relationship. External transfer places
predicted aTAL2-like cells mainly in inflammatory and injured states, but this
cellular context does not explain the donor-level difference through subtype
proportions. The result does not establish cell-intrinsic induction; it shows
that a between-donor expression contrast and a within-donor context
association are distinct quantities.

frTAL shows that this phenomenon is not unique to aTAL2, although its effect is
weaker. The within-donor disease-minus-Reference direction was negative across
all three comparisons, but only the combined-disease comparison passed FDR.
The appropriate conclusion is therefore graded rather than binary: aTAL2
shows the strongest and most consistent discordance, whereas frTAL provides
secondary support for the same general problem.

These findings argue against treating one significant result as definitive.
When two defensible units disagree, the disagreement itself is informative
about the biological question being asked and about which inference is
supported by the data.

### Complete data and donor-aware design reduce major artefacts

Complete TAL coverage reduces the risk that the main conclusion is an artefact
of the derivative subset that was used in the earlier analysis. The complete
resource contained 106,851 TAL cells and differed substantially from the
4,446-cell matrix, especially for aTAL2 and frTAL. The subset and complete
pseudobulk analyses agreed in 48/60 directions, but 12 directions changed,
which shows that coverage can alter gene-level conclusions even when broad
trends remain. This motivates using the complete resource for all primary
inference.

The donor-aware design reduces, but does not eliminate, confounding by donor
composition and disease mix. Directed comparisons were supported by
age- and sex-adjusted models, leave-one-donor analysis, donor bootstrap,
cortex/medulla stratification, and neighborhood purity adjustment. These
checks support the main directions while preserving the aTAL2 contrast, but
they remain internal stability analyses; independent-cohort support is reported
separately below.

### Independent cohorts support donor-level association but not exact subtype replication

The two non-KPMP cohorts add evidence without erasing the boundary of the
study. In GSE211785, the donor-level PAPPA2 disease elevation was reproduced in
both C_TAL and M_TAL, and leave-one-donor-out estimates retained the positive
direction. In GSE254185, TAL Injured showed positive donor-level and
neighborhood-abundance directions but negative neighborhood-expression and
within-donor directions. This is the same qualitative problem that motivated
the discovery analysis: different defensible units answer different questions.

The independent data do not provide an exact replication of the aTAL2 and
frTAL annotations. GSE211785 contains C_TAL and M_TAL rather than the four
discovery subtypes, and GSE254185 contains injured and inflammatory TAL states
in an obstructive-nephropathy setting. The independent within-donor estimates
were also not statistically significant. The correct conclusion is therefore
that donor-level PAPPA2 disease elevation and the general analysis-unit problem
are independently supported, whereas the exact subtype-specific effect remains
to be tested in a cohort with matching annotations.

### Kidney eQTL provides context but does not convert the result into genetic validation

The NephQTL2 result adds an independent human kidney regulatory observation:
PAPPA2 is strongly regulated in tubulointerstitial tissue and is not a
glomerular eGene. This compartment pattern is concordant with the TAL focus and
supports the plausibility of a tubule-compartment signal. The result is
nevertheless a regulatory association in kidney tissue, not proof that PAPPA2
causes TAL injury, kidney dysfunction, or disease progression.

The boundary evidence is equally important. PAPPA2 was not among the
FDR-significant NephQTL2 PTWAS genes for eGFR or urinary
albumin-to-creatinine ratio, and the aTAL2 marker panel was not broadly
enriched for tubulointerstitial eGenes. The GSE211785 snATAC analysis also did
not show higher disease accessibility for PAPPA2. These findings prevent the
manuscript from describing the case gene as genetically validated. They show
instead that a strong kidney eQTL context can coexist with uncertainty about
the disease direction, cellular source, and causal role, which is itself
consistent with the analysis-unit problem.

### PAPPA2 is a case gene, not the mechanism

PAPPA2 should be interpreted as a case gene that illustrates the analysis-unit
problem, not as the central biological hypothesis. Its expression was
directionally positive across TAL subtypes and comparisons, and its positive
fraction was associated with lower eGFR in KPMP. This makes PAPPA2 useful for
testing whether a gene-state association survives a change in the analysis
unit.

The clinical interpretation of PAPPA2 nevertheless remains bounded. The KPMP
eGFR association is internal and based on binned eGFR categories. GSE137570
provided only a small and non-significant adjusted estimate. Nephroseq
provided multi-dataset directional support but consists of precomputed
analyses that may overlap, and GSE183456 belongs to the same KPMP consortium
as the discovery data. PAPPA2 is therefore a disease-associated TAL state
marker hypothesis, not a validated biomarker, therapeutic target, or causal
mediator.

### Negative pathway findings constrain the mechanistic interpretation

The negative pathway results prevent the analysis-unit finding from being
converted into an unsupported mechanism. JAG1-NOTCH3-HES1 and EPHB2-STAT6 did
not form consistent TAL programs across subtypes and comparisons. Reporting
these results is important because it shows that isolated correlations do not
establish an upstream regulator, downstream effector, or direct interacting
partner.

### Comparison with prior work

Prior human kidney atlases established extensive tubular heterogeneity and
disease-associated cell states [5,17]. The contribution of the present study
is not to rediscover TAL heterogeneity, but to show that conclusions about
that heterogeneity depend on the analysis unit. This extends atlas-based
descriptions from a question of which states exist to a question of whether a
state association survives a change in inference.

PAPPA2 has been studied mainly in the context of insulin-like growth factor
regulation and has not been established as a major determinant of TAL disease
[8,9,10]. Our analysis therefore places PAPPA2 within a subtype-resolved
analysis-unit framework rather than presenting it as a standalone mechanism.

### Data and methodological limitations

The generalizability of the findings is limited by the discovery resource and
by the reuse of the atlas-provided embedding and subtype labels. The analysis
was based on KPMP, and the UMAP coordinates were reused rather than
independently rebuilt. Milo used logMS normalization, and one AKI donor
without neighborhood representation was excluded from that comparison. These
choices can influence neighborhood boundaries and should be tested in
independent datasets, particularly with alternative graph construction,
normalization, and subtype annotation strategies.

External validation is now available but remains partial.
The discovery data and its spatial component (GEO superseries
GSE183279 [18]; Visium GSE183456 [19]) share KPMP consortium
provenance, so GSE183456 cannot be treated as an independent spatial
cohort. GSE211785 independently supported donor-level PAPPA2 disease
elevation, but its TAL annotations are C_TAL and M_TAL rather than aTAL2 and
frTAL, and its within-donor estimates were not significant. GSE254185 provided
an independent qualitative cross-unit signal in injured TAL, but its disease
context and annotations also differ. Donor-stratified signature transfer
supported a donor-level aTAL2-like direction in both cohorts, but the external
cells were not exact aTAL2 labels and the within-donor estimates were
underpowered. NephQTL2 provided strong tubulointerstitial eQTL support for
PAPPA2, but PAPPA2 was not among the FDR-significant PTWAS associations, and
the same-cohort snATAC analysis did not show higher disease accessibility.
Only 24 GSE137570 [20] donors had public
eGFR, and the adjusted result was not significant. Nephroseq analyses were
precomputed and may contain overlapping datasets. Sample-level GFR annotations
for GSE115098 were not publicly available, so that planned validation could not
be performed. Mendelian randomization was not feasible because no usable
PAPPA2 cis-pQTL was identified; this is an absence of a valid instrument rather
than a negative causal test. No perturbation, protein, or functional experiment
was performed.

### Implications and next steps

The practical implication is that single-cell disease studies should state the
analysis unit explicitly and test whether conclusions survive a change in that
unit. This recommendation is especially relevant for rare cell subtypes, where
small changes in donor or cell coverage can materially alter the result. The
next steps are to replicate the exact aTAL2 and frTAL annotations in an
independent kidney cohort, use continuous eGFR outcomes, and test the
PAPPA2-associated state with protein or functional experiments (Fig. 8).

[[FIGURE:Fig8]]

## Conclusions

We developed and validated a unit-aware framework showing that the analysis
unit is not a neutral choice in human kidney TAL single-cell data. Donor pseudobulk, Milo
neighborhood abundance, neighborhood-expression association, and within-donor
cell-rank association were not interchangeable: they yielded different
subtype-level conclusions, with aTAL2 as the strongest cross-unit discordance
example and frTAL as a weaker related pattern. Complete TAL coverage also
changed selected gene and subtype conclusions relative to the smaller
derivative matrix. A prespecified 400-gene panel showed that
donor-versus-within reversal was more frequent in biologically selected TAL
genes than in expression-matched controls. A conditional effect-pairing null
showed that this excess was explained by the selected effect-size
distribution, whereas a smaller neighborhood-versus-within excess remained
threshold sensitive. PAPPA2 acts as a case gene with directionally consistent
disease association across subtypes. Its donor-level disease elevation was
independently supported in GSE211785, and an injured-TAL cross-unit signal was
observed in GSE254185. Donor-stratified signature transfer reproduced the
donor-level direction of an aTAL2-like state but did not provide exact subtype
replication or a significant within-donor effect. NephQTL2 placed PAPPA2 in a
tubulointerstitial regulatory context, without glomerular eQTL support,
whereas the same-cohort snATAC analysis did not show a disease-associated
accessibility increase. The present evidence does not establish causality,
direct intercellular communication, genetic colocalization, or a defined
molecular mechanism. These findings support a bounded, human-data-based
conclusion: in renal single-cell studies, the analysis unit should be stated
explicitly, its influence tested before a disease-state association is
interpreted as a general biological property, and effect-aware calibration
should be reported alongside expression-matched null sets.

## List of abbreviations

AKI: acute kidney injury

AUC: area under the receiver operating characteristic curve

BH: Benjamini-Hochberg

CI: confidence interval

CKD: chronic kidney disease

C-TAL-A, C-TAL-B: KPMP human kidney atlas cortical thick ascending limb
subtype labels

C_TAL: cortical thick ascending limb label in GSE211785

DKD: diabetic kidney disease

eGFR: estimated glomerular filtration rate

eQTL: expression quantitative trait locus

FDR: false discovery rate

frTAL: KPMP human kidney atlas thick ascending limb subtype label

GEO: Gene Expression Omnibus

GFR: glomerular filtration rate

HKD: hypertensive kidney disease

KPMP: Kidney Precision Medicine Project

logFC: log2 fold change

LODO: leave-one-donor-out

M_TAL: medullary thick ascending limb label in GSE211785

MR: Mendelian randomization

MSigDB: Molecular Signatures Database

NephQTL2: human kidney eQTL resource

OR: odds ratio

PAPPA2: pappalysin 2

pQTL: protein quantitative trait locus

PTWAS: probabilistic transcriptome-wide association study

snATAC-seq: single-nucleus assay for transposase-accessible chromatin sequencing

snRNA-seq: single-nucleus RNA sequencing

SpatialFDR: spatial false discovery rate

TAL: thick ascending limb

UMAP: uniform manifold approximation and projection

## Declarations

### Ethics approval and consent to participate

This study used publicly available or controlled-access, de-identified human
kidney data. No new human samples were collected, and no participant
recruitment or intervention was performed. Ethical approval and informed
consent were not required for the present secondary analysis; data access and
use were governed by the terms of the original studies and repositories.

No animal experiments were performed in this study.

### Consent for publication

Not applicable. No identifiable individual-level data are presented.

### Clinical trial number

Clinical trial number: not applicable.

### Availability of data and materials

The datasets supporting the conclusions of this article are available in
the Gene Expression Omnibus repository, GSE183279
[18]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183279), which
includes single-nucleus RNA sequencing data (GSE183277) and spatial
transcriptomics data [19].
GSE137570 is available in the Gene Expression Omnibus repository
[20]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137570). The independent
validation datasets are GSE211785
[21,22]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785) and GSE254185
[23,24]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185).
NephQTL2 gene-level and PTWAS summary data were obtained from the
supplementary information of [16]. Nephroseq v5 was queried
through its public interface [13]
(https://www.nephroseq.org).

Access to the Gene Expression Omnibus datasets is unrestricted, including
for non-academic users. Nephroseq v5 is accessible through its public web
interface, and its use is subject to the terms of use published on the
Nephroseq website [13].

Project name: pappa2-tal-analysis-unit-sensitivity

Project home page:
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity

Archived version: pending release of the v0.2.0 analysis archive

Operating system(s): Platform independent

Programming language: Python; R

Other requirements: see requirements.txt in the project repository

License: MIT for the unitaware software; project data follow the repository
data license

Any restrictions to use by non-academics: none

Unit-aware software package:

Project name: unitaware v0.2.0

Project home page:
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity/tree/main/unitaware

Archived version: included in the archived analysis release cited above

Operating system(s): Platform independent

Programming language: Python

Other requirements: Python 3.10 or higher; numpy, pandas, scipy

License: MIT

Any restrictions to use by non-academics: none

The processed analysis tables supporting the conclusions of this article are
included in Additional files 1-25 and 28-29.

Analysis code is available at
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity. The
versioned archive containing the signature-transfer, regulatory-support,
composition-expression, and prespecified analysis-unit panel is available at
the project repository. The final exact-version DOI will be inserted after the
corrected v0.2.0 archive is released.

### Competing interests

The authors declare no competing interests.

### Funding

This research received no specific grant from any funding agency in the
public, commercial, or not-for-profit sectors.

### Authors' contributions

P.Z. conceived and designed the study, performed the formal analysis and
data curation, and wrote the first draft. X.Y. contributed to validation
and data curation, and revised the manuscript. C.Z. contributed to the
software and data curation, and revised the manuscript. All authors read
and approved the final manuscript.

## Additional files

Additional file 1. File format: Comma-separated values (.csv). Title:
Permutation null distributions for the within-donor
cell-rank analysis. Description: Null distributions from 1,000 within-donor
expression permutations for each gene, subtype, and comparison.

Additional file 2. File format: Comma-separated values (.csv). Title:
Complete-data Milo neighborhood summary. Description: Neighborhood counts,
directions, and median logFC for each subtype and disease comparison in the
complete 106,851-cell TAL analysis.

Additional file 3. File format: Comma-separated values (.csv). Title:
Subset-level Augur separability summary. Description: Mean AUC values and
subsample intervals for disease versus Reference separability in the
4,446-cell derivative subset.

Additional file 4. File format: Comma-separated values (.csv). Title:
Complete versus subset pseudobulk direction comparison. Description:
Direction of effect for each gene, subtype, and comparison in the complete
analysis and in the 4,446-cell subset, including the 12 comparisons with
discordant direction.

Additional file 5. File format: Comma-separated values (.csv). Title:
Neighborhood-expression association robustness. Description: Spearman
correlations between neighborhood expression and Milo logFC, with
leave-one-donor and donor-stratified bootstrap stability.

Additional file 6. File format: Comma-separated values (.csv). Title:
Donor-stratified bootstrap replicates for neighborhood
association. Description: Per-replicate Spearman correlations used to derive
the donor-stratified bootstrap intervals.

Additional file 7. File format: Comma-separated values (.csv). Title: Joint
analysis-unit concordance table. Description: Donor pseudobulk, neighborhood
association, and within-donor cell-rank directions joined for each gene,
subtype, and comparison.

Additional file 8. File format: Comma-separated values (.csv). Title:
Within-donor cell-rank association results. Description: Donor-level
Spearman correlations between gene expression and neighborhood exposure for
each gene, subtype, and comparison.

Additional file 9. File format: Comma-separated values (.csv). Title:
Complete donor pseudobulk results. Description: Donor pseudobulk effect
sizes, age- and sex-adjusted estimates, and false discovery rates for the
five analyzed genes.

Additional file 10. File format: Comma-separated values (.csv). Title:
Leave-one-donor-out replicates for neighborhood
association. Description: Per-replicate Spearman correlations computed with
one donor removed at a time.

Additional file 11. File format: Comma-separated values (.csv). Title: KPMP
ordinal eGFR models. Description: Ordinal logistic model estimates for the
association between PAPPA2-positive TAL fraction and baseline eGFR category.

Additional file 12. File format: Comma-separated values (.csv). Title:
GSE137570 directional eGFR association. Description: Unadjusted and adjusted
PAPPA2-eGFR association estimates for the 24 donors with public eGFR.

Additional file 13. File format: Comma-separated values (.csv). Title:
Nephroseq precomputed GFR associations. Description: Precomputed GFR
association results for PAPPA2 and four control genes, summarized by
dataset, tissue type, sample count, correlation coefficient, and p value.

Additional file 14. File format: Comma-separated values (.csv). Title:
Independent non-KPMP validation of PAPPA2 and analysis-unit direction.
Description: TAL state, donor comparison, effect size, permutation p value,
false discovery rate, direction, and number of donors for GSE211785 and
GSE254185.

Additional file 15. File format: Comma-separated values (.csv). Title:
Donor-stratified KPMP TAL marker panel. Description: Top 35 training-donor
markers for each of the four TAL subtypes, with mean expression, detection
fraction, and specificity.

Additional file 16. File format: Comma-separated values (.csv). Title:
Held-out KPMP signature-transfer validation. Description: Held-out prediction
coverage, accuracy, macro F1, precision, recall, and confusion matrices at the
prespecified marker-margin thresholds.

Additional file 17. File format: Comma-separated values (.csv). Title:
External aTAL2-like signature transfer. Description: High-confidence predicted
cells, backbone-state composition, and donor-level and within-donor PAPPA2
results for GSE254185 and GSE211785.

Additional file 18. File format: Comma-separated values (.csv). Title:
NephQTL2 gene-level eQTL results for PAPPA2 and fixed markers. Description:
Tubulointerstitial and glomerular gene-level FDR, eGene status, model size, and
prespecified marker annotations from [16].

Additional file 19. File format: Comma-separated values (.csv). Title:
GSE211785 snATAC gene-activity analysis. Description: Donor pseudobulk values,
exact permutation P values, leave-one-donor ranges, and false discovery rates
for PAPPA2 and fixed marker genes in C_TAL and M_TAL.

Additional file 20. File format: Comma-separated values (.csv). Title:
Composition and expression decomposition. Description: Overall TAL
contribution, subtype-expression component, subtype-composition component,
interaction/covariance remainder, donor-bootstrap intervals, and subtype-level
contributions for PAPPA2 and control genes.

Additional file 21. File format: Comma-separated values (.csv). Title:
Prespecified 400-gene panel and matched control pool. Description: Group
assignment, selection rank, expression, and detection statistics for the
target panel and the 2,000-gene external negative-control pool.

Additional file 22. File format: Comma-separated values (.csv). Title:
Gene-level analysis-unit metrics. Description: Donor pseudobulk,
neighborhood-expression, and within-donor effects with validity flags and
standardized reversal indicators for all panel and control genes.

Additional file 23. File format: Comma-separated values (.csv). Title:
Random-set analysis-unit sensitivity results. Description: Gene-level reversal
fractions for 1,000 random 400-gene and 300-gene sets from the expression-matched
negative-control pool, including expression-stratified sensitivity sets.

Additional file 24. File format: Comma-separated values (.csv). Title:
Unit-aware simulation runs. Description: Per-replicate donor, neighborhood, and
within-donor effects and reversal fractions across six simulation scenarios.

Additional file 25. File format: Comma-separated values (.csv). Title:
Unit-aware simulation summary. Description: Mean effect sizes and reversal
fractions by simulation scenario.

Additional file 26. File format: Tagged Image File Format (.tif). Title:
Analysis-unit robustness. Description: a, Expression-matched random-set null
for donor-within reversal with the 400-gene target and 300-gene biological
panel values. b, Observed reversal fraction and effect-aware conditional null
mean with 95% interval for the three gene-specific unit pairs. c, Target and
expression-matched reversal fractions across standardized-effect thresholds
from 0.5 to 2.0.

Additional file 27. File format: Tagged Image File Format (.tif). Title:
Unit-aware simulation benchmark. Description: a, Mean donor and within-donor
effects across six simulation scenarios. b, Donor-within reversal fraction by
scenario. The opposing scenarios separate from concordant and largely null
scenarios.

Additional file 28. File format: Comma-separated values (.csv). Title:
Analysis-unit threshold sensitivity. Description: Target-panel and
expression-matched random-set reversal fractions for three unit pairs at
standardized-effect thresholds from 0.5 to 2.0, including the 300 biological
genes and expression-stratified sets.

Additional file 29. File format: Microsoft Excel Worksheet (.xlsx). Title:
Effect-aware pairing permutation. Description: Summary and 10,000-run null
distributions for donor versus within-donor, donor versus
neighborhood-expression, and neighborhood-expression versus within-donor
reversal after preserving effect-size marginals and breaking gene-wise unit
pairing.

## References

1. Bikbov B, Purcell CA, Levey AS, Smith M, Abdoli A, Abebe M, et al. Global, regional, and national burden of chronic kidney disease, 1990-2017: a systematic analysis for the Global Burden of Disease Study 2017. Lancet. 2020;395(10225):709-33. doi:10.1016/S0140-6736(20)30045-3.

2. Hoste EAJ, Kellum JA, Selby NM, Zarbock A, Palevsky PM, Bagshaw SM, et al. Global epidemiology and outcomes of acute kidney injury. Nat Rev Nephrol. 2018;14(10):607-25. doi:10.1038/s41581-018-0052-0.

3. Ferenbach DA, Bonventre JV. Mechanisms of maladaptive repair after AKI leading to accelerated kidney ageing and CKD. Nat Rev Nephrol. 2015;11(5):264-76. doi:10.1038/nrneph.2015.3.

4. Bonventre JV, Yang L. Cellular pathophysiology of ischemic acute kidney injury. J Clin Invest. 2011;121(11):4210-21. doi:10.1172/JCI45161.

5. Lake BB, Menon R, Winfree S, Hu Q, Melo Ferreira R, Kalhor K, et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature. 2023;619(7970):585-94. doi:10.1038/s41586-023-05769-3.

6. Knepper MA, Kim GH, Fernández-Llama P, Ecelbarger CA. Regulation of thick ascending limb transport by vasopressin. J Am Soc Nephrol. 1999;10(3):628-34. doi:10.1681/ASN.V103628.

7. Brezis M, Rosen S, Silva P, Epstein FH. Selective vulnerability of the medullary thick ascending limb to anoxia in the isolated perfused rat kidney. J Clin Invest. 1984;73(1):182-90. doi:10.1172/JCI111189.

8. Oxvig C. The role of PAPP-A in the IGF system: location, location, location. J Cell Commun Signal. 2015;9(2):177-87. doi:10.1007/s12079-015-0259-9.

9. Dauber A, Muñoz-Calvo MT, Barrios V, Domené HM, Kloverpris S, Serra-Juhé C, et al. Mutations in pregnancy-associated plasma protein A2 cause short stature due to low IGF-I availability. EMBO Mol Med. 2016;8(4):363-74. doi:10.15252/emmm.201506106.

10. Conover CA, Boldt HB, Bale LK, Clifton KB, Grell JA, Mader JR, et al. Pregnancy-associated plasma protein-A2 (PAPP-A2): tissue expression and biological consequences of gene knockout in mice. Endocrinology. 2011;152(7):2837-44. doi:10.1210/en.2011-0036.

11. Asghari M, Sabo AR, Barwinska D, Melo Ferreira R, Ferkowicz MJ, Bowen WS, et al. Integration of spatial protein imaging and transcriptomics in the human kidney tracks the regenerative potential of proximal tubules. Sci Adv. 2025;11(33):eadv8918. doi:10.1126/sciadv.adv8918.

12. Kennedy C, Doyle R, Gough O, McEvoy C, McAnallen S, Hughes M, et al. A novel role for FERM domain-containing protein 3 in CKD. Kidney360. 2024;5(12):1799-812. doi:10.34067/KID.0000000602.

13. Nephroseq v5. University of Michigan, Ann Arbor, MI, USA. https://www.nephroseq.org. Accessed 29 September 2026.

14. Skinnider MA, Squair JW, Kathe C, Anderson MA, Gautier M, Matson KJE, et al. Cell type prioritization in single-cell data. Nat Biotechnol. 2021;39(1):30-4. doi:10.1038/s41587-020-0605-1.

15. Dann E, Henderson NC, Teichmann SA, Morgan MD, Marioni JC. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nat Biotechnol. 2022;40(2):245-53. doi:10.1038/s41587-021-01033-z.

16. Han SK, McNulty MT, Benway CJ, Wen P, Greenberg A, Onuchic-Whitford AC, et al. Mapping genomic regulation of kidney disease and traits through high-resolution and interpretable eQTLs. Nat Commun. 2023;14:2229. https://doi.org/10.1038/s41467-023-37691-7.

17. Cippà PE, Sun B, Liu J, Chen L, Naesens M, McMahon AP. Transcriptional trajectories of human kidney injury progression. JCI Insight. 2018;3(22):e123151. doi:10.1172/jci.insight.123151.

18. Kidney Precision Medicine Project. Single-nucleus RNA sequencing resource, GEO superseries GSE183279. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183279. Accessed 29 September 2026.

19. Kidney Precision Medicine Project. Spatial transcriptomics resource, GEO GSE183456. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183456. Accessed 29 September 2026.

20. Human kidney transcriptome cohort, GEO GSE137570. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137570. Accessed 29 September 2026.

21. Abedini A, Kloetzer K, Ma Z, Frederick J, Dhillon P, Balzer M, et al. Single-cell multi-omic and spatial profiling of human kidneys implicates the fibrotic microenvironment in kidney disease progression. Nat Genet. 2024;56(8):1712-24. doi:10.1038/s41588-024-01802-x.

22. Human kidney multi-omic atlas, GEO GSE211785. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785. Accessed 28 October 2026.

23. Reck M, Baird DP, Veizades S, Sutherland C, Bell RMB, Hur H, et al. Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype. Nat Commun. 2025;16(1):4745. doi:10.1038/s41467-025-59997-4.

24. Human kidney multiome resource with obstructive nephropathy, GEO GSE254185. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185. Accessed 28 October 2026.

## Figure legends

**Fig. 1. TAL location and analysis-unit study design.**
The left panel shows a nephron with the thick ascending limb (TAL)
highlighted across the cortex and medulla. The right panels summarize the
analysis. The discovery resource was the complete KPMP single-nucleus TAL
set, comprising 106,851 cells annotated as C-TAL-A, C-TAL-B, aTAL2 or
frTAL. A previously generated 4,446-cell derivative matrix was retained only
as a coverage and directional comparator. Four established units of
inference were compared: donor pseudobulk, neighborhood abundance,
neighborhood-expression association, and within-donor cell-rank association.
Donor-level robustness analyses (leave-one-donor refitting, donor-stratified
bootstrap, age- and sex-adjusted models, and cortex or medulla
stratification) were treated as stability checks rather than as separate
units of inference. PAPPA2 was analyzed as a pre-specified case gene. No
arrow in this schematic implies causality.

**Fig. 2. Subtype-resolved disease-related neighborhood remodeling.**
Complete-data Milo results for the four TAL subtypes across three
comparisons (combined disease, AKI, and CKD versus Reference). a Fraction of
annotated neighborhoods that were significant. b Median neighborhood logFC.
c Number of significant neighborhoods. Neighborhoods were built from
atlas-provided UMAP coordinates and tested with Milo using logMS
normalization and a fixed seed; the significance threshold was SpatialFDR <
0.05. One AKI donor without neighborhood representation was excluded from
the AKI comparison. Numbers of annotated neighborhoods differed between
comparisons; in the combined-disease comparison they were 4,012 (C-TAL-A),
2,805 (C-TAL-B), 1,465 (aTAL2), and 1,332 (frTAL). Panels show counts and
summary statistics without error bars.

**Fig. 3. Analysis-unit comparison design.**
Schematic definition of the four analysis units compared in this study, with
the question each unit addresses and the limitation that applies to it.
Donor pseudobulk asks whether mean expression differs between disease and
Reference donors and is limited by donor composition and disease mix.
Neighborhood-expression association asks whether neighborhoods with higher
gene expression also have higher disease logFC and is limited because Milo
neighborhoods overlap. Within-donor cell-rank association asks whether cells
with greater neighborhood exposure have higher expression inside the same
donor and describes a state association rather than causality. Donor
robustness analysis resamples donors to test dependence on individual
donors or tissue region and provides internal stability only. The schematic
illustrates that one gene can be positive in one unit and negative in
another; no data are plotted.

**Fig. 4. Complete TAL coverage and subset discordance.**
Comparison of the complete 106,851-cell TAL source with the 4,446-cell
derivative subset. a Cells per subtype in the complete source and in the
subset (log scale). b Number of Reference donors per subtype in both
sources. c Gene, subtype, and comparison tests in which the complete and
subset pseudobulk analyses agreed or disagreed in direction; 60 tests were
evaluated (five genes × four subtypes × three comparisons). d Discordant
tests by gene. The coverage difference was largest for the sparse subtypes:
aTAL2 contained 6,624 Reference cells from 83 donors in the complete source
versus 102 cells from 38 donors in the subset, and frTAL contained 6,656
Reference cells from 86 donors versus 142 cells from 48 donors. Bars show
counts; no error bars are shown.

**Fig. 5. PAPPA2 case analysis across TAL subtypes.**
PAPPA2 results across the four analysis units, shown for three comparisons
(rows) and four TAL subtypes (columns). a Donor pseudobulk mean expression
difference between disease and Reference donors. b Age- and sex-adjusted
linear model beta. c Neighborhood-expression association, expressed as the
Spearman correlation between neighborhood mean PAPPA2 expression and Milo
logFC across neighborhoods. d Within-donor cell-rank association, expressed
as the median Spearman correlation between PAPPA2 expression and
neighborhood exposure within each donor, using donors with at least five
cells. Values are printed in each cell. Panels a and b are donor-level
estimates; panels c and d use different units of inference and should not be
combined.

**Fig. 6. aTAL2 and frTAL cross-unit discordance.**
Cross-unit comparison for the two subtypes that showed direction conflicts.
a aTAL2 donor pseudobulk mean expression difference for the three
comparisons. b aTAL2 median within-donor cell-rank correlation for the same
comparisons. c Number of comparisons with a direction conflict between the
donor pseudobulk and within-donor cell-rank units, by subtype. In aTAL2,
donor-level association was positive while the within-donor association was
negative in combined disease and AKI (FDR = 0.0086 for each); the CKD
difference was also negative but did not pass FDR (0.224). frTAL showed the
same qualitative pattern with weaker support (combined disease FDR = 0.018).
Bars show point estimates; no error bars are shown.

**Fig. 7. Internal robustness of complete TAL pseudobulk effects.**
Internal stability analysis of the PAPPA2 donor pseudobulk effect across the
12 gene, subtype, and comparison groups (four subtypes × three comparisons).
a Full-data effect (points) together with the range obtained when each donor
was removed in turn (whiskers span the minimum to the maximum of the
leave-one-donor refits). b Donor-stratified bootstrap median with the 95%
interval from 1,000 resamples. c Cortex and medulla mean
disease-minus-Reference difference for the combined-disease comparison. The
PAPPA2 direction was preserved in all 12 leave-one-donor groups, and cortex
and medulla agreed in direction in 11 of 12 groups. These intervals describe
donor-level stability and are not independent-cell significance tests.

**Fig. 8. Subtype-specific TAL state model.**
Schematic summary of the interpretation supported by this study. The upper
part shows subtype-specific TAL neighborhood states (C-TAL-A, C-TAL-B, and
frTAL with positive neighborhood association, and aTAL2 as an exploratory
decoupling), the donor-aware framework used here (Milo, Augur, donor
pseudobulk, and within-donor cell rank), and the PAPPA2-associated TAL state
together with lower kidney function. The lower panels state what the data
support (subtype-specific state structure, exploratory aTAL2 decoupling, and
a PAPPA2 association) and what they do not support (causality, direct
intercellular communication, and a unified JAG1-NOTCH3-HES1 or EPHB2-STAT6
mechanism). No arrow implies causal direction.
