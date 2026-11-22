# Unit-aware calibration of single-cell disease associations: analysis-unit, effect-size, and graph-stability audits in human kidney atlases

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

**Background:** Single-cell kidney studies routinely infer disease states from
donors, neighborhoods, or individual cells. These units answer different
questions but are often treated as interchangeable. Expression-matched nulls can
mistake effect-size selection for biological discordance; neighborhood graphs
can introduce rarely audited instability.

**Methods:** unitaware specifies donor pseudobulk, Milo neighborhood abundance,
neighborhood-expression association, and within-donor cell rank as separate
estimands. It applies effect-size-aware conditional nulls and audits calibration
across graph resolutions, neighborhood proportions, seeds, and embeddings. We
applied the workflow to 106,851 thick ascending limb (TAL) cells from 223
Kidney Precision Medicine Project (KPMP) donors and a prespecified 400-gene
panel. We then evaluated three kidney cohorts, matched non-kidney cohorts, and
PAPPA2 as a benchmark.

**Results:** Expression matching suggested excess donor-within reversal (19.0%
versus 12.7%). The donor-conditioned null removed this excess in
all three kidney cohorts (0/3 significant). Calibrated panel residuals varied
across graph resolutions, proportions, seeds, and embeddings. aTAL2/PAPPA2
remained negative across UMAP settings but attenuated in panel PCA30, despite
positive donor pseudobulk. PAPPA2 donor elevation was positive in all 12
subtype comparisons and was not explained by composition. PAPPA2-positive TAL
fraction tracked lower estimated glomerular filtration rate (eGFR) in KPMP. A
frozen 11-gene PAPPA2-associated TAL-injury score tracked lower eGFR in an
independent cohort (n = 87), while control-gene comparisons did not establish
gene-specific superiority.

**Conclusions:** The advance is a unit-aware calibration standard that
separates genuine cross-unit signal from effect-size selection and graph- or
embedding-dependent inference. unitaware makes this audit reproducible across
cohorts; PAPPA2 remains a benchmark, not causal or clinical evidence.

**Keywords:** single-cell genomics, analysis unit, calibration stability,
kidney disease, thick ascending limb, donor-aware inference, software framework

## Introduction

Chronic kidney disease and acute kidney injury remain major clinical burdens
[1,2], and maladaptive tubular repair links acute injury to progressive kidney
disease [3,4]. Human kidney atlases now resolve injury, adaptation, and repair
states across nephron segments [5]. These atlases place the central inferential
question at the level of the analysis unit rather than the cell alone. The same
cells can support a between-donor expression contrast, a neighborhood abundance
contrast, a neighborhood-expression association, or a within-donor cell-rank
estimate. These summaries represent different biological and statistical
quantities, and the same atlas can produce different answers without changing a
single cell. Their agreement is an empirical finding, not a property
guaranteed by shared data. Current workflows often assemble such
summaries as convergent evidence for one disease effect.

Current single-cell disease studies commonly combine donor pseudobulk, Milo
neighborhood abundance, neighborhood-expression association, and within-donor
cell rank as complementary views of one disease state. Donor pseudobulk
estimates a between-donor difference in mean expression. Milo neighborhood
abundance estimates whether tissue neighborhoods are more abundant in disease.
Neighborhood-expression association tests whether high-expression neighborhoods
also carry higher disease log2 fold change (logFC). Within-donor cell-rank
association tests whether neighborhood exposure tracks expression inside the
same donor. This design is attractive because each unit uses the same cells and
can appear to broaden the evidence base. The implicit assumption is that
agreement across units validates a shared disease effect.

That assumption fails because the units differ in replication, exposure, and
standardization, and because gene selection can create cross-unit effects
without biological coupling. Donor pseudobulk is dominated by between-donor
variation, whereas neighborhood abundance reflects tissue composition.
Neighborhood-expression association depends on graph-defined exposure, whereas
within-donor rank removes donor-level differences and asks a local question. A positive result in
one unit does not imply a positive result in another. Selection adds a
further failure mode: expression-matched random sets preserve expression and
detection but leave the joint distribution of effect sizes unconstrained. If
target genes were selected for large effects in one unit, their effects in
another may regress or reverse even without biological coupling.

The unresolved barrier is to distinguish biological cross-unit discordance from
effect-size selection while testing whether the remaining signal survives graph
construction. Expression matching alone cannot make this distinction because it
does not preserve the relevant marginal effect distributions. A valid
comparison must condition on effect size and then ask whether the observed panel
departs from that conditional expectation. Neighborhood analyses add another
requirement because Milo abundance and association depend on graph resolution,
neighborhood proportion, and sampling seed. External cohorts add another
requirement because different annotations and approximate exposures can support
portability without reproducing the same estimand. These constraints require an
audit that reports the unit, conditioning direction, and graph or embedding
stability together.

We developed and applied a unit-aware audit to the complete KPMP thick
ascending limb (TAL) resource, comprising 106,851 cells from 223 donors across
four annotated subtypes. TAL cells provide a sensitive disease model because
active sodium, potassium, and chloride transport makes them metabolically
demanding and vulnerable to ischemic and inflammatory injury [6,7]. We
prespecified a 400-gene panel and evaluated effect-aware calibration across
three kidney cohorts, matched non-kidney boundaries, Milo graph settings,
neighborhood proportions, seeds, and embeddings. We selected PAPPA2 as a
benchmark case because it has established systemic insulin-like growth factor
regulation but a less defined renal state role [8,9,10]. PAPPA2 served as a
prespecified benchmark for testing whether the workflow could retain a
cross-modally supported signal. The resulting workflow makes donor,
neighborhood, and within-donor analyses separate estimands, evaluates both
conditioning directions, and treats graph and embedding sensitivity as part of
the analysis specification. Together, these elements define a reusable workflow
for deciding which single-cell disease associations merit interpretation at a
given inferential scale.

## Methods

### Study design and data sources

The study developed and evaluated a unit-aware workflow to determine whether
subtype-level disease conclusions remained stable when the unit of inference
changed (Fig. 1). Public or controlled-access human kidney data were used, and
no new human samples were collected. No new ethical approval was required for
secondary analysis of de-identified data. The study was designed as a
framework evaluation rather than a search for a new molecular mechanism.

[[FIGURE:Fig1]]

Data sources were selected to provide one discovery resource and separate
directional evaluations. The discovery resource was the KPMP human kidney
single-nucleus RNA sequencing (snRNA-seq) dataset [5], which contained 1,388,643 cells and
35,455 features. We extracted all cells annotated as C-TAL-A, C-TAL-B, aTAL2,
or frTAL, yielding 106,851 TAL cells. A previously generated 4,446-cell
derivative matrix was retained only as a coverage and directional comparator.
Spatial cross-modality analysis used the KPMP GSE183456 Visium resource, part
of the KPMP Gene Expression Omnibus (GEO) superseries GSE183279, which contained 23 sections from 22
donors [5,11]. Because GSE183456 shares KPMP consortium provenance with the
discovery data, it was treated as same-consortium support rather than an
external cohort. Separate directional estimated glomerular filtration rate
(eGFR) testing used GSE137570 [12], and multi-dataset glomerular filtration
rate (GFR) association testing used Nephroseq v5 [13]. Independent
directional evaluation used the non-KPMP GSE211785 multi-omic kidney atlas and
the independent GSE254185 kidney multiome resource. These two non-KPMP datasets
were used to evaluate portability of the calibration framework and principal
directions, not to validate a specific disease mechanism. GSE284715/GSE284742
were not analyzed as quantitative validation cohorts. Their independently
published `TAL_Injury_h1` marker list was used only to define a frozen clinical
score in GSE175759, as described below.

### TAL subtypes and the four analysis units

We defined four TAL subtypes and four analysis units before any disease
comparison, so that each estimand could be interpreted separately. The subtypes
were C-TAL-A, C-TAL-B, aTAL2, and frTAL. All primary analyses used the complete
106,851-cell TAL set and the atlas-provided uniform manifold approximation and
projection (UMAP) coordinates.

We compared the four analysis units because each addresses a distinct question
about disease-related TAL state. Donor pseudobulk asks whether mean expression
differs between disease and Reference (control) donors. Milo neighborhood
abundance asks whether neighborhood abundance differs between disease and
Reference. Neighborhood-expression association asks whether neighborhoods with
higher expression also have higher disease logFC. Within-donor cell-rank
association asks whether cells with greater neighborhood exposure have higher
expression inside each donor. Composition-expression decomposition was included
as a descriptive extension of donor pseudobulk, not as a fifth gene-specific
analysis unit. These estimates reuse the same cells and are not independent
validations of one another.

Donor-level robustness was treated as a stability check, not as a fifth analysis
unit, to avoid interpreting resampling stability as an independent biological
measurement. The primary Milo threshold was spatial false discovery rate
(SpatialFDR) < 0.05. Milo used logMS
normalization and a fixed random seed of 20260927. One AKI donor with zero
neighborhood counts was excluded from that comparison. Augur was used only as a
secondary subset-level analysis of disease versus Reference separability [14],
with 20 subsamples, 30 cells per subsample, and 3-fold cross-validation.

### Prespecified 400-gene panel

The prespecified 400-gene panel was frozen before outcome analysis to define the
biological and control gene sets used for panel-level calibration. To prevent
information leakage during gene selection, the 223 KPMP donors were randomly
divided into 156 training and 67 held-out donors using a fixed seed of 20261028.
PAPPA2 coexpression genes were selected only in the training donors; held-out
donors were not used for gene selection. The panel contained
100 TAL subtype markers and 100 Molecular Signatures Database (MSigDB)
Hallmark genes from injury, inflammatory,
hypoxic, or metabolic pathways. It also contained 100 genes coexpressed with
PAPPA2 in the 156 training donors and 100 internal expression-matched negative
controls. A separate 2,000-gene pool was matched to the target mean expression
and detection distribution in ten strata. From this pool, we drew 1,000 random
400-gene sets without replacement within each set.

Only the three gene-specific units were compared across this panel: donor
pseudobulk, neighborhood-expression association, and within-donor cell-rank
association. Donor-versus-within reversal was the primary endpoint;
donor-versus-neighborhood and neighborhood-versus-within reversal were
secondary.

### Effect-aware calibration and null-model construction

The analysis-unit audit tested whether donor-versus-within reversal was more
frequent among biologically selected TAL genes than among expression-matched
controls. It also tested whether any excess persisted after conditioning on the
observed effect-size marginals. For each gene, subtype, and disease comparison, we
standardized the donor effect, neighborhood-expression effect, and within-donor
cell-rank effect using the pooled target and negative-control gene distribution.
The primary donor-within reversal indicator was defined as:

`R_i = 1` if `sign(zD_i) × sign(zW_i) < 0` and
`min(|zD_i|, |zW_i|) >= 0.5`; otherwise `R_i = 0`.

The 0.5 threshold was fixed a priori as a moderate standardized-effect boundary.
Threshold sensitivity was evaluated from 0.5 to 2.0.

We averaged `R_i` across evaluable subtype-comparison strata for each gene, then
averaged the gene-level fractions across the prespecified panel. The
expression-matched null used 1,000 random 400-gene sets drawn without
replacement within ten expression and detection strata from the 2,000-gene
control pool.

To condition on effect-size selection, we held each target gene's left-unit
effect fixed and sampled right-unit effects without replacement from the
same-stratum right-effect pool. This procedure was repeated 10,000 times within
each subtype and comparison. Every permutation preserved the left
effect, the observed marginal distribution of right-unit effects, and all valid
strata while breaking only the gene-wise pairing between units. We recomputed
`R_i`, the gene-level reversal fraction, and the panel statistic in every
permutation. One-sided empirical P values were calculated as
`(1 + number of null statistics >= observed) / (permutations + 1)`.

The expression-matched comparison was prespecified, whereas the effect-aware
permutation was a data-adaptive secondary calibration. Because this null
conditions on the donor effect, we repeated it with the donor and within-donor
units exchanged to assess order dependence. Paired outputs from both
conditioning directions are provided in Additional file 1. The panel,
threshold, units, subtypes, and comparisons were frozen before both analyses.
The 400-gene panel
was used for panel-level calibration of cross-unit reversal. A targeted five-gene
set comprising PAPPA2 and four control genes was used for detailed subtype and
comparison estimates across 60 tests.

### Donor-level pseudobulk expression

Donor pseudobulk served as the reference analysis unit because donor-level mean
expression is the most common way to relate cell-state genes to disease.
For the candidate-gene analysis, raw counts for PAPPA2, CXCL12, NOTCH3, HES1,
and EPHB2 were extracted directly from the sparse source matrix. Counts were
normalized using the full per-cell
RNA count and transformed with log1p. Expression was averaged within each
donor and TAL subtype. Disease groups were compared with Reference donors using
two-sided Mann-Whitney tests for unadjusted group contrasts. Donor-level linear
models adjusted for age midpoint and sex provided the covariate-adjusted
analyses. Both unadjusted and covariate-adjusted estimates were prespecified.
Benjamini-Hochberg correction was applied across
the 60 gene, subtype, and comparison tests. The 4,446-cell subset result was
retained only as a comparator, and its original 120-test correction was not
mixed with the complete-TAL analysis.

### Composition and expression decomposition

We used composition-expression decomposition to determine how much of the
overall four-subtype TAL PAPPA2 difference reflected subtype composition versus
subtype-specific expression. For each donor, the
overall TAL mean was the cell-weighted sum of subtype proportions and
subtype-specific mean expression. Group differences were separated into
composition at reference means, expression at reference proportions, and an
exact remainder containing interaction and within-group covariance terms. The
decomposition was:

`ΔĒ = Σ_k (ΔP_k × E_{k,ref}) + Σ_k (P_{k,ref} × ΔE_k) + remainder`,

where `k` indexes TAL subtype. `ΔP_k` is the disease-minus-Reference subtype
proportion, and `ΔE_k` is the disease-minus-Reference subtype mean expression.
Missing subtype means were set to zero contribution when the corresponding
proportion was zero. Confidence intervals were obtained from 2,000 donor-level
bootstrap resamples within disease and Reference groups.

### Milo neighborhood abundance

Milo neighborhood abundance tested whether disease was associated with shifts
in local TAL neighborhoods rather than only with donor-level average
expression. Neighborhoods were rebuilt from the complete TAL set using the
atlas UMAP coordinates, k = 30, a neighbor ratio of 0.1, and refined sampling.
Disease combined, AKI, and CKD were tested against Reference [15]. Results were
interpreted as differences in TAL-associated neighborhood state, not as evidence
of direct cell-cell contact.

### Neighborhood-expression association

Neighborhood-expression association tested whether neighborhoods with higher
gene expression were also those with higher disease logFC. Each Milo
neighborhood was represented by its mean PAPPA2 expression and its Milo logFC,
and Spearman correlation was calculated across neighborhoods for each
comparison and subtype. Stability was assessed with 500 donor-stratified
resamples, leave-one-donor-out analysis, and a partial Spearman correlation
controlling for neighborhood purity. These resampling intervals describe
donor-level stability and were not interpreted as independent cell-level P
values.

### Within-donor cell-rank association

Within-donor cell-rank association tested whether cells with greater
neighborhood exposure had higher expression inside the same donor. For each
cell, neighborhood exposure was defined as the mean Milo logFC of the
neighborhoods assigned to that cell's TAL subtype. Within each donor, Spearman
correlation was calculated between PAPPA2 expression and neighborhood exposure,
so each donor contributed one estimate. Statistical support was assessed with
1,000 within-donor permutations. In each permutation, gene-expression ranks
were shuffled across cells within a donor while neighborhood-exposure ranks were
held fixed; the same shuffled cell order was applied across the five focus
genes. The primary analysis required at least 5 cells per donor, and sensitivity
analyses required at least 3 and 10 cells. Benjamini-Hochberg (BH) correction
was applied across 60 gene, subtype, and comparison tests. Because neighborhoods overlap, conventional
correlation P values were treated as descriptive rather than as independent-sample
inference. The full permutation null distributions are archived in the
Zenodo analysis archive; per-test summaries are provided in Additional file 1.

### Analysis-unit comparison and internal robustness

Analysis-unit comparison tested whether the same gene, subtype, and disease
comparison led to the same direction across the four units. For each gene,
subtype, and comparison, we recorded the direction in donor pseudobulk, Milo
neighborhood abundance, neighborhood-expression association, and within-donor
cell-rank association. Concordance was defined as identical signs of the
standardized unit-specific effects, and discordance as opposite signs. Nominal
statistical significance was not required for direction classification. A
1,000-replicate donor bootstrap was used for complete-TAL pseudobulk effects.
The complete-TAL pseudobulk analysis, the previous subset pseudobulk analysis,
the neighborhood-level resampling result, and the within-donor permutation
result were joined into a single analysis-unit table. This table was used to
evaluate stability across units rather than to treat the most significant unit
as definitive. The target-panel hash and frozen gene lists are provided in
Additional file 1.

Internal robustness tested whether the donor-level directions were dependent
on individual donors, age or sex, or tissue region. We repeated the
pseudobulk analysis after removing one donor at a time, resampling donors
within condition strata, adjusting for age and sex, and stratifying by cortex
or medulla. These analyses were treated as stability checks and not as
additional discovery units.

Graph sensitivity used actual Milo runs rather than a fixed-kNN substitute.
The KPMP TAL neighborhoods were rebuilt with k = 15, 30, and 50 at neighborhood
proportion 0.1, and added a second k = 30 run with an independent Milo seed.
These graph and seed analyses were performed post hoc after the primary
k = 30 analysis.
At fixed k = 30 and seed 20261122, the Milo
neighborhood proportion across 0.05, 0.10, and 0.20. The prop = 0.10 run
served as the same-seed comparator.
To test representation dependence, a panel PCA30 space was built from the
frozen 2,400-gene expression cache. Genes were standardized, and Milo was
repeated at k = 30, prop = 0.10, and seed 20261122. This space represents a panel-PCA sensitivity
space, not a reproduction of a full-transcriptome PCA.
For each graph, the frozen 400-gene donor-within panel, the
effect-aware conditional null, and the full-TAL aTAL2 PAPPA2 cell-rank
statistic. Donor pseudobulk and composition-expression decomposition were not
re-estimated across graphs because they depend on donor and subtype
aggregation rather than on neighborhood membership.
Per-run cells, donors, neighborhood counts, significance counts, seeds, and
neighborhood proportions are reported in the graph-stability workbook.

### Unit-aware software framework

The analytical workflow was implemented as `unitaware` v0.3.0 to make estimand
definitions, reversal scoring, and calibration reproducible across analyses.
The platform-independent Python package provides functions for donor pseudobulk,
neighborhood-expression association, within-donor cell-rank association,
standardized reversal scoring, composition-expression decomposition, and
simulation. It also provides a conditional effect-pairing permutation for
calibrating reversal estimates against the observed effect-size structure and
a stability summary for graph, seed, or other calibration runs.
Analysis-unit definitions and thresholds are explicit inputs rather than
hidden defaults. The package includes a command-line interface (CLI), ten unit
tests, application programming interface (API) documentation, and an MIT
software license, and is distributed with the archived analysis code. The
versioned analysis record containing the package is archived under digital
object identifier (DOI) 10.5281/zenodo.24042041.

### Simulation benchmark

Simulation tested whether the framework distinguished genuine opposing
effects from null and concordant settings. Six prespecified scenarios were
simulated:
scenarios: null, donor-positive only, within-positive only, concordant positive,
donor positive with within negative, and donor negative with within positive.
Each scenario used 25 Reference and 25 disease donors, 30-80 cells per donor,
eight genes, donor random intercepts, 50 exposure-defined neighborhoods, and
100 replicates with seed 20261031. Donor random intercepts and residual noise
were drawn from normal distributions with standard deviations of 0.5 and 0.8,
respectively, and exposure was standard normal. Scenario-level donor effects
ranged from -2.0 to 2.0, with donor-effect standard deviations of 0.15-0.8.
Within-donor coupling terms ranged from -1.8 to 1.8, with standard deviations
of 0.10-0.5. The 100-replicate benchmark was designed to separate
scenario-level direction patterns, not to estimate extreme-tail error rates.
The same reversal threshold used in the empirical panel was applied.

### KPMP eGFR analysis

The KPMP eGFR analysis tested whether the donor-level PAPPA2-positive TAL
fraction was associated with baseline kidney function. PAPPA2-positive fraction
was defined as the proportion of TAL cells with at least one detected PAPPA2
transcript (raw count > 0). The public analytic release provided binned
baseline values rather than continuous eGFR. Baseline eGFR was therefore treated
as an ordinal outcome with four categories: 0-29, 30-59, 60-89, and 90 or more
mL/min/1.73 m². The primary ordinal logistic model included
PAPPA2-positive fraction per 10 percentage points, age, sex, and disease
category. Additional models added C-TAL-A, C-TAL-B, aTAL2, and frTAL fractions,
followed by mean detected genes and assay platform. Odds ratios below 1
indicate lower eGFR with higher PAPPA2-positive fraction.

### Candidate-gene specificity

Candidate-gene specificity tested whether the PAPPA2 direction was distinctive
relative to prespecified control genes. CXCL12, NOTCH3, HES1, and EPHB2 were
selected a priori as candidate comparators representing chemokine or stromal,
Notch, and ephrin-related programs. They were not expression-matched to PAPPA2;
the analysis therefore tested directional consistency rather than a
matched-null contrast. These genes were analyzed with the same donor-aware and
eGFR models used for PAPPA2, and the comparison assessed BH-corrected support
and tissue specificity.

### KPMP GSE183456 spatial cross-modality analysis

Spatial cross-modality analysis tested whether the PAPPA2-associated TAL state
co-occurred with a repair-failure spatial program within the same KPMP
consortium resource. The signatures were taken from prior KPMP spatial analyses
[5,11] and applied without re-selection or re-tuning. The frozen analysis used
nine main-branch genes:
THSD4, TEX41, EGF, CA10, CACNB2, CABP1, CXCL12, PRKG1, and NOS1. The
repair-failure signature used TMEM207, KLRD1, RNF144B, JAG1, TMEM52B, DACH1,
OSBPL10, LRMDA, and KCNJ10. For each section, genes were normalized by counts
per 10,000, transformed with log1p, and converted to within-section z scores.
Branch scores were the mean z scores of the signature genes, and high-state
thresholds were fixed at the 75th percentile to define the top expression
quartile within each section. The primary metric was the odds
ratio for co-occurrence of high main-branch and high repair-failure states.
The same parameters were applied to all 23 sections without re-tuning.

### GSE137570 public eGFR testing

GSE137570 was treated as a secondary external directional sensitivity analysis
of the PAPPA2-eGFR relationship because of its small sample size. Cohort 1
contained 24 donors with public eGFR, age, sex, and fibrosis. PAPPA2 expression
was extracted from the normalized read-count matrix. Spearman correlation was
calculated between PAPPA2 and eGFR, and a linear model was adjusted for age,
sex, and fibrosis. The same model was applied to the four control genes. The
analysis was treated as directional support, not definitive external
validation.

### GSE175759 locked clinical-association analysis

GSE175759 was used to test whether the frozen 11-gene TAL-injury score was
associated with continuous eGFR in an independent primary kidney-disease cohort
[16]. The public series contained 90 manually microdissected tubulointerstitial
samples. Diagnoses included 46 IgA nephropathy and 22 living-donor controls.
The remaining samples comprised nine minimal change disease, four membranous
nephropathy, three focal segmental glomerulosclerosis, three diabetic
nephropathy, and three lupus nephritis. Three submitter-marked technical
outliers were excluded before scoring, leaving 87 samples, including 43 with
IgA nephropathy. Age and sex were not available in the public annotations,
so no age- or sex-adjusted model was possible; available covariates were
diagnosis, library size, and TAL content.

The frozen 11-gene score was constructed without reweighting from the
independently published `TAL_Injury_h1` marker list (Supplemental Table S4)
associated with GSE284715/GSE284742 [17]. The score used the first ten markers
after excluding `NYAP2`, which was absent from GSE175759, plus `PAPPA2` as the
prespecified case gene. No GSE175759 outcome information was used for marker or
weight selection. The 11 genes were `PAPPA2`, `TMPRSS4`, `LIX1`, `WNT10A`,
`KRT7`, `CRIP1`, `S100A14`, `LAMC2`, `LAMB3`, `SOD3`, and `B3GALT5`; all 11
were detected in every original sample. Counts were converted to counts
per million (CPM) using the total count per sample, transformed as
log2(CPM + 1), and z-scored
within the 87-sample analysis cohort. The score was the equal-weight mean of
the gene-wise z scores.

The prespecified primary model was `eGFR ~ score` with
heteroskedasticity-consistent type 3 (HC3) standard errors. Prespecified secondary and
sensitivity models added diagnosis, library size, a PAPPA2-only score, and
leave-one-gene-out scores. A disease-only sensitivity model excluded
the 22 control samples. These models tested association, not causality or
clinical utility.

A locked secondary score-specificity protocol specified three analyses before
their execution. We added a fixed TAL-content score based on `UMOD`, `SLC12A1`,
`KCNJ1`, `CASR`, and `CLDN16` to the clinical model. We then drew 10,000
control 11-gene sets from 14,876 eligible genes after stratification into 20
bins by mean expression; each set was scored identically to the observed
signature. We also compared a PAPPA2-only model with a model containing PAPPA2
plus the remaining 10 score genes. These analyses tested specificity and
incremental information, not causality.

### Nephroseq GFR screening analysis

Nephroseq analysis tested whether the PAPPA2-GFR direction was observed across
multiple precomputed kidney datasets. Nephroseq v5 was queried through its
public API using the GFR filter. The query applied a fixed screen of
pValue < 0.05, abs(rValue) ≥ 0.5, and abs(foldChange) ≥ 1.5. The same screen
was applied to PAPPA2 and the four control genes. It prioritized large
precomputed correlations and was not treated as an inferential model or as an
estimate of association prevalence. The screen could therefore miss weaker
associations. Each returned analysis was summarized by dataset, tissue type,
GFR method, sample count, correlation coefficient, and P value. Analyses were
not treated as independent cohorts because several analyses can originate from
the same dataset. Direction and consistency were interpreted without pooling
analyses or treating their number as independent evidence.

### Independent non-KPMP cohort analyses

GSE211785 and GSE254185 were used for independent directional evaluation of
donor-level PAPPA2 disease elevation outside KPMP. GSE211785 was generated by
the Susztak laboratory and contains control kidney,
diabetic kidney disease (DKD), and hypertensive kidney disease (HKD) samples.
The non-KPMP merged object was used, and analyses were restricted to cells annotated
as cortical TAL (C_TAL) or medullary TAL (M_TAL) in the single-cell and
single-nucleus RNA libraries. After requiring at least
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
provided independent directional evaluation, not experimental or causal
validation.

### Cross-dataset benchmark of calibration methods

A harmonized benchmark was constructed across KPMP, GSE211785, and GSE254185
using the same frozen genes, reversal definition, and null-model hierarchy.
The benchmark compared three decision procedures: expression-matched
control-pool random sets, effect-aware pairing conditioned on the left unit,
and effect-aware pairing with the donor and within-donor units exchanged.
The prespecified primary endpoint was donor-versus-within reversal; the other
two unit pairs supplied secondary context.

For GSE211785, the public PreSCVI H5AD (Hierarchical Data Format 5 AnnData)
object was used, including its raw-count layer, public UMAP coordinates,
cell-state annotations, and sample labels. For GSE254185, submitter-provided
count matrices and cell metadata were used. Quality control, clustering, and
UMAP embeddings were not re-estimated. Counts were normalized to counts per 10,000 using
library-level totals and transformed with log1p. Donor pseudobulk required at
least five cells per donor, and within-donor analysis required at least 20
cells per donor. A state-comparison group required at least three disease and
three Reference donors.

Of the 400 target genes, 385 were available in all cohorts; of the 2,000
control-pool genes, 1,915 were available. One GSE254185 target gene lacked an
evaluable donor-within effect, leaving 384 evaluable target genes for that
endpoint. Genes absent from a cohort were excluded and were not reassigned to
another panel group.

For external cohorts, comparison-specific neighborhood exposure was calculated
within each library from a fixed 30-nearest-neighbor graph in public UMAP space,
excluding self-neighbors. For each cell, exposure was
`log2[((n_d + 0.5)/(N_d + 1))/((n_c + 0.5)/(N_c + 1))]`. Here `n_d` and
`n_c` were comparison-specific disease and Reference labels among the 30
neighbors, and `N_d` and `N_c` were the corresponding library totals.
Neighborhood-expression association was computed between this exposure and
each gene's neighborhood mean expression. This approach provides a kNN approximation of the
KPMP neighborhood unit, not an exact Milo reproduction. The external k = 15,
30, and 50 sensitivity analysis was retained, with k = 30 as the prespecified
primary scale.

The primary reversal threshold was an absolute standardized effect of at least
0.5 in both units. Expression-matched inference used 1,000 random sets of up
to 385 genes drawn from the 1,915 control genes. Effect-aware inference used
10,000 conditional permutations that preserved each target gene's left-unit
effect and the marginal distribution of the right-unit effect within each
state and comparison. The same permutation was repeated with the two units
exchanged. Threshold sensitivity used 0.5, 0.75, 1.0, and 1.25 with 5,000
permutations per threshold. The benchmark compared method decisions across
cohorts; it did not treat different TAL annotations as equivalent estimands.

### Matched cross-tissue portability analysis

To test whether the calibration behavior was kidney-specific, a
non-kidney portability analysis using the public GSE243917 synovial snRNA-seq
resource. Analyses were restricted to fibroblast-lineage nuclei. The contrast
was high lymphocytic inflammation (three donors) versus moderate inflammation
(eight donors); the single osteoarthritis donor was excluded. This contrast
represents a within-arthritis state comparison, not a kidney disease-control cohort.

The four cohorts were restricted to the same 166 target and 566 control genes
available in all datasets. One GSE254185 target gene lacked an evaluable
donor-within effect. The expression-matched random-set benchmark,
the donor-conditioned effect-aware null, and the order-swapped within-donor-
conditioned null using the same 0.5 threshold, 1,000 random sets, and 10,000
permutations. Thresholds of 0.5, 0.75, 1.0, and 1.25 were evaluated with 5,000
permutations. The analysis tested portability of the calibration procedure,
not biological equivalence between synovium and kidney.

As another non-kidney boundary test, the independent GSE138669 skin
case-control cohort was processed, comprising 12 systemic sclerosis and 10
healthy donors. Cell Ranger matrices were filtered to cells with at least 500 counts
and 200 detected genes, normalized to counts per 10,000, and log1p
transformed. A total of 3,000 highly variable genes were selected in a
sample-aware manner and adjusted for sample with ComBat. PCA, a
20-nearest-neighbor graph, UMAP, and Leiden clusters were then constructed at resolution 1.0. This produced
60,870 cells, 27 unsupervised states, and 380 target plus 1,884 control genes.
The cohort was used as an independent disease-control boundary analysis, not
as a positive replication requirement.

### Donor-stratified TAL signature transfer

To test whether the aTAL2-associated transcriptional program transferred
outside the discovery annotation, markers were derived for all four KPMP TAL
subtypes. Donors were divided into 156 training and 67 held-out donors. The top
35 markers per subtype were selected only in the training donors. Each cell
received the mean within-subtype marker z score, and the difference between the
top and second-ranked subtype score was defined as the marker margin. The
prespecified high-margin threshold was a marker margin of at least 0.2. We use
`high-margin aTAL2-like assignment` for cells above this threshold. The term
denotes a within-rule classifier result, not validated subtype identity.

Held-out KPMP cells were used to estimate coverage, accuracy, macro F1, and
aTAL2 precision and recall. The same frozen marker panel was then projected to
GSE211785 and GSE254185. In each external cohort, donor pseudobulk and
within-donor analyses were repeated in high-margin aTAL2-like assignments. This
analysis tests transfer of a transcriptional signature and its donor-within
behavior. It does not establish that an external cohort contains the
KPMP-defined aTAL2 subtype or that the labels are equivalent.

### Kidney regulatory and chromatin context

This analysis was performed as secondary biological-plausibility context. It did
not test calibration performance or establish a disease mechanism. Kidney
expression quantitative trait locus (eQTL) and chromatin context were evaluated to determine whether PAPPA2 had a
tubulointerstitial regulatory basis and whether the aTAL2 marker program showed
coordinated regulatory support. Gene-level eQTL results from NephQTL2, a
human kidney eQTL resource, were obtained from 240 glomerular and 311
tubulointerstitial human kidney samples [18]. PAPPA2 and the fixed
top-35 aTAL2 markers were evaluated without re-selecting genes. Gene-level
false discovery rate (FDR) values were reported for tubulointerstitial and
glomerular samples. PAPPA2 was then examined for significant
associations from a probabilistic transcriptome-wide association study (PTWAS)
for eGFR or urinary albumin-to-creatinine ratio. Absence from the significant
PTWAS list was recorded as non-significance, not as a negative association.

GSE211785 single-nucleus assay for transposase-accessible chromatin sequencing
(snATAC-seq) nuclei annotated as C_TAL or M_TAL were analyzed using the
gene-activity representation as a chromatin-accessibility proxy. Donor means
were calculated after requiring at least five snATAC-seq nuclei per donor and
state. Disease and control donors were compared with exact donor-label
permutation when feasible. Seven control and four disease donors contributed
to both states. Benjamini-Hochberg correction was applied across the fixed gene
and state tests. GSE211785 RNA and snATAC-seq are from the same cohort, so the
snATAC-seq result is an orthogonal-modality analysis rather than an independent
cohort.

### Reproducibility and software

Analysis parameters were fixed, and each reported result was linked to versioned code
and summary tables so that the reported analyses could be reproduced. Analysis
scripts, frozen analysis protocols, validation tables, and result hashes were
prepared for release in the public repository. The fixed-seed Milo analysis,
GSE183456 spatial analysis, and Nephroseq queries were rerun with fixed
parameters. Random seeds and versioned outputs are listed in the Availability
of data and materials section.

## Results

### Complete-data Milo revealed subtype divergence that masked opposing disease signals

Disease-related TAL neighborhoods differed by subtype in both abundance and
direction. The primary substrate was the complete set of 106,851 TAL cells from
the four annotated subtypes. A previously generated 4,446-cell derivative
matrix was retained as a comparator for the later coverage analysis.

In the combined disease comparison, C-TAL-A had
1,797 significant neighborhoods among 4,012 annotated neighborhoods, and
C-TAL-B had 715 among 2,805. aTAL2 had 632 among 1,465, and frTAL had 1,023
among 1,332. The directions were opposite: frTAL had 996 positive and 27
negative significant neighborhoods, whereas C-TAL-A had 333 positive and 1,464
negative neighborhoods. Their median neighborhood logFCs were +0.95 and -0.50,
respectively (Fig. 2b).

Disease-specific analyses showed that direction also depended on disease
context. In AKI, aTAL2 had 842 significant neighborhoods, 822 positive and 20
negative, with a median logFC of +1.36. In CKD, aTAL2 had 282 significant
neighborhoods, 34 positive and 248 negative, with a median logFC of -0.41.
C-TAL-A remained negative in AKI (-1.07) and CKD (-0.42), whereas CKD frTAL was
positive (+1.14; 956 positive and 26 negative significant neighborhoods).

These patterns show that a single neighborhood summary would mask opposite
directions across subtypes and disease contexts (Fig. 2a-c; Additional file 2).
A secondary subset-level Augur separability summary is provided in Additional
file 2 and was not used as a primary analysis unit. Fig. 3 summarizes the four
estimands and their interpretation boundaries. The next analysis tested whether coverage
differences alone explained the donor pseudobulk directions (Fig. 4).

[[FIGURE:Fig2]]

[[FIGURE:Fig3]]

### Complete versus subset coverage changed selected pseudobulk directions

The derivative subset preserved broad pseudobulk trends but changed selected
gene and subtype conclusions. In the complete analysis, 20 of 60 gene, subtype,
and comparison tests had FDR < 0.05 using unadjusted Mann-Whitney tests. After
adjustment for age and sex, 25 of 60 remained significant. Directions agreed
between the subset and complete analyses in 48 of 60 comparisons, whereas 12
comparisons changed direction. These 12 discordant tests were concentrated in
CXCL12, HES1, NOTCH3, and EPHB2, and none involved PAPPA2 (Fig. 4; Additional
file 2). Coverage alone changed one in five pseudobulk conclusions.
The next analysis compared all four analysis units directly in the complete TAL set.

[[FIGURE:Fig4]]

### Cross-unit concordance at the neighborhood level broke down within donors

PAPPA2 showed aligned donor-level and neighborhood-level signals that did not
extend to within-donor ranking. Across 12 comparison and subtype groups,
neighborhood PAPPA2 expression was positively correlated with Milo logFC, with
Spearman correlations from 0.029 to 0.729.
Donor-stratified intervals excluded zero in 11 of 12 groups. The exception was
AKI aTAL2, where the correlation was weak (rho = 0.029; 95% interval -0.086 to
0.255). Partial correlations controlling for neighborhood purity were similar
to the primary estimates (Fig. 5c; Additional file 2).

[[FIGURE:Fig5]]

This concordance did not extend to within-donor ranking: in combined disease,
aTAL2 donor pseudobulk was +0.411 and the neighborhood correlation was +0.145.
In contrast, the disease-minus-Reference difference in within-donor correlation
was -0.206 (FDR = 0.0086). In AKI, the corresponding within-donor difference
was -0.286 (FDR = 0.0086). frTAL showed the same qualitative reversal with a
weaker difference in combined disease (-0.090; FDR = 0.018). These estimates
answer different questions and cannot be treated as interchangeable
(Fig. 5a-d; Additional file 2).

The direct opposition raised a specific question. Was the excess reversal
frequency in biologically selected genes a property of those genes, or a
consequence of effect-size selection? This question was addressed with the
prespecified 400-gene panel.

### Effect-aware calibration exposed effect-size selection and graph-dependent residual

Expression matching alone overstated donor-versus-within reversal, whereas
effect-aware calibration removed that excess and exposed graph-dependent
residual. The observed fraction in the prespecified 400-gene panel was
19.0%, compared with 12.7% in expression-matched random sets. The empirical
P value was 0.000999, the minimum attainable value with 1,000 sets. Restricting to the
300 biological genes increased the observed fraction to 20.8% versus 12.7%. An
expression-stratified sensitivity analysis gave 19.0% versus 12.8% (empirical
P = 0.000999). Under expression matching alone, the panel appeared to
identify a biological excess.

The internal negative controls did not explain this pattern. Their
donor-versus-within reversal fraction was 13.6%, indistinguishable from random
100-gene sets from the external control pool (12.6%; two-sided empirical
P = 0.485). Reversal fractions differed across biological groups: 24.0% for
PAPPA2 coexpression genes, 21.6% for TAL subtype markers, and 16.8% for injury,
inflammatory, and metabolic genes.

Donor-versus-neighborhood reversal was also higher than in expression-matched
sets (21.9% versus 14.5%; empirical P = 0.000999). Neighborhood-versus-within
reversal showed a smaller full-panel excess (21.0% versus 19.6%; empirical
P = 0.005), but the biological 300-gene subset did not (P = 0.189).

The effect-aware null changed the donor-versus-within conclusion in the
prespecified panel. The observed fraction was 19.0%, while the donor-conditioned
null was 20.8% (one-sided empirical P = 0.999). The observed value was
below the expectation generated by preserving the effect-size marginals. The
apparent excess over expression-matched genes disappeared, and expression
matching materially underestimated the reversal expected from effect-size
structure (12.7% versus 20.8%). The order-swapped null was not equivalent:
conditioning on within-donor effects produced an expected reversal fraction of
14.9%, and the observed 19.0% exceeded it (one-sided empirical
P = 0.0001). This result makes the calibration conditional on donor-level effects.
Donor-versus-neighborhood reversal showed the same loss of evidence under the
donor-conditioned null (observed 21.9%; expected 21.2%; P = 0.111).

The next analysis tested whether this panel-level result depended on the Milo graph. The
expression-matched excess persisted in every actual-Milo run: target-panel
reversal ranged from 19.0% to 22.6%, compared with 12.7% to 13.7% in matched
sets. The donor-conditioned effect-aware null ranged from 20.8% to 21.1%. The
residual after calibration was not stable (Fig. 8b). k = 15 gave +1.8
percentage points (empirical P = 0.0008), and k = 50 gave +1.2 points
(P = 0.019). The original k = 30 run gave -1.8 points (P = 0.999), while a
second k = 30 seed gave +0.7 points (P = 0.125). Expression matching consistently
underestimated effect-size-driven reversal, whereas the conclusion of a
universal calibrated null was not graph- or seed-invariant. The panel supports
effect-aware calibration as a required sensitivity standard, not as proof that
all kidney panel excess disappears.

The actual-Milo graph and seed perturbations were restricted to KPMP. External
kidney cohorts used the fixed 30-neighbor exposure approximation described in
the Methods; they test portability of the calibration problem, not replication
of the complete graph-stability envelope. This boundary prevents the external
benchmark from being interpreted as independent confirmation of KPMP graph
robustness.

Neighborhood-versus-within reversal retained a small conditional excess
(21.0% versus 19.4%; P = 0.0017). This result was threshold-sensitive and was
not significant above the 0.5 threshold. Individual gene metrics, random-set
results, conditional null distributions, and threshold sensitivity are
provided in Additional file 1, with calibration plots in Supplementary Figs. S1-S2.

The panel analysis exposed two distinct inferential failures in cross-unit interpretation. The
conventional expression-matched null overinterpreted an apparent biological
signal, while the residual after effect-aware calibration was not stable enough
to support a universal negative conclusion. The resulting standard is not one
fixed P value or one fixed graph. It requires reporting the observed reversal,
the effect-aware expectation, the conditioning direction, and the range across
prespecified graph and seed variants. This panel-level conclusion does not
erase the descriptive aTAL2 direction contrast. Simulation next tested
whether known opposing effects could be recovered under the framework.

### Simulation benchmark supported detection of opposing versus concordant unit effects

The simulation benchmark recovered strong opposing effects while retaining an
explicit null expectation. Donor-versus-within reversal was detected in 95.5%
of simulations when the true donor effect was positive and the within-donor
effect negative, and in 97.3% of reciprocal simulations. Concordant positive
effects produced a 2.3% reversal fraction,
donor-only and within-only scenarios produced intermediate fractions, and the
null scenario produced 20.1%. The null fraction reflects the prespecified
effect-size threshold and finite sampling, not a true opposing effect. This
simulation separated scenario-level patterns and was not designed to estimate
extreme-tail error rates. The framework detected strong opposing
patterns while retaining explicit null calibration (Additional file 3). The next
examined subtype-level conflicts that were not resolved by the
panel-level conclusion.

### aTAL2 provided a UMAP-defined gene-level conflict beneath the calibrated panel

aTAL2 contained a stable gene-level direction conflict beneath the
panel-level result. PAPPA2 donor-level and within-donor directions diverged in
the same subtype.

Complete donor pseudobulk showed higher PAPPA2 expression in aTAL2 for all
three comparisons. Mean differences were +0.411 for combined disease, +0.265
for AKI, and +0.484 for CKD, with FDR < 0.05 for each. The same positive
direction was retained in the age- and sex-adjusted models.

Neighborhood-expression association was positive but less stable: the full-data
Spearman correlation was 0.145 for combined disease and 0.029 for AKI, and the
AKI donor-stratified bootstrap interval included zero. This result
indicates that neighborhood expression and donor pseudobulk were aligned in
direction but not equally strong.

The within-donor estimate moved in the opposite direction: the
disease-minus-Reference rank difference was -0.206 in combined disease and
-0.286 in AKI, with FDR = 0.0086 for each. In CKD, the difference was also
negative (-0.089) but did not pass FDR (0.224).

The aTAL2 conflict remained negative across UMAP k, seed, and proportion
sensitivity. Disease-minus-Reference cell-rank differences were -0.210 at
k = 15 and -0.206 at the original k = 30. They were -0.161 at a second k = 30
seed and -0.163 at k = 50, with label-permutation P = 0.000999 for each. The fixed-k
proportion runs for this cell-rank difference were also negative and
significant. PAPPA2 donor pseudobulk
remained positive across all settings, and the composition-expression
decomposition is graph-invariant by construction because it uses donor subtype
proportions and subtype means (Table 1).

**Table 1. Calibration-stability audit across UMAP actual-Milo settings.**

| Endpoint | Effect range across four settings | Significant settings | Stability class |
|---|---:|---:|---|
| Expression-matched target excess | +0.063 to +0.092 | 4/4 | Stable positive support |
| Donor-conditioned panel residual | -0.018 to +0.018 | 2/4 | Direction-unstable |
| aTAL2/PAPPA2 cell-rank conflict | -0.210 to -0.161 | 4/4 | Stable negative support |

**Table 1 legend.** The four UMAP settings are post hoc actual-Milo runs at
k = 15, the original k = 30, a second k = 30 seed, and k = 50. For the first two
endpoints, effects are differences in donor-within reversal fractions; +0.0167
therefore equals +1.67 percentage points. Positive values indicate higher
observed reversal than the comparator. The expression-matched endpoint is
target reversal minus the expression-matched expectation. The panel residual
is target reversal minus the donor-conditioned effect-aware null. The case
endpoint is the disease-minus-Reference within-donor cell-rank difference.
Significant settings denote empirical P < 0.05. The stability class is a
descriptive summary of direction and decision consistency, not a combined P
value. Underlying tables are retained in the archived analysis outputs.

The fixed-k proportion runs gave the same qualitative boundary. At k = 30 and
seed 20261122, prop = 0.05, 0.10, and 0.20 produced donor-conditioned
residuals of +0.0167, +0.0068, and +0.0080. Empirical P values were 0.0030,
0.125, and 0.089, respectively. All residuals were positive, but the significance
decision was not stable. In the same runs, the aTAL2/PAPPA2 cell-rank
differences remained negative (-0.144, -0.161, and -0.160) and were significant
for all three proportions. Proportion changed the calibrated
decision boundary without reversing the stable case-level conflict.

Embedding space introduced a stronger boundary (Fig. 8c,d). At the same
settings, panel PCA30 changed the donor-conditioned residual. The estimate moved
from +0.0068 (P = 0.125) in UMAP to -0.0077 (P = 0.915), a direction reversal. The
aTAL2/PAPPA2 cell-rank difference remained negative but attenuated from
-0.161 (P = 0.000999) in UMAP to -0.032 (P = 0.435) in PCA30. The aTAL2
conflict is UMAP-graph-stable and proportion-stable, but not
embedding-invariant (Table 2). This does not affect PAPPA2 donor pseudobulk or the
composition-expression decomposition, which are not defined by the
neighborhood graph.

**Table 2. Embedding sensitivity at fixed k, prop, and seed.**

| Endpoint | UMAP2D estimate | Panel PCA30 estimate | Stability class |
|---|---:|---:|---|
| Expression-matched target excess | +0.080 | +0.073 | Stable positive support |
| Donor-conditioned panel residual | +0.0068 (P = 0.125) | -0.0077 (P = 0.915) | Direction-unstable |
| aTAL2/PAPPA2 cell-rank conflict | -0.161 (P = 0.000999) | -0.032 (P = 0.435) | Direction-stable, decision-unstable |

**Table 2 legend.** Both runs used k = 30, prop = 0.10, and seed 20261122.
For the first two endpoints, effects are differences in donor-within reversal
fractions; the third endpoint is the within-donor cell-rank difference. UMAP2D
uses the atlas UMAP coordinates. Panel PCA30 uses a 30-component
principal-component representation of standardized expression from the frozen
2,400 genes. It is a representation sensitivity analysis, not a
full-transcriptome PCA.
Direction-stable, decision-unstable means that the sign was preserved but
statistical significance changed.

Fig. 6c provides context for this choice. Across five genes and three
comparisons, 15 tests were possible per subtype. Cross-unit conflicts occurred
in 12 of 15 aTAL2 tests, 10 of 15 C-TAL-A tests, 10 of 15 C-TAL-B tests, and
11 of 15 frTAL tests. Conflict count alone did not distinguish aTAL2.
Its PAPPA2 case was distinguished by the larger within-donor effect and
significant support in two of three disease comparisons (Fig. 6; Additional
file 2). The next analysis asked whether frTAL showed the same qualitative pattern.

[[FIGURE:Fig6]]

### frTAL showed a parallel but statistically weaker conflict

frTAL showed the same qualitative direction conflict as aTAL2, with weaker
statistical support. Donor pseudobulk and neighborhood-expression association
were positive. In contrast, within-donor disease-minus-Reference rank
differences were negative in all three comparisons. Estimates were -0.090 in
combined disease (FDR = 0.018), -0.117 in AKI (FDR = 0.146), and -0.056 in CKD
(FDR = 0.224). Thus, frTAL supports the generality of analysis-unit sensitivity
without matching the magnitude or FDR support of the aTAL2 conflict. PAPPA2
then provided a benchmark case across all four TAL subtypes.

### PAPPA2 donor elevation provided a graph-invariant benchmark with subtype-specific expression

PAPPA2 provided a benchmark case that remained directionally positive within
every subtype and comparison. Its donor pseudobulk mean difference was positive
in all 12 groups, and 11 passed FDR < 0.05. AKI frTAL remained positive but did
not pass FDR (0.080). In combined disease, the mean differences were 0.894 for
C-TAL-A, 0.477 for C-TAL-B, 0.411 for aTAL2, and 0.370 for frTAL. The same
positive direction was retained in the age- and sex-adjusted models (Additional
file 2). Leave-one-donor neighborhood-expression correlations for PAPPA2 are
provided in Additional file 2.

The direction remained positive in leave-one-donor refits and donor bootstrap
resampling (Fig. 7a-b). These panels summarize the same donor-level estimate
under different stability checks and are not independent validations. Cortex
and medulla mean differences agreed in direction in 11 of 12 groups. The
exception was AKI frTAL in medulla, a negative stratum estimate based on two
disease donors; the corresponding complete pseudobulk result was not
significant.

Composition-expression decomposition could not explain the aggregate difference
through subtype mixing. In combined disease, the observed difference was +0.563
(95% bootstrap interval +0.483 to +0.639). The subtype-expression component was
+0.641 (95% bootstrap interval +0.560 to +0.721). The composition component was
+0.010 (-0.003 to +0.023), and the interaction/covariance remainder was -0.088
(-0.133 to -0.044). Expression accounted for 86.7% of the absolute component
magnitude, compared with 1.3% for composition.

The same pattern was present in disease-specific comparisons. In AKI, observed,
expression, and composition components were +0.379, +0.522, and -0.009,
respectively. In CKD, the corresponding values were +0.650, +0.698, and +0.019
(Additional file 2). Higher PAPPA2 in disease donors was not
explained by a shift toward subtypes with higher baseline PAPPA2. The
decomposition identifies the source of the aggregate difference, not a
cell-intrinsic mechanism. The next analysis tested whether the PAPPA2-positive fraction
was associated with kidney function.

[[FIGURE:Fig7]]

### Higher PAPPA2-positive fraction tracked lower eGFR in KPMP

In KPMP donors, higher TAL PAPPA2-positive fraction was associated with lower
eGFR in ordinal models. The primary model included 202 KPMP donors and adjusted
for age, sex, and disease category. Each 10 percentage-point increase in PAPPA2
fraction was associated with an odds ratio of 0.819
(95% confidence interval [CI] 0.682-0.985, P = 0.0335). Adding TAL subtype fractions moved the OR farther below 1
(0.696, 95% CI 0.530-0.915), and adding technical covariates gave OR 0.661
(95% CI 0.500-0.875). The confidence intervals for the primary and
subtype-adjusted models overlapped, so the numerical strengthening alone was
not interpreted as a formal suppressor effect. The association direction
remained unchanged across these models.

Disease-specific models showed that the signal was concentrated mainly in CKD.
The primary model in disease donors
included 130 donors and produced an OR of 0.789 (95% CI 0.648-0.960). CKD
included 94 donors and produced an OR of 0.733 (95% CI 0.582-0.924). AKI
included 36 donors and did not show a clear association. These results support
an internal adjusted association between PAPPA2 and kidney function (Additional
file 4). The next comparison assessed PAPPA2 against prespecified control genes.

### Control-gene comparisons contextualize PAPPA2 specificity

PAPPA2's within-donor support did not exceed that of CXCL12, so specificity
rested on donor-level and tissue-context evidence. In the current within-donor
rank analysis, PAPPA2 had five FDR-supported comparisons, compared with six for
CXCL12 and none for NOTCH3, HES1, and EPHB2. CXCL12 also showed a stronger
negative-direction tendency. The PAPPA2 case rests on its donor-level consistency across all 12
subtype comparisons, its external donor-level direction, and its
tubulointerstitial-negative GFR pattern rather than on a superior within-donor
comparison count.

In Nephroseq, PAPPA2 had 10 significant GFR analyses, of which nine were
negative. In the same screen, all seven tubulointerstitial analyses for PAPPA2
were negative. HES1 also had
six negative tubulointerstitial analyses but one positive tubulointerstitial
analysis. EPHB2 had three negative tubulointerstitial analyses and additional
mixed-direction associations in glomeruli, cortex, and blood. NOTCH3 had two
negative tubulointerstitial analyses, and CXCL12 had one. The external screen
supports a tubulointerstitial-negative context for PAPPA2, while the
within-donor analyses do not establish PAPPA2-specific control-gene superiority.
The next analysis asked whether the PAPPA2-associated programs were spatially coordinated.

### Spatial cross-modality analysis placed the PAPPA2-associated program in tissue context

The PAPPA2-associated branch programs were spatially coordinated but not
disease-specific. We projected the frozen signatures onto the KPMP Visium
resource GSE183456. It belongs to the same KPMP consortium as the discovery
snRNA-seq data and provided cross-modality context rather than an independent
cohort. Across 23 sections, 20 had a co-occurrence odds ratio
above 1. All 11 DKD sections and all 10 deduplicated DKD donors remained above
1. This pattern indicates spatial coordination of the two PAPPA2-associated
programs within kidney tissue.

The co-occurrence was not disease-specific: five of six Reference sections also
had an odds ratio above 1, and the DKD versus Reference comparison had P = 0.078.
The result supports a spatially coordinated tissue-state relationship
within KPMP. It does not establish disease specificity or external spatial
validation. The next analysis evaluated the PAPPA2-GFR direction in external kidney
datasets.

### External kidney datasets provided directional PAPPA2-GFR support

External kidney datasets supported a tubulointerstitial-negative PAPPA2-GFR
direction. In the small secondary sensitivity dataset GSE137570, 24 donors had
public eGFR. PAPPA2 showed a negative unadjusted correlation with eGFR
(Spearman rho = -0.365), but the adjusted model did not reach significance
(P = 0.351; Additional file 4).

Nephroseq provided the broader directional assessment across 50 precomputed GFR
analyses. PAPPA2 had 10 significant results, of which nine were negative.
All seven tubulointerstitial analyses were negative, with correlations from
-0.506 to -0.837. The sole positive result
was from a cadaveric transplant kidney analysis (Additional file 4).

The clustering of negative results in tubulointerstitial datasets is consistent
with the tubule-centered PAPPA2 signal. The analyses are precomputed and may
overlap, so they provide multi-dataset association support rather than
independent cohort replication. The frozen TAL-injury score was then evaluated
in an independent primary kidney-disease cohort.

### A frozen TAL-injury score was specifically associated with lower eGFR in an independent cohort

The frozen 11-gene TAL-injury score was associated with lower eGFR in an
independent cohort. In GSE175759, all 11 score genes were detected in all 90
original samples, and 87 samples remained after exclusion of three technical
outliers. Higher score was associated with lower eGFR: -7.98 points per
1-standard deviation (SD) increase
(95% CI -12.29 to -3.66; P = 0.00050). Diagnosis adjustment preserved the
association (-8.71; 95% CI -13.79 to -3.63; P = 0.0012), and the disease-only
estimate remained negative (-9.27; P = 0.0027; Table 3).

The association was not explained by bulk TAL content. After adjustment for the
fixed TAL-content score and diagnosis, the estimate was -10.13 points per SD
(95% CI -15.70 to -4.56; P = 0.00062). Across 10,000 expression-matched random
11-gene sets, the observed score coefficient was at the 0.01 percentile of the
null distribution (one-sided empirical P = 0.0002; null mean -3.25). This places
the clinical association beyond both bulk TAL content and expression-matched
gene-set composition.

PAPPA2 alone was also associated with lower eGFR (-9.06; 95% CI -14.08 to
-4.03; P = 0.00067). Adding the remaining 10 score genes produced a residual
coefficient of -4.10 (95% CI -8.11 to -0.10; P = 0.0478). The nested F test was
not significant (P = 0.089), and R2 increased from 0.169 to 0.197. These data
do not provide clear evidence of incremental clinical information from
the remaining 10 genes. The 11-gene score is interpreted as a state-level summary, not as
evidence that each component contributes independently. The IgA-only analysis
did not reach significance, while library-size adjustment and all
leave-one-gene-out models retained a negative direction (Additional file 4).

**Table 3. Independent clinical association of the frozen TAL-injury score
with eGFR in GSE175759.**

| Analysis | n | eGFR change per 1-SD score | 95% CI | P |
|---|---:|---:|---:|---:|
| Primary, unadjusted | 87 | -7.98 | -12.29 to -3.66 | 0.00050 |
| Diagnosis-adjusted | 87 | -8.71 | -13.79 to -3.63 | 0.0012 |
| Disease samples only | 65 | -9.27 | -15.10 to -3.45 | 0.0027 |
| TAL-content and diagnosis adjusted | 87 | -10.13 | -15.70 to -4.56 | 0.00062 |

**Table 3 legend.** Estimates are changes in eGFR (mL/min/1.73 m²) per
1-standard-deviation increase in the frozen TAL-injury score and are from
ordinary least-squares models with heteroskedasticity-consistent type 3 (HC3)
standard errors. n denotes analyzed samples, CI confidence interval, and SD
standard deviation. The cohort provides independent external
clinical-association support in primary kidney disease. Age and sex were
unavailable. Additional disease-specific, IgA-only, expression-matched null,
and leave-one-gene-out models are provided in Additional file 4.

### Non-KPMP single-cell cohorts supported donor-level PAPPA2 elevation without resolving within-donor coupling

Donor-level PAPPA2 elevation was reproduced in both non-KPMP cohorts, whereas
within-donor coupling remained unresolved. In GSE211785, donor-level PAPPA2
expression was higher in disease than control donors in both TAL compartments.
C_TAL showed a difference of +0.283 (nominal
P = 0.0315; FDR = 0.0630), and M_TAL showed +0.200 (P = 0.105; FDR = 0.133).
DKD-specific contrasts were larger than HKD-specific contrasts, and
leave-one-donor estimates remained positive. Within-donor injury-score
estimates were negative in C_TAL and M_TAL, whereas neighborhood-logFC
estimates were negative in C_TAL and near zero positive in M_TAL. None reached
statistical support.

GSE254185 provided another test in cTAL1, cTAL2, mTAL, TAL Injured, and
TAL Inflammatory states. Donor pseudobulk PAPPA2 was higher in mTAL (+0.341;
nominal P = 0.0121; FDR = 0.0604) and TAL Injured (+0.168). The
neighborhood-expression association in TAL Injured was negative
(rho = -0.096; FDR = 7.1e-05), and within-donor differences were negative using
injury-marker and neighborhood-logFC exposure. The donor-level direction was
consistent across annotation systems.
Within-donor estimates were negative or near zero, and none were statistically
resolved (Table 4).

**Table 4. PAPPA2 directional support in external single-cell cohorts.**

| Cohort and context | TAL annotations | Directional result |
|---|---|---|
| GSE211785: control, diabetic kidney disease, and hypertensive kidney disease | C_TAL and M_TAL | Donor level positive in both compartments, with nominal support in C_TAL; within-donor injury-score estimates negative and the M_TAL neighborhood-logFC estimate near zero positive; none statistically supported |
| GSE254185: normal kidney and obstructive nephropathy | cTAL1, cTAL2, mTAL, TAL Injured, and TAL Inflammatory | Donor level positive in mTAL and TAL Injured, with nominal support in mTAL; within-donor estimates negative but not statistically supported |

**Table 4 legend.** Directional results are summarized from the main text.
Nominal support denotes unadjusted P < 0.05; no donor-level result in this table
survived FDR < 0.05. Effect sizes, permutation P values, false discovery rates,
and donor counts are provided in Additional file 5. TAL annotations are
cohort-defined and differ between datasets, so the table evaluates directional
consistency rather than label equivalence. The within-donor results do not
distinguish absence of an effect from limited precision.

### Cross-dataset benchmark showed direction-dependent calibration

Calibration decisions depended on conditioning direction across all three
cohorts. The same frozen 385-gene target panel and 1,915-gene control pool were
applied to KPMP, GSE211785, and GSE254185. KPMP contributed 12
state-comparison groups, GSE211785 contributed six, and GSE254185 contributed
two. At the prespecified standardized-effect threshold of 0.5, expression
matching identified a target-panel excess in all three cohorts (3/3, empirical
P < 0.05; Table 5 and Fig. 8a).

Effect-aware pairing changed the decision in a direction-dependent manner.
When donor-level effects were held fixed, the observed donor-within reversal
fraction did not exceed the conditional expectation in any cohort (0/3,
P >= 0.49). KPMP showed 19.0% observed versus 20.6% expected
(P = 0.9975). GSE211785 showed 19.0% versus 19.8% (P = 0.8816), and
GSE254185 showed 18.8% versus 18.7% (P = 0.4937). Exchanging the conditioning
direction retained an excess in KPMP (19.0% versus 15.0%, P < 0.0001) and
GSE254185 (18.8% versus 15.7%, P = 0.0053). It did not retain an excess in
GSE211785 (19.0% versus 18.3%, P = 0.155).

The direction contrast was stable across thresholds of 0.5, 0.75, 1.0, and
1.25. Across the 12 cohort-threshold combinations, expression matching was
significant in 12/12, donor-conditioned calibration in 0/12, and
within-donor-conditioned calibration in 8/12. The conditioning direction is
therefore part of the estimand rather than a nuisance sensitivity analysis.

Secondary unit pairs were not uniformly positive. In KPMP, donor-versus-
neighborhood expression showed an expression-matched excess that was absent
under donor conditioning, whereas neighborhood-versus-within retained
conditional support. In GSE211785, neither secondary pair showed a conditional
excess. In GSE254185, neighborhood-versus-within but not
donor-versus-neighborhood retained conditional support. Full cohort, gene,
random-set, permutation, and threshold tables are provided in Additional
file 5.

**Table 5. Three-cohort calibration benchmark for donor-within reversal.**

| Cohort | Cells and groups | Observed | Expression P | Donor P | Within-donor P |
|---|---:|---:|---:|---:|---:|
| KPMP | 106,851 cells; 12 groups | 19.0% | 0.0010 | 0.9975 | <0.0001 |
| GSE211785 | 29,847 cells; 6 groups | 19.0% | 0.0010 | 0.8816 | 0.155 |
| GSE254185 | 9,912 cells; 2 groups | 18.8% | 0.016 | 0.4937 | 0.0053 |

**Table 5 legend.** Observed values are the proportions of evaluable
gene-comparison pairs showing donor-within direction reversal. Groups are the
cohort-specific TAL subtype and disease-comparison units used for calibration.
P values are one-sided empirical values. Expression P compares the target panel
with 1,000 random sets from the 1,915-gene control pool. Donor P and
within-donor P compare the observed fraction with donor-conditioned and
within-donor-conditioned (order-swapped) effect-aware nulls, respectively. Each
null uses 10,000 permutations. The donor-conditioned null preserves donor
effects and the marginal distribution of within-donor effects; the
within-donor-conditioned null preserves within-donor effects and the marginal
distribution of donor effects.
Up to 385 target and 1,915 control genes were evaluable; GSE254185 had 384
evaluable target genes for this endpoint. The benchmark tests portability of the
calibration procedure and its conditioning-direction dependence, not universal
biological enrichment. Actual-Milo graph and seed sensitivity for KPMP is
reported in Additional file 6; the donor-conditioned result was not
graph-invariant.

[[FIGURE:Fig8]]

### Matched cross-tissue benchmark retained a context-specific signal

Effect-aware calibration retained a context-specific synovial signal while
removing the kidney panel excess. All four cohorts were restricted to the 166
target and 566 control genes available in every dataset. Expression matching
again suggested a target excess in all four cohorts (4/4, P < 0.05).
Donor-conditioned calibration was not significant in the three kidney cohorts
(0/3), but retained a significant excess in synovial fibroblasts (20.5%
observed versus 14.3% expected; P = 0.0067). The
order-swapped within-donor-conditioned test was significant in KPMP and
synovium but not in GSE211785 or GSE254185 at the 0.05 threshold
(Table 6; Supplementary Fig. S3).

The synovial signal was stable across thresholds of 0.5, 0.75, 1.0, and 1.25:
donor-conditioned and within-donor-conditioned P values remained below 0.05 in
all four analyses. The kidney donor-conditioned results remained
non-significant across the same thresholds. Thus, effect-aware calibration did
not simply remove all cross-unit structure. It removed an apparent kidney
panel excess while retaining a non-kidney tissue-state signal supported under
both conditioning directions.

**Table 6. Matched cross-tissue calibration benchmark for donor-within reversal.**

| Cohort | Cells and groups | Target genes | Observed | Expression P | Donor P | Within-donor P |
|---|---:|---:|---:|---:|---:|---:|
| KPMP | 106,851 cells; 12 groups | 166 | 22.6% | 0.0010 | 0.6995 | <0.0001 |
| GSE211785 | 29,847 cells; 6 groups | 166 | 20.0% | 0.0080 | 0.7257 | 0.0839 |
| GSE254185 | 9,912 cells; 2 groups | 165 | 22.1% | 0.048 | 0.4264 | 0.0517 |
| GSE243917 synovium | 8,254 fibroblasts; 1 contrast | 166 | 20.5% | 0.014 | 0.0067 | 0.0103 |

**Table 6 legend.** All cohorts were restricted to a common set of 166 target
and 566 control genes; 165 target genes were evaluable in GSE254185. Observed
values are the proportions of evaluable gene-comparison pairs showing
donor-within direction reversal. The synovial contrast was high versus moderate
lymphocytic inflammation among arthritis donors, not a kidney case-control
comparison. Expression P compares the target panel with control-pool random
sets. Donor P and within-donor P are one-sided empirical values from
donor-conditioned and within-donor-conditioned effect-aware nulls,
respectively, each using 10,000 permutations. Donor conditioning holds donor
effects fixed; within-donor conditioning holds within-donor effects fixed. The
table tests portability of the calibration procedure, not tissue equivalence.

An independent skin case-control cohort provided a complementary boundary.
GSE138669 contained 12 systemic sclerosis and 10 healthy donors, yielding
60,870 cells and 27 unsupervised states. In this cohort, donor-within reversal
was 7.3% observed versus 8.0% expression-matched (P = 0.949). It was 7.3%
donor-conditioned (P = 0.641) and 8.9% under swapped conditioning (P = 1.0).
Donor-versus-neighborhood also showed no excess (observed 4.6%;
expression-matched P = 0.998; both conditional P = 1.0).
Neighborhood-versus-within had an expression-matched excess (17.9% versus
16.6%, P = 0.002). This was not retained under primary conditioning
(17.3%; P = 0.069) or swapped conditioning (18.9%; P = 0.980). Thus, the
framework did not manufacture a calibrated positive signal in an independent
non-kidney case-control cohort. Together with the synovial severity contrast,
this shows that non-kidney portability is context dependent rather than
uniformly positive.

The next analysis asked whether the frozen aTAL2 marker program transferred outside the
discovery annotations and retained donor-within behavior.

### Donor-stratified signature transfer supported a transferable aTAL2-like program

The frozen aTAL2-like program transferred at the donor level but did not
resolve within-donor coupling. The four-subtype marker panel assigned 86.66% of
held-out KPMP cells at the prespecified margin threshold. Among assigned cells,
accuracy was 93.90%,
macro F1 was 93.17%, aTAL2 precision was 83.28%, and aTAL2 recall was 93.93%
(Additional file 5). These held-out metrics indicated that the marker
margin was not an arbitrary threshold and supported cautious projection to
non-KPMP cohorts.

In GSE254185, the frozen rule yielded 1,176 high-margin aTAL2-like assignments,
predominantly from TAL Inflammatory and TAL Injured states. Their donor-level
PAPPA2 difference was +0.554 (permutation P = 0.0043), whereas the
within-donor injury-score direction was negative but not statistically
supported.

In GSE211785, 3,751 high-margin aTAL2-like assignments were mostly derived
from M_TAL, with a smaller C_TAL contribution. The donor-level PAPPA2
difference was +0.570 (P = 0.0235). The within-donor neighborhood-logFC
direction was negative but not statistically supported (Additional file 5).

The external assignments supported transfer of the donor-level
direction and a consistent negative within-donor trend. The within-donor
effects were not resolved, and the assignments did not establish external
aTAL2 identity or exact label equivalence. The high-margin assignments in
GSE254185 were enriched for inflammatory and injured states, not a clean
cortical aTAL2 population. Kidney eQTL and chromatin context were then evaluated
for PAPPA2 and the aTAL2 marker program.

### Regulatory context supported a tubulointerstitial origin without detected disease chromatin change

PAPPA2 showed a tubulointerstitial regulatory context without detected disease
chromatin change. In a secondary biological-plausibility analysis, PAPPA2 was a
strong tubulointerstitial eGene in NephQTL2 (gene-level FDR = 1.64e-14; model
size 1.035). It had no glomerular eQTL support (FDR = 0.526) [18] (Additional
file 5). This compartment-restricted regulatory association is consistent with
the TAL anatomical origin. It does not establish a disease effect.

PAPPA2 was not among the FDR-significant NephQTL2 PTWAS associations for eGFR or
urinary albumin-to-creatinine ratio. This indicates an absence of PTWAS support,
not evidence against the expression association. The fixed aTAL2 marker set was
also not enriched for tubulointerstitial eGenes: 11 of 35 markers were eGenes,
with an exploratory odds ratio of 0.58.

In GSE211785 snATAC data, PAPPA2 gene activity did not show higher disease
accessibility. The disease-minus-control difference was -0.029 in C_TAL and
-0.016 in M_TAL (FDR = 1.0 for both; Additional file 5). This orthogonal
analysis from the same cohort did not detect a disease-associated chromatin
increase. It does not exclude regulatory effects in other cell states, assays,
or disease contexts. Together, these data support a tubulointerstitial
regulatory context without a detected chromatin-level disease change. The next
tested whether candidate pathways had unified cross-unit support.

### Candidate pathways lacked unified cross-unit support

Candidate JAG1-NOTCH3-HES1 and EPHB2-STAT6 axes lacked unified cross-unit
support. Consistent direction was required across subtypes, comparisons, units,
and external GFR analyses, but this condition was not met. NOTCH3 was significant
only in CKD frTAL donor pseudobulk, HES1 showed mixed directions across
subtypes and GFR analyses, STAT6 showed no significant donor-level change, and
EPHB2 had subtype-dependent directions.

These boundary tests did not promote sporadic single-unit signals into a
pathway model. The findings do not exclude
other mechanisms, but they leave the upstream regulators of PAPPA2 undefined.
These boundary tests complete the Results and lead to the framework-level
Discussion.

## Discussion

### From cell states to estimands: analysis units define distinct targets

The new insight is that donor, neighborhood, and within-donor associations are
distinct estimands whose apparent agreement or conflict cannot be interpreted
without effect-aware calibration and a declared graph/embedding specification.
In the KPMP TAL resource, the same cell population produced positive donor and
neighborhood-expression signals and a negative within-donor signal in aTAL2
(FDR = 0.0086). Treating analysis units as interchangeable changes the
biological scale assigned to an association. This conclusion concerns
interpretation; it does not establish which unit is biologically primary.

### Effect-aware calibration is directionally conditional and graph-sensitive

Effect-aware calibration showed that cross-unit reversal is conditional on
effect-size structure, graph construction, and representation. The prespecified
400-gene panel addressed this limitation through a multi-gene design.
Expression-matched random sets preserved gene expression and detection but not
the joint distribution of effects across units. The target panel initially
showed 19.0% donor-versus-within reversal,
compared with 12.7% in expression-matched sets. Under the effect-aware null,
the expected reversal fraction was 20.8% and the observed value was below that
expectation (one-sided empirical P = 0.999). The calibration turns an apparent
coupling signal into a falsification test. The conventional null overinterpreted
effect-size selection, while the donor-conditioned null showed no excess
reversal in the prespecified k = 30 graph. The
order-swapped null behaved differently:
conditioning on within-donor effects produced an expected reversal fraction of
14.9%, and the observed 19.0% exceeded this expectation (one-sided empirical
P = 0.0001). The calibrated result is therefore conditional on donor-level
effects and is not a symmetric cross-unit test.

Actual-Milo graph and seed sensitivity qualified the strength of this
conclusion without removing the need for calibration. Expression matching
underestimated reversal in every graph, whereas the donor-conditioned residual
ranged from -1.8 to +1.8 percentage points and was not uniformly null. The
supported claim is that expression-only matching can overinterpret selected
panels and that effect-aware calibration is required. The stronger claim that
calibrated kidney excess is absent across all reasonable graph choices is not
supported. Graph construction and the Milo sampling seed should
be prespecified and reported alongside the conditioning direction.
The fixed-k proportion analysis varied prop from 0.05 to 0.20. This preserved
the positive donor-conditioned residual direction but changed its significance
decision, while the aTAL2 cell-rank conflict remained negative and significant
at all three proportions.
Embedding sensitivity was stronger: panel PCA30 reversed the sign of the
donor-conditioned panel residual and attenuated the aTAL2 cell-rank conflict
to non-significance. Graph construction, seed, proportion, and embedding
representation therefore all belong to the calibration specification.

The harmonized three-cohort benchmark made the direction dependence explicit.
Expression matching identified a target excess in all three cohorts, whereas
donor-conditioned calibration did so in none. Within-donor conditioning
retained an excess in KPMP and GSE254185 but not GSE211785. The same 3/3, 0/3,
and 2/3 pattern was observed at four thresholds, corresponding to 12/12,
0/12, and 8/12 significant cohort-threshold results. The harmonized benchmark
is stronger than a single-dataset sensitivity analysis. It shows that the
apparent target-panel excess and the effect-aware correction are reproducible,
but that the effect-aware conclusion is not symmetric across conditioning
directions.

The matched cross-tissue benchmark reinforced that interpretation. When all
cohorts were restricted to the same genes, the synovial high-inflammation
contrast remained significant under both conditioning directions at all four
thresholds, whereas the kidney donor-conditioned result remained non-significant.
Effect-aware calibration did not impose a universal
negative answer. It retained a context-specific non-kidney signal while
removing the apparent kidney panel excess.

An independent systemic sclerosis skin case-control cohort supplied the
opposite boundary. It showed no expression-matched donor-within excess and no
conditional excess in any direction, while an apparent neighborhood-within
excess under expression matching disappeared after effect-aware calibration.
The same framework therefore retained a synovial signal where one was present
and returned a null result where the independent case-control data did not
support one.

More generally, expression matching alone is insufficient for cross-unit
claims when selected genes differ in their marginal effect distributions. A
single effect-aware test is also insufficient unless the conditioning unit is
reported. The benchmark supports a dual-reporting standard: preserve
and report both donor-conditioned and within-donor-conditioned results.

The secondary unit pairs also showed context-dependent behavior across cohorts. Neighborhood-versus-
within retained conditional support in KPMP and GSE254185 but not GSE211785,
whereas donor-versus-neighborhood did not retain conditional support in the
primary cohorts. Analysis-unit disagreement is recurrent at the
gene level, but broad evidence for a particular pair is not universal. The
calibration therefore explains why selected panels can appear enriched; it does
not show that every kidney cross-unit signal is null.

### aTAL2 and frTAL exemplify the non-mergeability of analysis units

aTAL2 provides the clearest empirical case that donor-level disease elevation
and within-donor state coupling are distinct estimands. Donor pseudobulk and
neighborhood-expression association were positive, whereas the
disease-minus-Reference within-donor difference was negative in combined
disease and AKI. These opposing signs describe different dimensions of tissue
remodeling: average expression shifts between donors and local
cell-microenvironment relationships inside donors. Combining them into one
direction would discard information that is directly relevant to interpretation.
The negative aTAL2 cell-rank difference persisted across UMAP k, seed, and
proportion settings, but it attenuated to non-significance in panel PCA30.
The case is UMAP-graph-stable rather than embedding-invariant.

Composition-expression decomposition addresses the main compositional
alternative for the aggregate PAPPA2 difference. The subtype-expression
component dominated, whereas the composition component was small and its
interval included zero. The aggregate disease difference was not
explained by expansion of a high-PAPPA2 subtype in this decomposition. This
separates the source of the aggregate contrast from an inferred cell-intrinsic
mechanism.

frTAL provides a weaker but directionally consistent example: the within-donor
difference was negative across all three comparisons, with FDR support in the
combined comparison. Together, aTAL2 and frTAL show that cross-unit discordance
is graded across TAL subtypes rather than confined to one annotation state.
This graded pattern does not establish equivalence between the two annotations.

### Complete-data and donor-aware design provide orthogonal protection against common artifacts

The complete 106,851-cell TAL resource is central to the analysis. The earlier
4,446-cell subset preserved broad trends but changed directions for multiple
gene-level tests, especially in sparse subtypes. Coverage is an
analytic variable. The complete resource
provides the appropriate primary substrate and avoids carrying a known
coverage-dependent conclusion into the unit comparison.

The donor-aware design aligns resampling and inference with the biological unit
that contributes independent observations. Age- and sex-adjusted models,
leave-one-donor analysis, donor bootstrap, cortex/medulla stratification, and
neighborhood-purity adjustment tested different sources of dependence. The
donor-level main directions remained visible across these checks. The aTAL2
within-donor conflict was stable to donor resampling within UMAP analyses but
not to an alternative embedding. These analyses establish internal stability
only within the declared representation; external cohorts and experimental
systems address the separate question of generalization.

### External single-cell cohorts show cross-annotation consistency of the unit contrast

The two non-KPMP cohorts extend the donor-level PAPPA2 signal across different
disease etiologies and annotation systems. GSE211785 showed donor-level
elevation in both C_TAL and M_TAL, and GSE254185 showed donor-level elevation
with opposing neighborhood-expression and within-donor directions in injured
TAL states. These findings show that the donor-level direction extends beyond
KPMP. They also show that the direction of the unit contrast travels across
annotation systems, even where the within-donor estimate remains unresolved.

The external cohorts use different TAL labels: GSE211785 uses C_TAL and M_TAL,
whereas GSE254185 contains injured and inflammatory TAL states in obstructive
nephropathy. Their value is cross-annotation consistency of the directional
contrast and donor-level support. Exact subtype equivalence between these
annotations remains a separate question.

### An independent clinical cohort supports a PAPPA2-associated TAL-injury state

The PAPPA2-associated TAL-injury state showed an independent clinical
association with lower eGFR. In GSE175759, 87 primary kidney-disease samples
were analyzed. Higher TAL-injury score was associated with lower eGFR
after diagnosis adjustment (-8.71 points per SD; 95% CI -13.79 to -3.63;
P = 0.0012). The direction persisted among disease-only samples and after
library-size adjustment. The association also remained after adjustment for
TAL content and exceeded the effect generated by expression-matched random
11-gene sets. These analyses address two major alternatives: differences in
bulk TAL content and arbitrary selection of gene sets.

PAPPA2 alone carried a comparable clinical association (-9.06 points per SD;
P = 0.00067). The remaining 10 score genes produced a nominally significant
residual model, but the nested comparison did not demonstrate clear incremental
value (P = 0.089). The 11-gene score serves as a state-level summary, while
PAPPA2 alone reached -9.06 eGFR points per SD in this cohort. This evidence
remains limited to one independent clinical-association cohort.

### Kidney eQTL supports a compartment-specific regulatory context for PAPPA2

PAPPA2 has a compartment-specific regulatory context but no detected
disease-associated chromatin change. In the NephQTL2 resource reported by Han
et al. [18], PAPPA2 was among the strongest tubulointerstitial eGenes
(FDR = 1.64e-14), with no glomerular support (FDR = 0.526). This compartment
pattern aligns with the TAL anatomical origin and provides the clearest
regulatory evidence for PAPPA2 in this study. It supports a
tubulointerstitial regulatory context.

The same analysis did not provide PTWAS support for an eGFR association, and
the fixed aTAL2 marker set was not enriched for tubulointerstitial eGenes.
snATAC gene activity also did not show a disease-associated increase in PAPPA2
accessibility. These results bound the regulatory interpretation: PAPPA2 has a
strong compartment-specific eQTL context, whereas disease-specific chromatin
evidence and genetic causality remain unresolved.

### PAPPA2 as a benchmark case for unit-aware inference

PAPPA2 is the prespecified biological benchmark for the framework. Its donor-level
elevation was positive in all 12 subtype-comparison groups, and its direction
extended to non-KPMP cohorts. It also showed a tubulointerstitial eQTL context
and an association with lower eGFR in one independent clinical cohort. These
separate lines of evidence make PAPPA2 more than an arbitrary illustrative
gene: it is a cross-modally supported TAL-state candidate that survived
calibration at several scales.

The case illustrates the methodological contribution of the framework, which
requires explicit analysis units and separate reporting of each estimand. It
also requires paired expression-matched and effect-aware nulls that account for
expression depth and effect-size structure.
Candidate signals should then be validated at the scale at which they will be
interpreted. The workflow is implemented as `unitaware` v0.3.0 for reuse in
datasets with comparable exposure definitions.

PAPPA2 is best interpreted as a candidate TAL-state marker. The present data
support disease association and cross-modal context; biomarker,
therapeutic-target, and causal claims require perturbation and protein-level
validation. Its value here is to show that the framework can retain a
cross-modally supported signal at the scale supported by the data.

### Evidence boundaries define a higher standard for mechanism

Candidate pathway evidence did not meet a unified cross-unit standard. A
JAG1-NOTCH3-HES1 or EPHB2-STAT6 axis would need consistent direction across
subtypes, units, and comparisons, but this condition was not met. NOTCH3 was
significant only in CKD frTAL donor pseudobulk, HES1 had mixed directions,
STAT6 showed no significant donor-level change, and EPHB2 was subtype-dependent.
These results prevent isolated single-unit associations from being assembled
into a pathway model and leave the upstream regulators of PAPPA2 undefined.

This boundary clarifies the relationship to prior work. KPMP [5] and Cippà
et al. [19] established that human kidney tubules contain heterogeneous injury
and disease-associated states. The present findings agree on the existence of these states, but
show that their inferred direction can change with the analysis unit. The
conflict is with unit-agnostic interpretation, not with the states
themselves.

Oxvig [8], Dauber et al. [9], and Conover et al. [10] established PAPP-A2/PAPPA2
as a regulator of insulin-like growth factor availability and systemic growth.
The present data do not contest that established systemic endocrine role. They add a renal TAL-state context and show
that a systemic endocrine function should not be transferred directly to a
cell-state disease interpretation. A candidate should survive unit change,
calibrated nulls, cross-modal context, and validation at the intended
inferential scale. This standard defines the evidence required for mechanism-level interpretation.

### Boundary conditions and methodological limitations

The calibration framework rests on three assumptions and has a defined scope.
The model assumes that gene labels are exchangeable after conditioning on
marginal effect distributions. It also assumes comparable exposure definitions
across the units being compared. Each subtype-comparison stratum must contain
enough donors and cells for stable effect estimation. Under these conditions,
the estimand is
a residual reversal fraction after effect-size structure has been removed. The
framework explains calibration behavior within the tested panel, graph, and
representation; it does not identify causal cell-state transitions. It can fail
when strata are sparse, exposures are not comparable, or unmeasured
gene-specific structure violates exchangeability.

The discovery resource reused KPMP subtype labels and atlas-provided UMAP
coordinates. One AKI donor lacked neighborhood representation, and Milo used
logMS normalization. Alternative Milo graphs were directly evaluated at
k = 15, 30, and 50 and a second k = 30 seed. The aTAL2/PAPPA2 cell-rank
conflict remained directionally stable, but the residual after donor-conditioned
panel calibration did not. This makes graph construction and Milo sampling seed
part of the analysis specification. Framework interpretation also requires
comparable exposure definitions and adequate stratum size, and the conditional
pairing null assumes gene-label exchangeability after fixing marginal effects.
The null is directional: reversing the donor and within-donor units changes the
expected reversal fraction from 20.8% to 14.9%. The central empirical findings
rest on the complete KPMP resource and the prespecified calibration panel,
where analysis-unit non-equivalence and effect-aware calibration were measured
directly.

External cohorts use different TAL annotations and did not resolve within-donor
effects; the spatial resource shares KPMP provenance, and the external panel uses
a kNN approximation of neighborhood exposure. Sensitivity at k = 15, 30, and
50 preserved the external calibration pattern, but this remains an approximation
rather than a full external Milo replication. GSE175759 provides one
clinical-association cohort, without age or sex and with mixed diagnoses.
Nephroseq analyses are precomputed, and no perturbation or protein experiment
was performed, so these secondary data limitations constrain mechanistic and
clinical-utility claims.

Exact TAL states were also explored in the PREMIERE component of an integrated
adult kidney atlas. Among 31 PREMIERE donors whose labels did not occur in
KPMP, higher aTAL2 fraction tracked lower eGFR category (Spearman rho = -0.456;
FDR = 0.023). The association did not survive adjustment for diagnosis
(P = 0.115). Because low eGFR was concentrated in CKD, this analysis is a
disease-linked exploratory boundary and not an independent clinical validation
(Additional file 6).

### Implications for single-cell disease inference

The framework provides a reporting and validation template for single-cell
disease studies. An association is calibration-stable only when its direction
remains consistent across prespecified analysis units, both conditioning
directions, and a declared graph/seed/embedding ensemble. Its calibrated
residual must also not cross the prespecified decision boundary.
Studies should report:

1. the analysis unit and the independent biological replicate;
2. every unit-specific estimand separately, rather than treating agreement as validation;
3. the expression-matched and effect-aware null results in both conditioning directions;
4. the residual range across prespecified graph, seed, proportion, and embedding variants;
5. validation at the same inferential scale as the intended claim.

For rare cell subtypes, coverage and effect-size selection can materially alter
interpretation. It also prevents
a single graph-specific or seed-specific calibration result from being reported
as a general biological conclusion.

The next evidence threshold is exact subtype replication in an independent
kidney cohort, continuous eGFR outcomes, and perturbation or protein-level
analysis of the PAPPA2-associated state (Fig. 9). These studies will test the
biological mechanism while preserving the methodological conclusion established
here.

[[FIGURE:Fig9]]

## Conclusions

The new knowledge from this study is that single-cell disease associations
change meaning across donor, neighborhood, and within-donor units. Effect-aware
calibration plus graph/embedding reporting is required before those results are
combined as one disease effect. The framework treats donor pseudobulk,
neighborhood abundance, neighborhood-expression association, and within-donor
cell rank as distinct estimands. It applies effect-aware nulls and reports
graph, seed, proportion, and embedding stability as part of the inference.

The strongest evidence comes from the calibration benchmark and the PAPPA2 case.
Expression matching produced a target-panel excess, whereas donor-conditioned
calibration did so in none of the three kidney cohorts (3/3 versus 0/3). The
calibrated panel residual varied across Milo graphs, proportions, seeds, and
embeddings. aTAL2/PAPPA2 retained a positive donor-level signal and a negative
within-donor direction across UMAP settings, while panel PCA30 attenuated the
within-donor conflict. PAPPA2 also showed positive donor elevation in all 12
subtype comparisons, external directional support, a tubulointerstitial eQTL
context, and association with lower eGFR in one independent cohort.

The framework supports calibration and association, not causal inference. Its
conclusions apply to the tested KPMP TAL resource, prespecified panel, graph
and embedding variants, and external cohorts with comparable exposures.
PAPPA2 is a candidate TAL-state marker; causal, biomarker, therapeutic, and
clinical-utility claims require perturbation, protein-level, and further cohort
validation.

## List of abbreviations

AI: artificial intelligence

AKI: acute kidney injury

API: application programming interface

BH: Benjamini-Hochberg

CI: confidence interval

CKD: chronic kidney disease

CLI: command-line interface

CPM: counts per million

C-TAL-A, C-TAL-B: KPMP human kidney atlas cortical thick ascending limb
subtype labels

C_TAL: cortical thick ascending limb label in GSE211785

DKD: diabetic kidney disease

DOI: digital object identifier

eGFR: estimated glomerular filtration rate

eQTL: expression quantitative trait locus

FDR: false discovery rate

frTAL: KPMP human kidney atlas thick ascending limb subtype label

GEO: Gene Expression Omnibus

GFR: glomerular filtration rate

H5AD: Hierarchical Data Format 5 AnnData

HC3: heteroskedasticity-consistent type 3

HKD: hypertensive kidney disease

KPMP: Kidney Precision Medicine Project

logFC: log2 fold change

MSigDB: Molecular Signatures Database

M_TAL: medullary thick ascending limb label in GSE211785

NephQTL2: human kidney eQTL resource

OR: odds ratio

PTWAS: probabilistic transcriptome-wide association study

SD: standard deviation

snATAC-seq: single-nucleus assay for transposase-accessible chromatin sequencing

snRNA-seq: single-nucleus RNA sequencing

SpatialFDR: spatial false discovery rate

TAL: thick ascending limb

UMAP: uniform manifold approximation and projection

## Declarations

### Ethics approval and consent to participate

This study involved computational analysis only. It used publicly available or
controlled-access, de-identified human kidney data. No new human samples were
collected, no participants were recruited, and no clinical intervention was
performed. No human or animal experiments were conducted, and no experimental
protocol or approval number applies. Ethical approval and informed consent were
not required for this secondary analysis; data access and use were governed by
the terms of the original studies and repositories.

### Consent for publication

Not applicable. No identifiable individual-level data are presented.

### Clinical trial number

No clinical trial was conducted; therefore, no clinical trial registration
number applies.

### Availability of data and materials

The datasets supporting the conclusions of this article are available in
the Gene Expression Omnibus repository: GSE183279 [20]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183279), which
includes single-nucleus RNA sequencing data (GSE183277) and spatial
transcriptomics data [21]; GSE137570 [22]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137570); GSE211785
[23,24]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785); GSE254185
[25,26]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185); and GSE175759
[16]
(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175759), which includes
public sample-level eGFR and processed expression matrices. The frozen
`TAL_Injury_h1` marker list used to construct the score was obtained from
[17]; GSE284715/GSE284742 were not reanalyzed in this study.
GSE243917 synovial single-nucleus RNA sequencing data [27] are available at
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243917.
GSE138669 skin case-control single-cell RNA sequencing data [28] are
available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138669.
NephQTL2 gene-level and PTWAS summary data were obtained from the
supplementary information of [18]
(doi:10.1038/s41467-023-37691-7). Nephroseq v5 was queried through
its public interface [13]
(https://www.nephroseq.org).

Access to the Gene Expression Omnibus datasets is unrestricted, including
for non-academic users. Nephroseq v5 is accessible through its public web
interface, and its use is subject to the terms of use published on the
Nephroseq website [13].

The processed and derived datasets supporting the conclusions of this article
that are not available in public repositories are included within this article
and its Additional files 1-6 and Supplementary Figs. S1-S3. The full
within-donor permutation matrix and other large source tables are retained in
the Zenodo archive to keep each submission attachment below the journal file
size limit. Analysis data are released under CC BY 4.0.

### Code availability

All analysis code is available at
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity. The
`unitaware` v0.3.0 source and calibration-stability scripts are available at
https://github.com/ZhiPeng1022/pappa2-tal-analysis-unit-sensitivity/tree/main/unitaware
and are included in the submission package. The archived v0.3.0 release is
available at https://doi.org/10.5281/zenodo.24441900. The archived v0.2.0
baseline is available at https://doi.org/10.5281/zenodo.24042041.

The code runs on a platform-independent Python environment (Python 3.10 or
higher; numpy, pandas, and scipy) and uses R for the Milo analyses. Code is
released under the MIT license. No restrictions apply to non-academic use.

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

### Use of artificial intelligence

OpenAI Codex CLI was used to assist with literature summarization, code
organization, and English-language editing. All AI-assisted outputs were checked
against source data and references, reviewed, and edited by the authors, who
retain full responsibility for the manuscript. No scientific data, results, or
conclusions were generated by AI.

### Acknowledgements

The authors thank the participants and investigators of the Kidney Precision
Medicine Project, NephQTL2, Nephroseq, and the GEO datasets used in this study
for generating and sharing the data.

## Additional files

Additional file 1. File format: Microsoft Excel Worksheet (.xlsx). Title:
Analysis-unit calibration and frozen panel. Description: Frozen 400-gene panel
and 2,000-gene control pool; gene-level donor, neighborhood, and within-donor
metrics; expression-matched random-set results; threshold sensitivity;
effect-aware conditional-null summary and runs; and per-test summaries of the
within-donor permutation nulls. The full within-donor permutation matrix is
retained in the Zenodo archive.

Additional file 2. File format: Microsoft Excel Worksheet (.xlsx). Title:
Single-cell neighborhood and within-donor analyses. Description: Complete-data
Milo neighborhoods, subset-level Augur separability, complete-versus-subset
pseudobulk comparison, neighborhood-expression robustness and bootstrap
replicates, joint analysis-unit concordance, within-donor cell-rank results,
donor pseudobulk models, leave-one-donor neighborhood results, and
composition-expression decomposition.

Additional file 3. File format: Microsoft Excel Worksheet (.xlsx). Title:
Unit-aware simulation benchmark. Description: Per-replicate simulation runs and
scenario-level summaries for the six prespecified effect structures.

Additional file 4. File format: Microsoft Excel Worksheet (.xlsx). Title:
PAPPA2 clinical and composition analyses. Description: KPMP ordinal eGFR
models; GSE137570 directional association; Nephroseq precomputed GFR
associations; and the locked GSE175759 sample scores, gene coverage, primary
and sensitivity models, TAL-content adjustment, expression-matched nulls, and
PAPPA2 incremental models.

Additional file 5. File format: Microsoft Excel Worksheet (.xlsx). Title:
External validation and regulatory context. Description: Non-KPMP donor-level
results, frozen marker transfer and held-out metrics, external aTAL2-like
assignments, NephQTL2 gene-level eQTL summaries, GSE211785 snATAC summaries,
the k = 15/30/50 neighborhood sensitivity, the harmonized three-cohort
calibration benchmark, cohort and unit-pair results, method comparison,
donor-within threshold sensitivity, and the matched cross-tissue benchmark,
including synovial donor mapping and threshold results. It also includes the
independent GSE138669 skin case-control boundary analysis, donor composition,
and unit-pair results.

Additional file 6. File format: Microsoft Excel Worksheet (.xlsx). Title:
Graph, seed, and exact-state boundary analyses. Description: Actual-Milo
KPMP panel calibration at k = 15, 30, and 50 plus a second k = 30 seed; the
formal calibration-stability classification table; fixed-k neighborhood
proportion sensitivity at prop = 0.05, 0.10, and 0.20; UMAP-versus-PCA30
embedding sensitivity; per-run Milo provenance (cells, donors, neighborhoods,
significance counts, and seeds); full-TAL aTAL2 PAPPA2 cell-rank sensitivity;
graph-invariant PAPPA2 composition-expression decomposition; and exploratory
PREMIERE exact TAL-state eGFR-category associations before and after diagnosis
adjustment.

## Supplementary figures

Supplementary Fig. S1. File format: Tagged Image File Format (.tif). Title:
Analysis-unit robustness. Description: Expression-matched and effect-aware
calibration summaries, including the 400-gene target, the 300-gene biological
panel, and standardized-effect thresholds.

Supplementary Fig. S2. File format: Tagged Image File Format (.tif). Title:
Unit-aware simulation benchmark. Description: Mean donor and within-donor
effects and donor-within reversal fractions across the six simulation
scenarios.

Supplementary Fig. S3. File format: Tagged Image File Format (.tif). Title:
Matched cross-tissue calibration benchmark. Description: Donor-within
reversal fractions and empirical P values in three kidney cohorts and
GSE243917 synovial fibroblasts after restriction to the same target and
control genes. The synovial contrast is high versus moderate inflammation, not
a kidney disease-control comparison.

## References

1. Bikbov B, Purcell CA, Levey AS, Smith M, Abdoli A, Abebe M, Adebayo OM, Afarideh M, Agarwal SK, Agudelo-Botero M, et al. Global, regional, and national burden of chronic kidney disease, 1990–2017: a systematic analysis for the Global Burden of Disease Study 2017. Lancet. 2020;395:709–33.

2. Hoste EAJ, Kellum JA, Selby NM, Zarbock A, Palevsky PM, Bagshaw SM, Goldstein SL, Cerdá J, Chawla LS. Global epidemiology and outcomes of acute kidney injury. Nat Rev Nephrol. 2018;14:607–25.

3. Ferenbach DA, Bonventre JV. Mechanisms of maladaptive repair after AKI leading to accelerated kidney ageing and CKD. Nat Rev Nephrol. 2015;11:264–76.

4. Bonventre JV, Yang L. Cellular pathophysiology of ischemic acute kidney injury. J Clin Invest. 2011;121:4210–21.

5. Lake BB, Menon R, Winfree S, Hu Q, Melo Ferreira R, Kalhor K, Barwinska D, Otto EA, Ferkowicz M, Diep D, et al. An atlas of healthy and injured cell states and niches in the human kidney. Nature. 2023;619:585–94.

6. Knepper MA, Kim GH, Fernández-Llama P, Ecelbarger CA. Regulation of thick ascending limb transport by vasopressin. J Am Soc Nephrol. 1999;10:628–34.

7. Brezis M, Rosen S, Silva P, Epstein FH. Selective vulnerability of the medullary thick ascending limb to anoxia in the isolated perfused rat kidney. J Clin Invest. 1984;73:182–90.

8. Oxvig C. The role of PAPP-A in the IGF system: location, location, location. J Cell Commun Signal. 2015;9:177–87.

9. Dauber A, Muñoz-Calvo MT, Barrios V, Domené HM, Kloverpris S, Serra-Juhé C, Desikan V, Pozo J, Muzumdar R, Martos-Moreno GÁ, et al. Mutations in pregnancy-associated plasma protein A2 cause short stature due to low IGF-I availability. EMBO Mol Med. 2016;8:363–74.

10. Conover CA, Boldt HB, Bale LK, Clifton KB, Grell JA, Mader JR, Mason EJ, Powell DR. Pregnancy-associated plasma protein-A2 (PAPP-A2): tissue expression and biological consequences of gene knockout in mice. Endocrinology. 2011;152:2837–44.

11. Asghari M, Sabo AR, Barwinska D, Melo Ferreira R, Ferkowicz MJ, Bowen WS, Cheng YH, Gisch DL, Gulbronson C, Phillips CL, et al. Integration of spatial protein imaging and transcriptomics in the human kidney tracks the regenerative potential of proximal tubules. Sci Adv. 2025;11:eadv8918.

12. Kennedy C, Doyle R, Gough O, McEvoy C, McAnallen S, Hughes M, Sheng X, Crifo B, Andrews D, Gaffney A, et al. A novel role for FERM domain-containing protein 3 in CKD. Kidney360. 2024;5:1799–812.

13. Nephroseq v5. University of Michigan, Ann Arbor, MI, USA. https://www.nephroseq.org. Accessed September 29, 2026.

14. Skinnider MA, Squair JW, Kathe C, Anderson MA, Gautier M, Matson KJE, Milano M, Hutson TH, Barraud Q, Phillips AA, et al. Cell type prioritization in single-cell data. Nat Biotechnol. 2021;39:30–4.

15. Dann E, Henderson NC, Teichmann SA, Morgan MD, Marioni JC. Differential abundance testing on single-cell data using k-nearest neighbor graphs. Nat Biotechnol. 2022;40:245–53.

16. Park S, Lee H, Lee J, Lee S, Cho S, Huh H, Kim JY, Park M, Lee S, Kim Y, et al. RNA-seq profiling of tubulointerstitial tissue reveals a potential therapeutic role of dual anti-phosphatase 1 in glomerulonephritis. J Cell Mol Med. 2022;26:3364–77. Dataset: GEO GSE175759. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175759. Accessed November 8, 2026.

17. Pfefferkorn AM, Jahn L, Gauthier PT, Kulow VA, Roeles J, Müller-Bötticher N, Gerhardt LMS, Leiz J, Sarfraz S, Plumbom I, et al. Injured epithelial cell states impact kidney allograft survival after T-cell-mediated rejection. Nat Commun. 2026;17:1060.

18. Han SK, McNulty MT, Benway CJ, Wen P, Greenberg A, Onuchic-Whitford AC, Nephrotic Syndrome Study Network (NEPTUNE), Jang D, Flannick J, Burtt NP, et al. Mapping genomic regulation of kidney disease and traits through high-resolution and interpretable eQTLs. Nat Commun. 2023;14:2229.

19. Cippà PE, Sun B, Liu J, Chen L, Naesens M, McMahon AP. Transcriptional trajectories of human kidney injury progression. JCI Insight. 2018;3:e123151.

20. Kidney Precision Medicine Project. Single-nucleus RNA sequencing resource, GEO superseries GSE183279. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183279. Accessed September 29, 2026.

21. Kidney Precision Medicine Project. Spatial transcriptomics resource, GEO GSE183456. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183456. Accessed September 29, 2026.

22. Human kidney transcriptome cohort, GEO GSE137570. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137570. Accessed September 29, 2026.

23. Abedini A, Levinsohn J, Klötzer KA, Dumoulin B, Ma Z, Frederick J, Dhillon P, Balzer MS, Shrestha R, Liu H, et al. Single-cell multi-omic and spatial profiling of human kidneys implicates the fibrotic microenvironment in kidney disease progression. Nat Genet. 2024;56:1712–24.

24. Human kidney multi-omic atlas, GEO GSE211785. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211785. Accessed October 28, 2026.

25. Reck M, Baird DP, Veizades S, Sutherland C, Bell RMB, Hur H, Cairns C, Janas PP, Campbell R, Nam A, et al. Multiomic analysis of human kidney disease identifies a tractable inflammatory and pro-fibrotic tubular cell phenotype. Nat Commun. 2025;16:4745.

26. Human kidney multiome resource with obstructive nephropathy, GEO GSE254185. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254185. Accessed October 28, 2026.

27. Gupta A, Weinand K, Nathan A, Sakaue S, Zhang MJ, Accelerating Medicines Partnership RA/SLE Program and Network, Donlin L, Wei K, Price AL, Amariuta T, et al. Dynamic regulatory elements in single-cell multimodal data implicate key immune cell states enriched for autoimmune disease heritability. Nat Genet. 2023;55:2200–10.

28. Human skin single-cell RNA sequencing resource for systemic sclerosis, GEO GSE138669. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138669. Accessed November 21, 2026.

## Figure legends

**Fig. 1. Study design and analysis-unit evaluation.**
The upper row summarizes the complete Kidney Precision Medicine Project (KPMP)
single-nucleus RNA sequencing (snRNA-seq) thick ascending limb (TAL) resource
(106,851 cells from 223 donors), the 4,446-cell derivative subset used only as
a comparator, and the prespecified audit based on a 400-gene panel, matched
controls, and an effect-aware conditional null. The lower row defines four
analysis units: donor pseudobulk, Milo neighborhood abundance,
neighborhood-expression association, and within-donor cell-rank association.
Each unit is shown with the question it answers and its principal limitation.
All units use the same cells and embedding. Donor-level resampling is a
stability check, and PAPPA2 is the prespecified case gene.

**Fig. 2. Subtype-specific disease-related neighborhood remodeling.**
Complete-data Milo results evaluating the four TAL subtypes across combined
disease, acute kidney injury (AKI), and chronic kidney disease (CKD) versus
Reference. a, Fraction of annotated neighborhoods called significant per
subtype. b, Median neighborhood log2 fold change (logFC). c, Number of
significant neighborhoods.
Neighborhoods were defined using atlas uniform manifold approximation and
projection (UMAP) coordinates and tested with Milo using logMS normalization,
a fixed seed, and a spatial false discovery rate (SpatialFDR) threshold of
0.05. One AKI donor without neighborhood representation was excluded.

**Fig. 3. Conceptual comparison of the four analysis units.**
Schematic comparison of the question answered and principal interpretation
boundary for each analysis unit. Donor pseudobulk can be influenced by donor
composition. Milo neighborhoods overlap. Neighborhood-expression association
is a cross-neighborhood relation and does not resolve within-donor coupling.
Within-donor cell-rank association describes local state coupling without
implying causality. Donor-level robustness analyses are stability checks, not a
fifth analysis unit.

**Fig. 4. Coverage loss and directional consequences in the derivative TAL subset.**
Comparison of the complete 106,851-cell TAL resource with the 4,446-cell
derivative subset. a, Cells per subtype (log scale). b, Number of Reference
donors per subtype. c, Agreement or conflict in pseudobulk direction across 60
gene-subtype-comparison tests. d, Number of discordant tests by gene.
The severe depletion of rare subtypes, including aTAL2 and frTAL, is evident
and accompanies the 12 directional conflicts.

**Fig. 5. PAPPA2 case analysis across TAL subtypes.**
PAPPA2 results are arranged by comparison (rows) and TAL subtype (columns).
a, Mean disease-minus-Reference donor pseudobulk difference. b, Age- and
sex-adjusted linear-model beta. c, Spearman correlation across Milo
neighborhoods between mean PAPPA2 expression and neighborhood logFC.
d, Disease-minus-Reference difference in the median within-donor Spearman
correlation between PAPPA2 expression and neighborhood exposure, using donors
with at least five cells. Positive values in d indicate a higher within-donor
correlation in disease donors than in Reference donors. Panels a and b are
donor-level estimates; panels c and d represent different analysis units and
are not independent replications.

**Fig. 6. PAPPA2 cross-unit conflict in aTAL2 and subtype-wide conflict counts.**
a, aTAL2 donor pseudobulk mean expression difference for PAPPA2 across the
three comparisons. b, aTAL2 disease-minus-Reference difference in the median
within-donor cell-rank correlation across the same comparisons. The aTAL2
pattern combines positive donor-level association with a negative
within-donor difference. c, Number of conflicting gene-comparison tests for
any pair among donor pseudobulk, neighborhood-expression, and within-donor
units across five genes and three comparisons, shown for each subtype. The
panel scale is 15 possible tests per subtype. Conflict count alone does not
distinguish aTAL2.

**Fig. 7. Internal robustness of the PAPPA2 TAL pseudobulk effect.**
PAPPA2 donor pseudobulk effect across 12 subtype-comparison groups (four
subtypes by three disease comparisons). a, Full-data estimate,
leave-one-donor range, and donor-stratified bootstrap median with its 95%
interval. b, Descriptive cortex and medulla mean disease-minus-Reference
difference for the combined-disease comparison. Panel b is not an inferential
test.

**Fig. 8. Calibration stability and portability across cohorts, graph
settings, and embeddings.**
a, Observed donor-within reversal proportion compared with the
expression-matched control-pool expectation and the donor-conditioned
effect-aware expectation in KPMP, GSE211785, and GSE254185. Whiskers show the
2.5th-97.5th percentile of the corresponding null distributions. b, KPMP
donor-conditioned calibration residual in percentage points across UMAP k
values, a second k = 30 seed, and fixed-k neighborhood proportions. Filled and
open points denote P < 0.05 and P >= 0.05, respectively; asterisks and ns
summarize the same tests. c,d, Comparison of atlas UMAP2D and panel PCA30 at
fixed k = 30, prop = 0.10, and seed 20261122: calibrated panel residual (c) and
aTAL2/PAPPA2 within-donor cell-rank effect (d). Exact estimates and stability
classes are provided in Tables 1 and 2 and Additional file 6.

**Fig. 9. Subtype-specific TAL state model.**
Schematic summary of supported and unsupported interpretations. The upper part
shows subtype-specific thick ascending limb (TAL) neighborhood states, the
unit-aware analysis framework, a PAPPA2-high TAL state, and its association
with lower estimated glomerular filtration rate (eGFR) in one independent
cohort. The lower part distinguishes supported observations from unsupported
interpretations, including causality, stable multicohort clinical utility, and
a unified JAG1-NOTCH3-HES1 or EPHB2-STAT6 mechanism. The arrow between PAPPA2
and eGFR denotes association, not a causal direction or established clinical
utility.
