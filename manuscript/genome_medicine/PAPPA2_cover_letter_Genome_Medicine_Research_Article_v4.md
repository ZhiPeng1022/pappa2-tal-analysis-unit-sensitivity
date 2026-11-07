8 November 2026

Dear Editors,

**Article type.** We submit our manuscript, "Analysis-unit dependence and effect-aware calibration in human kidney thick ascending limb single-cell studies", for consideration as a Research Article.

**The problem.** Single-cell disease studies often describe a cell state as disease-associated without specifying whether the signal reflects donor composition, local neighborhood structure, or coupling within a donor. These analysis units are commonly treated as interchangeable, although they answer different biological questions.

**What we developed.** We developed unitaware v0.2.0, a reusable framework that reports donor pseudobulk, neighborhood-expression association, and within-donor cell-rank association as complementary estimands, standardizes direction reversal, decomposes donor differences into composition and expression effects, and calibrates observed signals with both expression-matched and effect-aware nulls. The implementation, command-line workflow, six unit tests, documentation, and MIT license accompany the manuscript.

**What we found.** In 106,851 KPMP human kidney thick ascending limb cells, donor-level and within-donor analyses could yield opposite subtype conclusions. In a prespecified 400-gene panel, donor-within reversal occurred in 19.0% of target genes versus 12.7% of expression-matched random sets. An effect-aware pairing null gave 20.8%, showing that the initial excess reflected the selected effect-size distribution rather than broad negative coupling. A smaller neighborhood-versus-within excess remained threshold sensitive. The same expression-matched versus effect-aware divergence recurred in two external cohorts. A simulation benchmark recovered true opposing effects in 95.5-97.3% of replicates, compared with 2.3% under a concordant scenario. In an independent 87-sample primary kidney-disease cohort, the frozen 11-gene score was associated with lower eGFR after diagnosis adjustment.

**Independent context.** Two non-KPMP kidney cohorts supported donor-level PAPPA2 disease elevation, but did not provide exact subtype replication or a significant within-donor effect. Kidney tubulointerstitial eQTL and snATAC analyses are presented as regulatory context and boundary evidence rather than as causal validation. The GSE175759 clinical association is reported as one independent cohort, not as stable multicohort validation or clinical utility.

**Why Genome Medicine.** The manuscript addresses a general measurement problem in genomic medicine rather than a single-gene association. It combines a large human kidney resource, external single-cell cohorts, an independent primary-kidney-disease clinical-association cohort, expression-matched and effect-aware null calibration, reproducible software, and a simulation benchmark. The negative result for broad panel-level coupling is reported explicitly and strengthens the methodological conclusion. The intended readership is computational biology, kidney atlas research, and translational teams selecting associations for independent validation.

**Editorial assessment.** We would be grateful for the editors' assessment of whether the manuscript meets the journal's requirements and merits external peer review.

**Declarations.** The manuscript is original, has not been published, and is not under consideration elsewhere. All authors have approved the submission. The authors declare no competing interests and no specific funding. The study used de-identified public or controlled-access data; no new human samples were collected. No clinical trial was performed. The manuscript is not intended for a special issue.

We have not requested the exclusion of any reviewer.

Thank you for considering our manuscript.

Yours sincerely,

Peng Zhi

Department of General Practice, Henan Provincial Chest Hospital

Zhengzhou, Henan, China

Email: zhiyoyo0409@163.com

ORCID: 0009-0003-2201-5131

On behalf of all authors
