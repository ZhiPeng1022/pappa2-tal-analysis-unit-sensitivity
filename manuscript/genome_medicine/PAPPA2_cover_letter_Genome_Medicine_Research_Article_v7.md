22 November 2026

Dear Editors,

**Article type.** We submit our manuscript, "Unit-aware calibration of single-cell disease associations: analysis-unit, effect-size, and graph-stability audits in human kidney atlases", for consideration as a Research Article.

**The problem.** Single-cell disease studies often describe a cell state as disease-associated without specifying whether the signal reflects donor composition, local neighborhood structure, or coupling within a donor. These analysis units are commonly treated as interchangeable, although they answer different biological questions.

**What we developed.** We developed unitaware v0.3.0, a reusable framework that reports donor pseudobulk, neighborhood-expression association, and within-donor cell-rank association as complementary estimands, standardizes direction reversal, decomposes donor differences into composition and expression effects, and calibrates observed signals with both expression-matched and effect-aware nulls. The implementation, command-line workflow, ten unit tests, documentation, and MIT license accompany the manuscript.

**Central contribution.** The paper defines calibration stability as a reporting standard. Analysis unit, expression matching, effect-aware conditioning direction, and graph/seed/embedding sensitivity are reported together, so a single graph-specific null is not mistaken for a general biological conclusion. PAPPA2 is a benchmark showing how the framework separates donor-level elevation from representation-sensitive within-donor coupling; it is not the sole finding.

**What we found.** In 106,851 KPMP human kidney thick ascending limb cells, donor-level and within-donor analyses could yield opposite subtype conclusions. In a prespecified 400-gene panel, donor-within reversal occurred in 19.0% of target genes versus 12.7% of expression-matched random sets. An effect-aware pairing null gave 20.8% in the prespecified graph, showing that expression matching underestimated effect-size-driven reversal. Actual-Milo graph and seed sensitivity then showed that the residual after calibration was not uniform: the aTAL2/PAPPA2 cell-rank conflict remained negative across k = 15, 30, and 50 and a second k = 30 seed, but it attenuated in a panel PCA30 graph while the panel-level residual changed direction. The supported methodological claim is therefore calibration-dependent rather than a universal negative result. In a harmonized three-cohort benchmark, expression matching was significant in 3/3 cohorts, donor-conditioned calibration in 0/3, and within-donor-conditioned calibration in 2/3; this direction contrast remained stable across four thresholds. In a matched-gene non-kidney synovial benchmark, the target signal remained significant under both conditioning directions, showing that calibration was not uniformly suppressive. An independent systemic sclerosis skin case-control cohort showed no donor-within expression-matched or calibrated excess, providing the expected null boundary. A simulation benchmark recovered true opposing effects in 95.5-97.3% of replicates, compared with 2.3% under a concordant scenario. In an independent 87-sample primary kidney-disease cohort, the frozen 11-gene score was associated with lower eGFR after diagnosis adjustment.

**Independent context.** Two non-KPMP kidney cohorts supported donor-level PAPPA2 disease elevation, but did not provide exact subtype replication or a significant within-donor effect. Kidney tubulointerstitial eQTL and snATAC analyses are presented as regulatory context and boundary evidence rather than as causal validation. The GSE175759 clinical association is reported as one independent cohort, not as stable multicohort validation or clinical utility. Actual-Milo graph, seed, proportion, and PCA30 embedding sensitivity is restricted to KPMP; external cohorts test portability of the calibration procedure, not full graph replication.

**Why Genome Medicine.** The manuscript addresses a general measurement problem in genomic medicine rather than a single-gene association. It combines a large human kidney resource, matched kidney and non-kidney method benchmarks, an independent primary-kidney-disease clinical-association cohort, expression-matched and effect-aware null calibration, reproducible software, and a simulation benchmark. Graph- and seed-dependent calibration limits are reported explicitly and strengthen the methodological conclusion. The intended readership is computational biology, kidney atlas research, and translational teams selecting associations for independent validation.

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
