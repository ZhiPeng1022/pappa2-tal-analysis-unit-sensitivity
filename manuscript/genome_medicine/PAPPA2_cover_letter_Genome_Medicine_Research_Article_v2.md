29 October 2026

Dear Editors,

**Article type.** We submit our manuscript, "Unit-aware analysis reveals opposing donor and within-donor disease signals in human kidney single-cell data", for consideration as a Research Article.

**The problem.** Single-cell disease studies often describe a cell state as disease-associated without specifying whether the signal reflects donor composition, local neighborhood structure, or coupling within a donor. These analysis units are commonly treated as interchangeable, although they answer different biological questions.

**What we developed.** We developed unitaware v0.2.0, a reusable framework that reports donor pseudobulk, neighborhood-expression association, and within-donor cell-rank association separately, standardizes direction reversal, decomposes donor differences into composition and expression effects, and calibrates observed signals with both expression-matched and effect-aware nulls. The implementation, command-line workflow, six unit tests, documentation, and MIT license accompany the manuscript.

**What we found.** In 106,851 KPMP human kidney thick ascending limb cells, donor-level and within-donor analyses could yield opposite subtype conclusions. In a prespecified 400-gene panel, donor-within reversal occurred in 19.0% of target genes versus 12.7% of expression-matched random sets. An effect-aware pairing null gave 20.8%, showing that the initial excess reflected the selected effect-size distribution rather than broad negative coupling. A smaller neighborhood-versus-within excess remained threshold sensitive. A simulation benchmark recovered true opposing effects in 95.5-97.3% of replicates, compared with 2.3% under a concordant scenario.

**Independent context.** Two non-KPMP kidney cohorts supported donor-level PAPPA2 disease elevation, but did not provide exact subtype replication or a significant within-donor effect. Kidney tubulointerstitial eQTL and snATAC analyses are presented as regulatory context and boundary evidence rather than as causal validation.

**Why Genome Medicine.** The manuscript addresses a general measurement problem in genomic medicine rather than a single-gene association. It combines a large human kidney resource, external cohorts, expression-matched and effect-aware null calibration, reproducible software, and a simulation benchmark. The negative result for broad panel-level coupling is reported explicitly and strengthens the methodological conclusion. The resulting framework can be applied to other tissues and disease states, and we believe it fits the journal's broad readership in genomic medicine and precision health.

**Declarations.** The manuscript is original, has not been published, and is not under consideration elsewhere. All authors have approved the submission. The authors declare no competing interests and no specific funding. The study used de-identified public or controlled-access data; no new human samples were collected. No clinical trial was performed. The manuscript is not intended for a special issue.

**APC.** The authors have no publication fund. We request consideration for a discretionary full APC waiver or the maximum available discount at submission, in accordance with the journal policy for case-by-case financial need.

We have not requested the exclusion of any reviewer.

Thank you for considering our manuscript.

Yours sincerely,

Peng Zhi

Department of General Practice, Henan Provincial Chest Hospital

Zhengzhou, Henan, China

Email: zhiyoyo0409@163.com

ORCID: 0009-0003-2201-5131

On behalf of all authors
