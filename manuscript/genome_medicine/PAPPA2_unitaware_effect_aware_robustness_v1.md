# PAPPA2 unitaware effect-aware robustness v1

Date: 2026-10-29

## Purpose

The original analysis-unit panel compared the 400-gene target panel with
expression-matched random gene sets. Expression matching does not match the
two effect-size distributions that determine whether an opposite-sign pair
passes the standardized-effect threshold. We therefore added a conditional
pairing permutation that preserves:

- the target gene's left-unit effect;
- the marginal distribution of the right-unit effect within each subtype and
  comparison;
- all missingness and group-validity rules.

Only the pairing between the two units was broken. The test used 10,000
permutations, a seed of `20261028`, and gene-balanced aggregation matching the
original panel endpoint.

## Main result

At the prespecified `|z| >= 0.5` threshold:

| Unit pair | Target panel | Expression-matched sets | Conditional null | Conditional P |
|---|---:|---:|---:|---:|
| donor vs within-donor | 19.0% | 12.7% | 20.8% | 0.999 |
| donor vs neighborhood-expression | 21.9% | 14.5% | 21.2% | 0.111 |
| neighborhood-expression vs within-donor | 21.0% | 19.6% | 19.4% | 0.0017 |

The expression-matched comparison shows that the selected panel contains more
large, opposite-sign pairs than randomly selected genes with similar mean
expression. The conditional test shows that this excess for donor versus
within-donor analysis is explained by the joint effect-size distribution:
observed pairing was not more discordant than independent re-pairing of the
same effects.

The neighborhood-expression versus within-donor pair retained a small
conditional excess. This result was attenuated when the analysis was limited
to the 300 biological genes (`P = 0.189` with expression-matched sets) and was
not significant when the effect threshold was increased to 0.75 or higher.

## Threshold sensitivity

The expression-matched excess for donor versus within-donor and donor versus
neighborhood-expression remained at every tested threshold from 0.5 to 2.0.
That pattern is consistent with selection of genes that have larger effects
across units, not with a threshold-specific artefact.

Neighborhood-expression versus within-donor did not remain significant above
the 0.5 threshold. It should therefore be described as a small, threshold-
sensitive excess rather than a general feature of the panel.

## Consequence for the manuscript

The following claims are not supported and must be removed:

- that donor-within reversal is enriched beyond an effect-aware null;
- that the 400-gene panel establishes a broad negative coupling between
  donor-level and within-donor disease effects;
- that the null simulation is irrelevant because it differs from the random-set
  comparison.

The following claims remain supportable:

- analysis units can disagree in direction for individual genes and states;
- aTAL2 remains a concrete cross-unit discordance example;
- expression-only random sets can overstate unit-specific coupling when they
  do not match effect-size structure;
- effect-aware calibration should accompany future unit-aware analyses.

## Outputs

Directory:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\outputs\unitaware_robustness
```

Files:

- `threshold_sensitivity.csv`
- `pairing_permutation_summary.csv`
- `pairing_permutation_null.csv`
- `robustness_manifest.json`

The analysis code is:

```text
C:\Users\Elsa\Documents\Codex\2026-10-15\ji\work\run_unitaware_robustness.py
```
