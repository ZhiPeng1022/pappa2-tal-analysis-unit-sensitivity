# PAPPA2 unit-aware simulation results v1

Date: 2026-10-29

## Purpose

Test whether the unit-aware framework distinguishes genuine opposing donor and
within-donor effects from null, donor-only, within-only, and concordant
scenarios.

The simulation is implemented in the standalone `unitaware` package and uses:

- 25 Reference and 25 disease donors per replicate;
- 30-80 cells per donor;
- eight genes per replicate;
- 100 replicates per scenario;
- donor random intercepts;
- gene-specific donor effects and within-donor coupling;
- 50 exposure-defined neighborhoods.

## Scenarios

- null;
- donor-positive only;
- within-positive only;
- concordant positive;
- donor positive with within negative;
- donor negative with within positive.

## Results

Mean donor-versus-within reversal fraction:

| Scenario | Reversal fraction |
|---|---:|
| Null | 20.1% |
| Donor positive only | 27.8% |
| Within positive only | 32.0% |
| Concordant positive | 2.3% |
| Donor positive / within negative | 95.5% |
| Donor negative / within positive | 97.3% |

The two genuine opposing scenarios were therefore separated from concordant and
mostly null scenarios. The residual 20% classification rate under the null
scenario reflects the pre-specified effect-size threshold and finite sampling;
it is a screening-threshold calibration, not a biological false-positive rate.
The empirical panel was therefore interpreted with both expression-matched
sets and a conditional effect-pairing null. The latter preserves effect-size
marginals and breaks gene-wise unit pairing.

## Software

Package:

```text
work\unitaware_framework\
```

Core functions:

- donor pseudobulk;
- neighborhood-expression association;
- within-donor cell rank;
- standardized reversal scoring;
- composition-expression decomposition;
- simulation generator;
- conditional effect-pairing permutation.

Unit tests:

```text
python -m unittest discover -s tests -v
```

Current test status: six tests passed.

## Files

```text
outputs\unitaware_simulation\simulation_runs.csv
outputs\unitaware_simulation\simulation_summary.csv
outputs\figures\Figure_S2_unitaware_simulation.tif
outputs\figures\Figure_S2_unitaware_simulation.png
```
