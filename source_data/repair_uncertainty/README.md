# Paired repair-gain uncertainty and configuration-level fold stability

This exploratory analysis supplements main III-C and supplementary S4. Its protocol was written before the new resampling run, after the reference results had been inspected. It does not change model fits, calibration thresholds, recorded labels, target rates, cohort membership or the original reference judgments.

## Reproduce

Tested with Python 3.13.14, NumPy 2.5.1 and SciPy 1.18.0 (bundled HiGHS 1.12.0):

```text
python run_bootstrap.py
```

Run on a copy to retain the original output files and timing. Verdicts, weights and percentiles are deterministic at seed 20261109; runtime and equally optimal MILP witnesses can vary. `prepare_inputs.py` documents extraction from the preserved pre-revision package and is not needed for portable replay using the included inputs.

## Design and estimand

The reference grid comprises 100 matched primary/repaired pairs: 2 cohorts × 5 fitted variants × 5 evaluation batches × 2 flag caps. At 20% imaging, the point estimate is the equal-configuration feasibility-rate difference, 14/100 minus 1/100 = 13 percentage points.

For each of 1,000 draws, sample G groups with replacement within each cohort/evaluation-batch stratum. The 10 strata contain 1,055 groups (602 MRA groups and 453 UQ participants). Keep all lesions from each sampled group, including excluded endpoints. Use identical multiplicities across the five variants, two caps and primary/repaired stages in that stratum. This preserves pairing and shared-case dependence across configurations; configurations are not treated as independent patients.

Recompute count margins, subtype/pooled detection requirements, both baselines' protected identities and relative workload requirements, imaging capacity and absolute flag capacity using the resampled N. Target percentages and numerical decision thresholds stay fixed. The 95% interval is the 2.5th–97.5th percentile of paired bootstrap differences, conditional on the frozen fits/decisions and the original evaluation-fold partition. It is not a refitting/recalibration or prospective-deployment interval. No hypothesis tests or p-values are used.

The single primary contrast is the aggregate 100-pair repair gain. Cohort/cap summaries are descriptive secondary outputs. `configuration_fold_stability.csv` gives primary/repaired pass counts and the five batch verdicts for all 20 fixed cohort×variant×cap combinations, addressing stability without pooling different configurations into an independent sample.

## Results and checks

The observed gain is 13 pp; its percentile interval is 5–19 pp. The bootstrap mean is 11.094 pp and is distinct from the observed-grid estimate. All 200 original count tables/verdicts were reconciled. The analysis computes 200,000 bootstrap table verdicts; 400 deterministically selected resampled tables agree with independently specified allocation MILPs and exact integer witness checks.

`results/` includes all 1,000 aggregate draws, all configuration verdicts, shared group weights, configuration-level fold counts, group counts per stratum, paired summaries and validation receipts. The current manuscript reports the overall 20%-imaging contrast; it does not compare 14/100 with its nested 10/50 subset as independent samples.
