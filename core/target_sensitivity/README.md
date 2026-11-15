# Frozen-table target sensitivity and system summaries

This module reproduces Supplementary Tables SII and SIVB in the 2026-11-08 target-sensitivity revision. `tables.json` contains 236 original frozen count tables: 36 source-S/P, 100 primary and 100 repaired. No patient identifiers or images are included.

Run with Python 3.13.14 and the tested dependencies:

```text
python -m pip install -r requirements-tested.txt
python analyze_targets.py
```

The tested environment was Python 3.13.14, NumPy 2.5.1, SciPy 1.18.0 and SciPy's bundled HiGHS 1.12.0. Install SciPy rather than an additional standalone solver. New elapsed times and equally optimal witnesses can vary across machines; verdicts, minima and summary counts must agree.

`PROTOCOL.json` records the exploratory grid written before these new profiles were run, after review of the original experiment. Every profile changes one requirement from the original complete criterion: pooled sensitivity 90 to 85%; each subtype floor down five percentage points; FP reduction 10 to 5%; flag reduction ten to five percentage points; absolute caps up ten percentage points; or imaging 20 to 25%. Both baselines' detected cancer identities remain protected. Source S/P retains its 100% absolute cap. Decisions, labels, table membership and fitting/calibration remain frozen. These are at-most-K, free-quota audits.

Each table/profile is solved twice: by generalized closed-form margins and an independently specified allocation MILP. Every feasible MILP witness is rechecked with exact integer inequalities. The script also checks 294 small-table/profile combinations by literal subset enumeration. A timeout or solver error aborts the run instead of producing an infeasible result.

## Output files

| File | Contents |
|---|---|
| `results/target_sensitivity_all.csv` | All 1,652 table/profile judgments, requirements, margins and minimum image counts |
| `results/target_summary.csv` | Exact family-level counts for the seven profiles |
| `results/system_cohort_summary.csv` | System/cohort/cap counts and feasible batch identifiers under the reference criterion |
| `results/denominator_reconciliation.csv` | Cohort and cap decomposition of the two 100-table grids |
| `results/allocation_witnesses.json` | Selected counts for all 309 feasible allocations |
| `results/AUDIT.json` | Validation counts, tested versions and SHA-256 provenance |

System IDs map to manuscript names as follows: `headonly` = Source heads; `laststage` = Source features; `multiscale` = Multiscale; `localheads` = Local heads; `localfeatures` = Local features. Source-S/P systems are the three backbones crossed with L1--L3. Table SIIA groups the three fusions per backbone; the CSV preserves all nine systems separately.

The 100 low-workload configurations are two cohorts by five variants by five batches by two original flag caps. Under reference targets, 4 of 50 tables at the 30% cap and 10 of 50 at the 50% cap are feasible: (4 + 10)/(50 + 50) = 14/100 across both cohorts. Counts are descriptive properties of this dependent grid. Existing policy outcomes and the original complete criterion remain the reference results.

## Reference-target batch summary (consistency revision)

Run `python summarize_reference_batches.py` to reproduce the added S4 sentence and `results/reference_batch_summary.csv`. Across evaluation batches 1 through 5, repaired feasible counts are 7, 0, 4, 3 and 0 out of 20 configurations per batch (two cohorts, five variants and two original flag caps), under reference targets. The resulting 0–35% range describes this fixed grid.
