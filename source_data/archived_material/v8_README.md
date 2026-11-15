# RC-VoA P2 completion release

The original scientific files are retained. Added modules:

- repair_uncertainty/: 1000 paired conditional group bootstraps at 20% imaging under reference targets; 100 matched configuration pairs, 1055 groups, observed gain 13 pp, percentile interval 5–19 pp. Includes configuration-level fold stability. Run `python repair_uncertainty/run_bootstrap.py` on a copy.
- threshold_transport/: restored existing S5 source data, fold thresholds/results, curves, 1728 uncertainty rows and prepared decisions. Run `python threshold_transport/check_transport_data.py`. It regenerates current Table SV using Delta FP = baseline minus final FP.

Tested scientific environment: Python 3.13.14, NumPy 2.5.1, SciPy 1.18.0, pandas 3.0.5; SciPy bundled HiGHS 1.12.0. Each module documents scope and commands. The paired repair interval conditions on frozen fits/decisions; it is not a new training experiment. P1 exclusion/label-mapping sensitivity remains deferred.

## Retained reference documentation

# RC-VoA Source Data

Updated for the 2026-11-08 target sensitivity revision. Python 3.13.14 with NumPy 2.5.1, pandas 3.0.5 and SciPy 1.18.0 was used; SciPy supplies HiGHS 1.12.0.

```text
python -m pip install -r requirements-tested.txt
python check_source_data.py
```

The checker recalculates the original 150 policy rows, quota/selection diagnoses and three minimum-volume witnesses, verifies the original Figure 4 data, then validates the new target/system summaries and provenance hashes. To rerun all new allocation MILPs, use `python target_sensitivity/analyze_targets.py` on a copy.

| Directory/file | Content |
|---|---|
| `source_data/` | Frozen outcomes for three policies across 50 repaired configurations; the ten oracle-feasible rows, main Figure 3 and Table III data |
| `policy_cells.csv`, `witness_cells.csv` | Cell-count inputs for independent policy/witness reconstruction |
| `exact_verification/` | Main Figure 4 and Supplementary Table SVI data; 20,000 matched table-order records and sequential timing ledgers |
| `target_sensitivity/` | Supplementary Table SII and SIVB data; 236 tables, seven target profiles, complete results, MILP witnesses and executable analysis |
| `archived_material/` | Former supplementary Figure S1 and retained detailed table/model-selection notes |

All original outcome and timing data are preserved. The original 100-table repaired grid combines 30% and 50% flag caps; the 50-table policy comparison is its 50%-cap subset. Original clinical-performance intervals are pointwise and conditional on fitted policies; configuration counts and query-order percentiles have their separately stated descriptive meanings.

## Consistency revision

The 50%-cap selector experiment uses repaired tables under reference targets. Its ten feasible configurations are the subset of the original 14/100, separate from 26/100 under the relaxed flag target. The mean 5.47 saved records uses all 1,400 repaired-feasible table-order runs.

Run `python target_sensitivity/summarize_reference_batches.py` to reproduce S4: reference-target repaired feasibility across batches 1–5 is 7, 0, 4, 3, 0 out of 20 configurations per batch. `check_target_data.py` validates this added summary as well as the original target/system counts. Original scientific inputs, experimental result rows and timing ledgers are preserved.
