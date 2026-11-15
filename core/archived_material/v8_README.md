# RC-VoA P2 completion release

The original scientific files are retained. Added modules:

- repair_uncertainty/: 1000 paired conditional group bootstraps at 20% imaging under reference targets; 100 matched configuration pairs, 1055 groups, observed gain 13 pp, percentile interval 5–19 pp. Includes configuration-level fold stability. Run `python repair_uncertainty/run_bootstrap.py` on a copy.
- threshold_transport/: restored existing S5 source data, fold thresholds/results, curves, 1728 uncertainty rows and prepared decisions. Run `python threshold_transport/check_transport_data.py`. It regenerates current Table SV using Delta FP = baseline minus final FP.

Tested scientific environment: Python 3.13.14, NumPy 2.5.1, SciPy 1.18.0, pandas 3.0.5; SciPy bundled HiGHS 1.12.0. Each module documents scope and commands. The paired repair interval conditions on frozen fits/decisions; it is not a new training experiment. P1 exclusion/label-mapping sensitivity remains deferred.

## Retained reference documentation

# RC-VoA reproducibility package

Reference experiment documentation from the preceding revision. Original policy and verification results are retained; `target_sensitivity/` adds the exploratory analyses and system/cohort summary.

## Tested environment and commands

Python 3.13.14, NumPy 2.5.1, pandas 3.0.5, SciPy 1.18.0, with SciPy's bundled HiGHS 1.12.0. Install `requirements-tested.txt` for the tested package versions. The original root `requirements.txt` records the older count-audit module's minimum dependencies.

```text
python -m pip install -r requirements-tested.txt
python reproduce.py --example
python reproduce.py --all --check
python replay_verification.py
python target_sensitivity/analyze_targets.py
python check_target_data.py
```

The exact-completion module has its own commands in `exact_completion/README.md`, including `validate.py`, `audit_boundaries_portable.py`, and sequential replay. Run experiments on a copy when preserving archived output/timing files.

## Contents and manuscript mapping

| Entry | Content |
|---|---|
| `tables.json` | 236 frozen cell-count tables: 36 source-S/P, 100 primary, 100 repaired |
| `reproduce.py` | Complete allocation audits, minimum-image witnesses and cumulative/single-removal ablations |
| `results/` | Original reference outputs, including all 5,664 ablation judgments |
| `verification_ledgers.csv` | 23,600 archived table-order audit judgments/costs |
| `exact_completion/inputs/` | Frozen query inputs and the 100 outcome-blind group orders for 100 primary plus 100 repaired tables |
| `exact_completion/` | 20,000 matched verification results, 5,600 stopping-boundary pairs and 560 sequential timing pairs |
| `target_sensitivity/` | 1,652 new table/profile judgments, 309 allocation witnesses, original-reference reconciliation and all system/cohort/cap counts |
| `archived_material/` | Previous supplementary Figure S1, previous table notes and retained model-selection details |

Target sensitivity uses at-most-K allocations and free quotas; its original complete-criterion row is the reference. The grid was added after review and changes one target at a time, retaining frozen decisions, labels and both baselines' detected cancer identities. It does not re-evaluate policy success under replacement targets.

`tables.json` cells are `[diagnosis, B, Dpair, R, count]`; Excluded cases use imaging/flag resources but not binary endpoint counts. Counts are sufficient for exact allocation under these constraints. Group-query structure for verification is supplied separately in `exact_completion/inputs/`.

Timing corresponds to in-memory sequential replay, including ledger/model/solver work. Repaired-feasible means are 0.0236 s for RC-VoA and 1.6292 s for Exact-Completion MILP. Original measurements used the CPU and software recorded in the exact-completion module; new runtimes can vary.

The 14 repaired feasible tables consist of 4/50 at a 30% flag cap and 10/50 at a 50% cap, across both cohorts. They are descriptive counts from a grid sharing cases, model fits and calibration. All 150 realized-policy rows for the 50%-cap subset are in the companion Source Data archive.

## Consistency revision

The 50%-cap selector experiment uses repaired tables under reference targets. Its ten feasible configurations are the subset of the original 14/100, separate from 26/100 under the relaxed flag target. The mean 5.47 saved records uses all 1,400 repaired-feasible table-order runs.

Run `python target_sensitivity/summarize_reference_batches.py` to reproduce S4: reference-target repaired feasibility across batches 1–5 is 7, 0, 4, 3, 0 out of 20 configurations per batch. `check_target_data.py` validates this added summary as well as the original target/system counts. Original scientific inputs, experimental result rows and timing ledgers are preserved.
