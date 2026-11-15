# RC-VoA review reproducibility package

This standalone archive accompanies the manuscript. It provides an executable
audit example and the complete post hoc constraint analysis without images,
patient identifiers, CNN weights or private workspace paths. It is an accompanying
review artifact; no public repository or DOI is claimed.

## Run

Use Python 3.10 or later in this directory:

```text
python -m pip install -r requirements.txt
python reproduce.py --example
python reproduce.py --all --check
python replay_verification.py
```

The example reconstructs two endpoint-informed minimum-volume witnesses:

| Frozen table | Minimum images | Required rescues | Required net FP reduction |
|---|---:|---:|---:|
| Primary UQ multiscale batch 3, 50% cap | 168 | 6 | 75 |
| Repaired UQ source heads batch 1, 50% cap | 164 | 4 | 76 |

The optimizer can choose a different equally optimal witness on another solver
version. Feasibility and minimum image counts should agree; individual identities
and all secondary witness metrics need not be unique. The full source archive
retains the original case-specific witnesses, fitting/calibration provenance,
folds, source probabilities, thresholds, endpoint queries and figure sources.

## Data dictionary

`tables.json` contains 236 tables: 36 source-S/P, 100 primary low-workload, and
100 repaired configurations. Each entry gives `table_id`, `family`, `cohort`,
`batch`, the original absolute flag percentage `cap`, lesion count `N`, and
`cells`. A cell is `[diagnosis, B, Dpair, R, count]`. Decisions are fixed binary
flags. Diagnosis is Melanoma, Invasive SCC, SCC NOS, BCC, SCCIS, Other cancer,
Benign, or Excluded. SCCIS implements the manuscript's IEC mapping. Excluded
cases contribute to image/flag budgets, but not sensitivity or FP counts.

`expected_feasible`, `expected_minimum`, and `expected_margins` are frozen
reference results used to check the replay; they are not optimization inputs.
Cell aggregation is exact for these allocation constraints because all cases
in a cell have identical coefficients. The counts do not support group-bootstrap
uncertainty or pathology-query order reconstruction.

The populations and model predictions are inherited from the source experiment;
multiple configurations share cases, fitting and calibration. Fourteen repaired
configurations are not fourteen independent clinical cohorts. The source-S/P
27/36 feasible tables require no additional cancer detections. Repaired positive
configurations require nonzero rescues and FP reductions. Neither result implies
that a clinical-only ranking will select a feasible allocation.

## Post hoc constraint analysis

`--all` uses at most 20% imaging with a free quota at every stage. It solves a
cumulative sensitivity-first order, a protection-first order, and removal of
each individual constraint from the complete set. Final full-criterion counts
are 27/36, 1/100 and 14/100. The source-S/P absolute cap is 100%; low-workload
tables retain their original 30%/50% caps. Source thresholds and decisions never
change. Earlier stages deliberately relax requirements and are not clinically
admissible alternatives to the complete criterion.

The eight constraint blocks are imaging, pooled sensitivity, subtype sensitivity,
protection of each current cancer detection, protection of each source-S cancer
detection, 10% FP reduction against both baselines, 10 percentage-point flag
reduction against both baselines, and the absolute flag cap. The first stage
permits zero acquisitions. Cumulative effects depend on order; single removals
can remain ineffective when another constraint independently rules out a table.
The original nonexclusive failure counts are retained and reproduced.

`results/all_configurations.csv` records all 5,664 table-stage judgments;
`summary.csv` gives counts; `cell_witnesses.json` stores cumulative witnesses;
`complete_margins.json` gives all final margins. `--check` validates 448 generated
small-table problems by exhaustive subset enumeration. `results/AUDIT.json`
records dependency versions and hashes. Original certificate-to-action checks
remain in the full archive under `analysis/frozen_landscape_action/`.

## Verification replay and scope

`verification_ledgers.csv` contains table-order costs and final judgments from
the existing source-S/P, primary, and repaired verification runs: 23,600 rows.
`replay_verification.py` recomputes judgment agreement and equal-table/order mean
checked fractions. This is a replay of archived ledgers, not a rerun of endpoint
queries. To rerun the original query algorithm, use `complete_verification.py`,
`verify_repaired_tables.py` and their frozen inputs in the full source archive.

Savings are checks avoided relative to retrieving all eligible records. They
are first normalized within each table, then averaged equally across tables and
the 100 query orders. The 74.16% and 67.19% values concern only queried infeasible
tables; 92.02% across all primary tables includes zero-query pre-screening.
Repaired feasible tables save about 0.016% and check about 99.984%.

## Analysis provenance and access

The constraint analysis, pool-specific calibration repair, local risk-protected
selector extensions and threshold-transport analyses are post hoc explorations.
Repair threshold selection uses calibration labels; evaluation decisions remain
frozen. Prior outcome inspection informed the choice of repair rule. The new
constraint analysis uses already frozen endpoint tables to diagnose feasibility;
it trains no model and evaluates no new prospective selector.

The full source archive includes provider citations and access instructions in
the manuscript references and existing data records. Original clinical and
dermoscopic images are not included in this replay package; access remains under
the original providers' terms. `PREPARATION.json` and `MANIFEST.json` record
provenance and hashes. No independent clinical efficacy claim is made by this
computational replay.


## Exact-completion addition (2026-11-07)

See `exact_completion/README.md` for the full matched Figure 4 comparison, integer completion models, boundary validation and sequential timings. Existing allocation examples and audits are unchanged.
