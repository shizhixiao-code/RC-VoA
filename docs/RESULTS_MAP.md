# Mapping to the 17 November 2026 manuscript

Use this mapping instead of historical table numbers inside the archived modules. Submission snapshot: `submission-2026-11-17`.

Current reading order is SI–SIV (supplement page 2), SV (page 3), SVI (page 4), and SVII then SVIII (page 5). The 17 November update changes presentation and indexing, not scientific inputs or reported results.

| Current item | Evidence / command |
|---|---|
| Results B, 1/100 to 14/100 and constraint bottlenecks | `core/tables.json`, `core/results/`, `core/target_sensitivity/`; `--suite core` |
| Main Figure 3 / Table III, clinical-input policy realization | `source_data/policy_cells.csv`, `source_data/witness_cells.csv`, `source_data/source_data/`; `--suite quick` checks original 150 policy rows and witnesses |
| S4 quota sweep/projection, 26,025 allocations | `quota/inputs/`, `quota/recovered/prepared/`, `quota/results/`; `--suite quota` |
| Figure 4 / Table SVI, endpoint verification | `core/verification_ledgers.csv`, `core/exact_completion/`, `source_data/exact_verification/`; `quick` summarizes archived costs, `verification-checks` tests boundaries, `verification-full` reruns query-order experiments |
| Table SVIII-A, source-training inclusion | `core/label_sensitivity/`; quick check writes the legacy-named `Table_SVII.csv`. Full numerical replay is `--suite source-inclusion` |
| Table SVIII-B, external endpoint alternatives | `endpoint_sensitivity/results/feasibility_summary.csv`; `--suite endpoints` reconstructs the allocation results, and outputs `endpoint_tables/Table_SVIII_B.csv` |
| Table SVII / S5, threshold transport | `core/threshold_transport/`; quick check reconstructs legacy-named `Table_SV_current.csv` from the original transport records |
| Table SV, endpoint provenance/exclusions | `endpoint_sensitivity/evidence/raw_field_slices/`, `results/raw_diagnosis_inventory.csv`; quick/endpoints output `Table_SV_counts.csv` and `Table_SV_recorded_bases.csv` |
| Repair uncertainty | `core/repair_uncertainty/`; `--suite bootstrap` |

The provenance summary lists the four external cohorts. The development MILK row in Table SV is documented by the source-training protocol and manuscript: 4,717 histopathology-confirmed components, with no binary-endpoint exclusions in that selected development population.

Endpoint settings are separate from source-training changes. `reference` retains excluded records for imaging/flag resources. `AK_noncancer` changes only the binary side of 92 explicitly diagnosed AK records. `endpoint_complete_case` removes every original excluded row, including AK, and updates the resource population. Do not mix these settings when comparing counts.

Query costs reported in the paper remain the reference-setting experiment. A complete-case result does not update the policy success denominator or imply new query-cost estimates.
