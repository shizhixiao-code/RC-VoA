# Mapping to the 15 November 2026 manuscript

Use this mapping instead of historical table numbers inside the archived modules.

| Current item | Evidence / command |
|---|---|
| Results B, 1/100 to 14/100 and constraint bottlenecks | `core/tables.json`, `core/results/`, `core/target_sensitivity/`; `--suite core` |
| Main Figure 3 / Table III, clinical-input policy realization | `source_data/policy_cells.csv`, `source_data/witness_cells.csv`, `source_data/source_data/`; `--suite quick` checks original 150 policy rows and witnesses |
| S4 quota sweep/projection, 26,025 allocations | `quota/inputs/`, `quota/recovered/prepared/`, `quota/results/`; `--suite quota` |
| Figure 4 / Table SV, endpoint verification | `core/verification_ledgers.csv`, `core/exact_completion/`, `source_data/exact_verification/`; `quick` summarizes archived costs, `verification-checks` tests boundaries, `verification-full` reruns query-order experiments |
| Table SVI-A, source-training inclusion | `core/label_sensitivity/`; quick check writes the legacy-named `Table_SVII.csv`. Full numerical replay is `--suite source-inclusion` |
| Table SVI-B, external endpoint alternatives | `endpoint_sensitivity/results/feasibility_summary.csv`; `--suite endpoints` reconstructs the allocation results, and outputs `endpoint_tables/Table_SVI_B.csv` |
| Table SVII / S5, threshold transport | `core/threshold_transport/`; quick check reconstructs legacy-named `Table_SV_current.csv` from the original transport records |
| Table SVIII, endpoint provenance/exclusions | `endpoint_sensitivity/evidence/raw_field_slices/`, `results/raw_diagnosis_inventory.csv`; quick/endpoints output `Table_SVIII_counts.csv` and `Table_SVIII_recorded_bases.csv` |
| Repair uncertainty | `core/repair_uncertainty/`; `--suite bootstrap` |

The provenance summary lists the four external cohorts. The development MILK row in Table SVIII is documented by the source-training protocol and manuscript: 4,717 histopathology-confirmed components, with no binary-endpoint exclusions in that selected development population.

Endpoint settings are separate from source-training changes. `reference` retains excluded records for imaging/flag resources. `AK_noncancer` changes only the binary side of 92 explicitly diagnosed AK records. `endpoint_complete_case` removes every original excluded row, including AK, and updates the resource population. Do not mix these settings when comparing counts.

Query costs reported in the paper remain the reference-setting experiment. A complete-case result does not update the policy success denominator or imply new query-cost estimates.
