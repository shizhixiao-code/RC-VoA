# RC-VoA

Reproducibility materials for **RC-VoA: Certified Allocation Auditing and Verification for Skin-Lesion Imaging**.

This release matches the **17 November 2026** manuscript revision (8-page main paper and 5-page supplement). It includes the latest **26,025-allocation quota experiment** and **complete-case endpoint analysis**. Supplementary tables are numbered SI–SVIII in order of appearance. The submission snapshot is tagged `submission-2026-11-17`. Paper identity and input archive hashes are recorded in [`docs/MANUSCRIPT_VERSION.json`](docs/MANUSCRIPT_VERSION.json).

中文说明：本仓库可用冻结预测/决策和终点数据重算论文中的分配审计、配额诊断与 complete-case 结果。无需原始图像或 GPU。快速检查主要核对已有结果；完整数值重算与原始模型重新训练是不同层次，见下方命令与范围说明。

## Start here

Use Python **3.13**. The tested numerical environment is NumPy 2.5.1, pandas 3.0.5 and SciPy 1.18.0. It uses SciPy's bundled HiGHS MILP solver; no commercial solver is required.

```bash
git clone https://github.com/shizhixiao-code/RC-VoA.git
cd RC-VoA
python -m venv .venv
# Windows PowerShell: .venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -r requirements.txt

# Verify release checksums, examples, archived results and current table mappings.
python run_reproduction.py --suite quick

# Recompute the main allocation, quota and endpoint results plus verification checks.
python run_reproduction.py --suite standard
```

The repository is **public**. Original images and pretrained weights are obtained separately under their providers' terms; see [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md).

Each command creates a new directory under `runs/`, copies only the needed modules there, and writes logs plus `REPLAY.json`. A successful run ends with `PASS` and the receipt path. Published inputs/results are not overwritten. Use `--work-dir PATH` to choose another run location. Do not run the legacy module commands directly inside the checked-in folders if you want to preserve their archived outputs.

## Reproduction suites

| Suite | What it actually does |
|---|---|
| `quick` | Checks allocation examples, summarizes 23,600 archived verification verdicts/costs, checks 150 policy rows and related source data, reconstructs source-inclusion/transport tables, and summarizes endpoint provenance. This is **not** a rerun of all query orders. |
| `core` | Recalculates complete-criterion allocations and ablations with MILP, checks small enumerated cases, and reruns target-sensitivity analyses. |
| `quota` | Replays all 26,025 admissible within-pool quotas for the three frozen rankings, then independently verifies TP/FP/flags/protected losses/verdicts using cumulative counts. |
| `endpoints` | Rebuilds 236 tables under three endpoint settings (1,524 table/setting/budget checks), reruns allocation MILPs and six lesion-level checks for the changed configuration, and generates current endpoint tables. |
| `verification-checks` | Runs hidden-completion/literal-allocation tests and the portable stopping-boundary audit. |
| `standard` | Runs `quick`, `core`, `quota`, `endpoints`, and `verification-checks`. It does **not** rerun the full timed verifier benchmark, bootstrap, or image training. |
| `verification-full` | Reruns all 20,000 primary/repaired table-order audits and the 560 sequential timing pairs, then summarizes them. Computationally heavier; timings depend on hardware. |
| `bootstrap` | Recomputes the original conditional paired group-bootstrap repair interval with frozen fits/decisions. |
| `source-inclusion` | Extracts the sealed expanded-source package and reruns allocations, conditional bootstrap and SCCIS/IEC equivalence under both updated-R and fixed-R contracts. It uses already fitted decisions; it does not retrain image networks. |

Examples for the two newest experiments:

```bash
python run_reproduction.py --suite quota
python run_reproduction.py --suite endpoints
```

## Main results to expect

| Analysis | Expected result and denominator |
|---|---|
| Original versus repaired allocation | 1/100 versus 14/100 configurations at reference targets; all 14 repaired passes are from UQ. |
| Policies at 20% imaging / 50% flag cap | Each policy passes 0/50 overall and 0/10 in the feasible subset. |
| Fixed-ranking quota diagnostic | 26,025 allocations; neither quota projection nor exhaustive quota choice makes any tested policy pass. |
| Explicit external AK coded noncancer | Repaired feasibility remains 14/100. All resource records are retained. |
| Complete-case | Repaired feasibility becomes 13/100; primary remains 1/100. Counts hold at both 15% and 20% imaging. |
| Expanded source AK/KA inclusion | Repaired feasibility is 7/100 under either updated or fixed source-S reference. |
| Verification ledger consistency | 23,600 archived table-order verdicts agree with full-record verdicts; 20,000 primary/repaired comparisons have the exact-completion comparator. |

**Complete-case means removal of all 199 reference endpoint exclusions**: 92 explicit AK and 107 other ambiguous/unspecified endpoints, leaving 8,969 external lesions. Retained labels, decisions and fold memberships remain fixed; resource-dependent quantities are recalculated. It is neither diagnostic adjudication nor a correction for missingness/selection bias.

The single changed configuration is UQ local-head batch 3 at the 50% flag cap. Removal reduces total downgrade support from 160 to 158, while the requirement decreases from 160 to 159. Four required rescues and 158 benign downgrades remain unchanged. The resulting deficit is one downgrade. This endpoint experiment did not rerun policy passage or verification query costs.

## Repository layout

```text
core/                   Frozen allocation/verification inputs, algorithms and archives
source_data/            Policy cells, witnesses, figure/table source data and checks
quota/                  Latest fixed-ranking quota replay and independent checker
endpoint_sensitivity/   Original endpoint rules, inputs, raw-field provenance and replay
archives/               Sealed expanded-source P1 numerical/training-receipt package
scripts/                Current-manuscript table summaries
docs/                   Results map, data/access boundaries and manuscript identity
validation/             Checks performed for this GitHub delivery
run_reproduction.py     Portable entry point; isolates generated outputs in runs/
MANIFEST.json           SHA-256 checksums for delivered scientific files
```

[`docs/RESULTS_MAP.md`](docs/RESULTS_MAP.md) is the authoritative mapping to **current** figures/tables. Some unchanged historical filenames and module READMEs refer to older table numbers or manuscript page counts. Their numeric inputs/results are preserved; those old labels are not the current manuscript index.

## Scope and access

This is a frozen-data numerical reproduction release. It includes coded lesion/group identifiers and endpoint labels needed for matching and grouped analyses, not only aggregate counts. No original lesion image pixels, neural-network weights, credentials, or patient contact records are included. The sealed P1 ZIP retains historical manuscript snapshots as provenance; these are not the current submission PDFs. The current manuscript is identified by hashes in `docs/MANUSCRIPT_VERSION.json`. See [`docs/DATA_ACCESS.md`](docs/DATA_ACCESS.md) for provider sources and limitations.

The verifier assumes previously established endpoint eligibility and correct revealed diagnoses. Read-count savings exclude eligibility ascertainment/curation; in-memory timings exclude real clinical retrieval time. Configurations reuse cases and fitted models; 100 configurations and 23,600 table-order evaluations are not independent patient counts.

Raw-image retraining requires provider-authorized data, preprocessing, pretrained weights and the original training environment. The sealed P1 archive retains its training recipes/receipts and `REPRODUCE.md`, but this repository alone is not a self-contained image-training dataset. No project-wide open-source license is newly assigned to third-party data or archived material by this upload.
