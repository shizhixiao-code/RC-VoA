# Submission update — 17 November 2026

`submission-20261117/` records the PASS result of the six-step quick suite after updating the manuscript identity and supplementary-table exports. Scientific input files and expected numeric results are unchanged. A clean build of the local submission source reproduced all 13 final PDF pages pixel-for-pixel, with 64 reciprocal links. The full solver/training suites were not rerun for these presentation-only changes.

# Delivery validation — 15 November 2026

Validated on Python 3.13.14 / Windows, NumPy 2.5.1, pandas 3.0.5, SciPy 1.18.0. See the receipts and individual subprocess logs for measured durations and exact commands. Absolute delivery paths are replaced with `<REPOSITORY>` and `<RUN_DIR>`.

- `standard`: PASS, all 16 subprocesses. Includes the independent 26,025-allocation quota check, endpoint MILP replay and changed-case checks, main allocation replay, source-data checks, and exact-verifier validation/boundary audits.
- `source-inclusion`: PASS, all 5 subprocesses. Includes updated-R and fixed-R allocation replays, their conditional bootstraps, and SCCIS/IEC equivalence.
- Boundary audit: 5,600 opposite-completion pairs, 10,232 independent MILPs, zero disagreements.
- `PACKAGE_SCAN.json`: Python parse checks and credential-pattern scan, including the sealed ZIP contents. This is not a formal privacy or data-license audit.

The full timed verifier suite, the separate original-repair bootstrap, and image-level training were not rerun for this delivery. Their original artifacts/recipes remain included. Expanded-source conditional bootstraps were rerun as stated above.

The root manifest covers the delivered package, including nested original manifests. It excludes itself, generated runs, Git metadata and these delivery receipts. Initial suite receipts report the manifest entry count at their execution time; subsequent documentation-only changes and broader manifest coverage do not change their scientific inputs.
