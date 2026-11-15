# Table SIV and Figure S2 data archive

The compact printed Table SIV is expanded in Table_SIV_full.csv. All records are retained in the clinical workload denominator; excluded diagnoses are not cancer/benign endpoints. Percentages use explicit counts/denominators. S_local rows pool five independently calibrated held-out thresholds; oracle_cohort is a single evaluation-label threshold, and oracle_foldwise is a separate retrospective like-partition comparator. Neither oracle is deployable.

Thresholds retain float64 precision. For threshold reproduction from CSV in pandas use float_precision='round_trip'; the prepared NPZ files retain native exact arrays. Do not round thresholds before applying the >= whole-tie rule.

Each clinical backbone is shared by L1–L3 paired systems. Table/fold/budget configurations are dependent; do not treat them as independent patients. evaluation_fold = -1 denotes a pooled cohort table; 0–4 are the archived five folds. Primary audit comparisons are matched within held-out folds; pooled 36-table comparisons connect to the original manuscript.

Uncertainty is a 1,000-draw paired grouped percentile bootstrap, conditional on frozen fits and thresholds. Absent denominators are undefined; degenerate intervals are flagged. Calibration coverage does not certify a population sensitivity floor. No target was relaxed.

All six audit margin blocks and their subtype components are provided. Negative margins may co-occur. minimum_images = -1 means support makes the criterion impossible even without the imaging cap. Exact-volume feasibility and actual ranker passage are distinct from at-most-budget feasibility. The original source-S reference remains protected.

Old detailed clinical/certificate tables remain in the original result files and the before-threshold-transport manuscript snapshot. No prior experimental record has been removed. New evidence does not alter the original 14/100 repair or 23,600 verification results.
