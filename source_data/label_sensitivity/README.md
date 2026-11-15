# Supplementary Table SVII: source-inclusion sensitivity

Run `python check_table.py` with standard-library Python to reconstruct `Table_SVII.csv` and check the printed LaTeX rows when manuscript sources are present.

The table selects the **reference external endpoint setting**, original clinical targets, free quotas, and 15%/20% imaging. Source S/P has 36 configurations; primary and repaired each have 100 (50 at each 30%/50% flag cap). Updated R uses refitted source-S decisions. Fixed R retains the original per-case reference and recalibrates repair on the disjoint calibration groups. The two contracts are reported separately.

The expanded source includes 99 AK as provider-negative and 162 KA as provider-positive in the SCC-family stratum, increasing 4,717 to 4,978 duplicate components. Source and low-workload model/calibration branches were rerun with the original recipes and folds; other exclusions remained. SCCIS remains an independent 90%-floor stratum. Original policy, verification and timing results are unchanged.

`pipeline_*` and `frozen_*` CSVs are unchanged outputs of the completed P1 experiment. The compact table does not combine these with alternative clinical-target experiments or the complete-case endpoint sensitivity.

For the full P1 numerical replay, extract `RCVoA_P1_Full_Results.zip` from this module in the reproducibility archive, preserving its two sibling directories, and follow its `REPRODUCE.md`. That sealed package contains the complete per-case decision inputs, allocation checks, conditional bootstrap and training receipts. Raw dataset files and neural weights remain separate. No new fitting is required to reproduce the table or its numerical audits.
