# RC-VoA selector Source Data

This companion to RCVoA_Reproducibility.zip reproduces the frozen policy comparisons in Results III-C, Figure 3 and Table III. Unzip the archive and run these commands in its top-level directory.

    python -m pip install -r requirements.txt
    python check_source_data.py

The checker derives final decisions from frozen clinical decisions and the acquisition indicator. It recomputes counts, the complete criterion, feasible negative-pool quotas at the full imaging volume, and the first failed audit level. It reconstructs the three original minimum-image witnesses used in Table III. A successful run prints a JSON report and exits with code zero.

- **policy_cells.csv** groups the frozen lesion-level acquisition records by configuration, policy, endpoint category, clinical decision B, paired decision P, reference R and acquisition indicator. Each row gives the number of lesions with those values. Batches and policies share cases.
- **witness_cells.csv** contains selected counts from the archived minimum-image oracle witnesses. The other cells retain their clinical decisions. These witnesses used 164, 168 and 160 images. The matched policies used 310, 311 and 307 images.
- **source_data/** contains all 150 published comparison rows, the ten Fig. 3 rows, both Table III tables, the three-policy summary and the previously archived displays. The original case illustration credits the UQ provider in its caption. Fig3_quota_selection.pdf is the published diagnostic figure.
- **frozen_audit.py** reuses the full-criterion and margin functions from the companion reproducibility archive. A subtype with zero observed cancers supplies no sensitivity evidence.
- **SHA256SUMS.txt** lists SHA-256 digests of the distributed files.

This is an aggregate replay of historically frozen decisions. The training and case-ranking stages precede these inputs. The manuscript describes grouped evaluation, calibration and prior cohort inspection. The separate reproducibility archive supplies the constraint ablation, numerical witness optimization and endpoint-ledger replay.


## Verification comparator revision

Figure 4 and Supplementary Table SV are included in `exact_verification/`. `python check_source_data.py` now checks the original Fig. 3/Table III displays and the added matched verification summaries.
