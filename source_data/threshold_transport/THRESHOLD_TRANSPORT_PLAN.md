# Threshold transport robustness analysis

Analysis specification fixed before computing new oracle/local evaluation results.
This is a retrospective addition to an existing study, not a prospective registration.

## Inputs and estimands

Reuse the three archived MILK-fitted clinical backbones and their L1/L2/L3 paired heads. The nine systems contain three distinct clinical score vectors per cohort. Retain every original external record, binary outcome, fine subtype, group, exclusion, and frozen source threshold. All excluded records remain in workload denominators but not cancer/benign denominators. No CNN or fusion head is fitted. Source S uses MILK MEL/SCCKA/BCC/AKIEC and pooled floors; source P uses only pooled 90%. Reproduce every archived threshold from the matching MILK score file, including its hash and whole-tie maximality.

External target floors are pooled 90%, melanoma/invasive SCC/SCC NOS 95%, BCC/SCCIS (IEC) 90%. Integer requirements are ceil(rate*N/100); zero-denominator targets are not evaluated. Decisions use score >= threshold, never split ties. Cohorts are MRA, HIBA, UQ, MCR-SL; all 12 backbone/cohort combinations and their mapping to all 36 paired systems are retained.

## Existing folds and separation

Recover the stored five-fold assignment from all 36 archived adaptive-action files and verify that copies agree. Reproduce the original label-blind greedy size balancing, seed 20261023 + cohort index in the order above. In rotation f, evaluation is f, calibration is (f+1) mod 5, and the remaining three folds are unused design data for this no-refitting experiment. Verify group disjointness and case identity. Use the archived grouping unit; do not upgrade image-linkage groups to confirmed participants.

## Local threshold selection (prepare stage)

Give the selector only calibration scores, diagnoses and fine subtypes. For each represented target, take the required descending cancer-score order statistic; the minimum of these cutoffs is the highest whole-tie threshold satisfying every represented target. Freeze these 60 clinical thresholds and the 12 cross-fitted clinical decision vectors before evaluation reporting or oracle analysis. There is no search over folds, seeds, caps, models, objectives or alternate thresholds.

List all five subtype denominators for every calibration fold. A missing subtype has no calibration evidence and its floor cannot be certified. If a represented subtype has fewer than 20 cases for a 95% floor, or fewer than 10 for a 90% floor, flag that integer rounding requires zero calibration misses; do not call this reliable population calibration. Report achieved counts and held-out outcomes even for incomplete/sparse folds; never weaken a target or omit a difficult fold. If a calibration fold contains no cancer, retain source S as an explicitly uncalibrated fallback. Full evaluation target status is determined only after freezing. Source S remains the protection reference throughout.

## Oracle and operating-point evaluation (measure stage)

Enumerate all distinct clinical scores with >= whole-tie decisions, including the all-negative endpoint. For each full cohort/backbone, select the highest threshold meeting every observed pooled/subtype floor. This is an empirical retrospective oracle within the class of one scalar clinical threshold, not a deployable setting or a bound over all classifiers. Also archive foldwise evaluation-label oracle thresholds as a like-partition benchmark. Do not use either oracle to change S_local.

Report source S, cohort-wide oracle and cross-fitted S_local: threshold(s), flags/N, FP/benign N, pooled TP/cancer N, each subtype TP/N, floor attainment, and excluded counts. S_local aggregates held-out decisions from five separately frozen thresholds and need not lie on a single-threshold ROC/frontier. Report each fold separately in the data archive. Descriptive AUROC and class-score distributions may contextualize discrimination; these comparisons do not uniquely identify a causal distribution-shift mechanism.

Uncertainty: 1,000 paired bootstrap draws of the archived groups, stratified by evaluation fold, seed 20261102 + cohort index. Source and local decisions use identical draws. Freeze thresholds and models during resampling, matching the study's conditional policy uncertainty; intervals do not include threshold-fitting uncertainty. Archive pointwise percentile 95% intervals for flags, FP per 1,000 lesions, pooled/subtype sensitivity and paired differences. For each metric report contributing cases/groups, valid resamples, and degenerate intervals. Missing-denominator draws are undefined and counted, not imputed. Small or homogeneous subtype counts cannot establish a population guarantee.

## Complete-criterion audit (audit stage)

Compare B=source S and B=S_local, with D_pair=original source-P paired decisions and R=original source S. Preserve both current- and reference-detected cancer identities. Use the unchanged source-S/P criterion: absolute flag cap 100%, >=ceil(0.1*N) fewer flags than each baseline, <=floor(0.9*FP) against each baseline, the same pooled/subtype floors, and imaging budgets 10%, 15%, 20% with K=floor(b*N).

Primary unit: each held-out fold (180 tables per setting, nine systems x four cohorts x five folds), matched against source S on exactly the same records. Also report 36 pooled cross-fitted tables per setting to connect to the existing full-cohort 36-table experiment; these are aggregate audit objects, not a single deployable scalar threshold or independent replications. Report both at-most-K feasibility and exact-volume feasibility. Keep all six margin blocks, per-subtype margins, minimum K, valid quotas, oracle witnesses, and actual prespecified Proximity/Clinical-risk outcomes. Proximity orders abs(clinical score - the independently calibrated threshold for that case); Clinical-risk orders decreasing clinical score; case ID breaks ties. Rankings never read paired scores or outcomes.

Re-derive applicability for arbitrary frozen B and R before interpreting outcomes. Compare margin decisions/minimum acquisition counts and exact-volume feasibility with an independent lesion-variable MILP that constrains final decisions directly. Validate selected policy quotas where support holds, and verify every returned witness by direct clinical counting. Keep complete raw and summary outputs and preserve all earlier evidence files byte-for-byte.

## Figure and reporting contract

Figure S2: quantitative grid, one panel per external cohort; three backbone-specific sensitivity/workload curves, source/oracle/local markers, pooled 90% reference, and the region of each curve satisfying every subtype floor. Display each cohort denominator. Distinguish the cross-fitted local point from a single-threshold curve. Python/matplotlib only, editable PDF/SVG, source CSV, >=5 pt rendered text, panel-alignment gate and rendered collision/visual checks. No favorable-cohort selection. The conclusion will follow the data, including any persistent high workload or loss of sensitivity.

Table number and final page allocation follow user clarification. Main-paper additions must remain brief; the original source-S/P experiment, repair evidence, proofs, and verification remain intact.
