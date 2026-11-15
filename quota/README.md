# Fixed-ranking quota diagnostic

Run `python run.py --replay` and then `python check_prefix.py` on a working copy.
The root runner creates that copy automatically. NumPy and pandas are required.

The 50 repaired configurations use a 20% image budget and 50% flag cap.
All three frozen within-pool rankings and the image volume are retained.
The exhaustive quota sweep contains 26,025 allocations; projection onto the
oracle-feasible quota interval is also evaluated. Each policy remains 0/50
overall and 0/10 in the feasible subset. `check_prefix.py` independently checks
TP, FP, flags, protected losses and verdicts using cumulative counts.

This is a post-hoc diagnostic, not a deployable selector or new model training.
The sole packaging change to run.py replaces the unused original-workstation
recovery branch with packaged relative paths. Inputs and evaluation logic are unchanged.
