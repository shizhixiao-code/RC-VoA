Current manuscript mapping: Figure 4 and Supplementary Table SVI. Original 2026-11-07 experiment and timings are retained.

# Exact-completion verification comparator

The frozen inputs, original 100 outcome-blind orders, results for all 20,000 table-order audits, 5,600 boundary pairs and 560 sequential timing pairs accompany the manuscript revision of 2026-11-07. The solver receives observed counts and public frozen decisions/eligibility only. The full outcome arrays are retained solely by the replay custodian and validation code.

Python 3.13.14, NumPy 2.5.1, pandas 3.0.5, SciPy 1.18.0 / HiGHS 1.12.0 were used for the reported experiment. Install `requirements.txt` and run from this directory:

```sh
python validate.py
python audit_boundaries_portable.py
python run_experiment.py pilot
python run_experiment.py full
python run_experiment.py online --orders 10
python summarize_portable.py
```

The portable boundary audit checks the packaged frozen-input hashes and reruns the independent allocation/completion models. The unchanged original `audit_boundaries.py` additionally checks the author's original workspace files; its signed provenance and original receipt are retained. `prepare.py` is the original provenance extraction script and is not required for a frozen-input replay. The portable summary differs only in the location of the included necessary-query bounds.

Run replays on a copy to preserve archived timings: new runtime values will depend on the machine. The offline bisection timings are not deployment timings. The 560 sequential pairs use the first ten archived orders, complete ledger/model/solver evaluation, and in-memory endpoints; no clinical retrieval latency is included. Default solver threading was retained. A timeout/error remains unresolved and triggers the next group; offline bracketing aborts if any prefix is unproved.

Reported checks use equal-table means after eligible-record normalization; P10--P90 is computed across the 100 matched order means. The 72 shared zero-query configurations per stage are excluded from retrieval-requiring strata. The complete benchmark is at most 20% imaging with free quota, not an exact-volume or fixed-quota comparator. Earliest stopping is conditional on each frozen order.
