Current manuscript mapping: Figure 4 and Supplementary Table SVI. Original 2026-11-07 experiment and timings are retained.

# Figure 4 and Supplementary Table SVI

`full_results.csv`: 20,000 matched table-order pairs, including the common zero-query screen. `summary_by_order.csv`: equal-table order means for the four retrieval-requiring strata. `summary.csv`: mean/P10/P90 and separate online timing summaries. `online_results.csv`: 560 sequential replay pairs (three methods, 1,680 rows). `fig4_cumulative_curves.csv`: all unsampled panel-a curve points. `positive_1400_bounds.csv`: archived necessary-query bound for repaired feasible tables. The lower bound and the exact completion certificate are distinct objects.

All judgments match full records. Infeasible retrieval-requiring audits stop at the same prefix for both methods (4,100/4,100). Repaired feasible audits use 99.984% versus 99.627% of eligible endpoints on average, with 419 earlier MILP stops. P10--P90 describes query-order variation, not population confidence intervals. Online timing includes ledger updates/model/solver work and excludes clinical retrieval latency.
