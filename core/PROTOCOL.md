# Retrospective constraint ablation — 4 November 2026

This analysis was specified in response to manuscript feedback, after inspection
of the original results. It is a post hoc diagnostic analysis, not a prospectively
registered experiment and not an acquisition-policy evaluation.

Freeze the existing 36 source-S/P, 100 primary low-workload and 100 repaired
decision tables. Do not retrain a model, change a cutoff, alter a group split,
or select configurations by their outcome. The main comparison uses the original
20% imaging budget; all original 30%/50% absolute flag caps are retained when
that constraint is active. The source-S/P absolute cap is 100%.

For each table, an integer variable counts acquisitions in each cell defined by
the diagnostic category and the fixed (B, Dpair, R) decisions. Final positive
counts are n*B + x*(Dpair-B). Cell variables range from zero to cell size, with
no implicit cancer protection when a protection constraint is switched off.
This aggregation is exact because cases in a cell have identical coefficients
and no constraints on group-level acquisitions are imposed.

Primary cumulative order: imaging cap; pooled cancer sensitivity; subtype
sensitivity; protection of cancers detected by B; protection of cancers detected
by R; FP reduction against both baselines; relative flag reduction against both
baselines; original absolute flag cap. All stages solve for existence and minimize
acquisitions, with free quota. They do not require spending the complete budget.

Sensitivity analyses: an alternative order adds the two protection constraints
before sensitivity constraints; leave-one-out analyses remove each constraint
from the complete set. Retain the original nonexclusive failure counts; neither
cumulative losses nor leave-one-out changes are additive causal contributions.

Validation: check returned integer witnesses directly; compare all final complete
verdicts with frozen published records and the complete-criterion margin formula;
check cumulative monotonicity; compare solver answers with exhaustive subsets
for small generated tables. Archive every configuration-stage result and selected
cell count. Release an image-free, patient-ID-free count-table replay example.
