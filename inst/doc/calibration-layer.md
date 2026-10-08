# External-null empirical calibration

La Terra T2B separates the scientific null generator from the generic
calibration engine. An `lt_branch_state_ensemble` records externally supplied
branch-by-replicate values. Every replicate column is one coherent branch-state
world applied to every gene. La Terra does not shuffle branch states, simulate
traits, choose a phylogenetic model, or certify the scientific null hypothesis.

`lt_calibrate_association()` freezes the observed T2A gene-specific domain
before reading null scores. Branches are reordered only through exact
scientific keys. A missing required null value produces
`NULL_DOMAIN_MISMATCH` for the affected gene. A binary replicate with one
group, or a continuous replicate with no predictor variation, produces
`NULL_CALIBRATION_INCOMPLETE`; no smaller replicate count or empirical p-value
is released silently.

For `R` complete replicates, greater, less, and two-sided tails use inclusive
ties and the plus-one rule. The minimum possible p-value is `1/(R+1)`. The
reported Monte-Carlo standard error is descriptive computational precision,
not a biological confidence interval.

The operator streams scores one null world at a time and does not require a
gene-by-replicate score matrix. For an `lt_log_relative`, the identical ordered
replicate ledger is consumed separately by ZERO_MASS and POSITIVE_LOG.
ZERO_MASS x Null A/B and POSITIVE_LOG x Null A/B are distinct hypotheses with
inclusive ties, plus-one empirical p-values, and fail-closed null
non-estimability. Components, Null A/B, and their p-values are never pooled,
averaged, minimized, or converted into a joint p-value.

Multiple-testing adjustment is a separate explicit step through
`lt_adjust_calibration()`. The default is none. One family is exactly one
trait by one rate representation by one domain mode by one statistic by one
alternative. BH is exposed without a universal arbitrary-dependence claim;
BY is the conservative arbitrary-dependence FDR option. Holm and Bonferroni
are family-wise error procedures.

Boundary BH and BY adjustments are separate families for each component and
null hypothesis. Calibration results certify only the computation. They do not certify the
external generator, establish causality, validate prediction, identify
adaptive genes, or select a preferred rate representation.

The generic `lt_calibrate_association()` route rejects historical finite
`lt_rate(representation = "logGBI")` inputs before any null score is computed.
Current boundary calibration is component-separated. Historical global-k10
calibration is not an ordinary La Terra V1 operation:
`lt_calibrate_logGBI_encoded()` and its adjustment route always fail closed.
Frozen evidence may be compared only through the statistic-free,
source-bound `validate_legacy_logGBI_replay()` route. Use GBI for V1 ordinary
inference.
