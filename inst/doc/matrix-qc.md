# Generic matrix-QC contract

La Terra standardizes the diagnostic procedure, not a universal scientific
threshold. The QC1 pipeline is strictly:

```
raw lt_matrix
  -> lt_diagnostic (diagnose only)
  -> explicit application decision
  -> Layer-2 measurement-fit eligibility
```

`lt_qc_matrix()` never alters matrix values, coordinate state, payload absence
reasons, or axes. Diagnostic flags are INFO/WARNING evidence and never become
exclusions on their own. Exact zero is retained without a pseudocount. Finite
extremes are summarized and may be warning-flagged, but remain measurements.

The named `standard` profile is threshold-free. It computes continuous
distribution, tail, rank/leverage, concentration, availability, zero, and
point-mass summaries without installing a binary outlier rule. It is always
diagnose-only.

The top-N point-mass table is accompanied by a separate global-extrema record.
For all finite available payloads it reports the exact minimum and maximum,
their exact counts, and their fractions of available payloads, even when an
extreme is absent from the frequency-ranked top-N table. Its interpretation is
always `exact_extreme_repetition_no_upstream_cause_inferred`; it does not infer
an optimizer cap, clipping, serialization, a biological boundary, or
measurement error.

Generic rules include none, quantile, IQR, MAD, absolute threshold, and an
identity-aligned user mask at whole-matrix, within-gene, or within-branch
scope. They are explicit optional `custom`-profile rules, not a recommended
generic recipe. Every threshold and rule is stored in the spec and diagnostic
provenance. Quantile rules must declare tie semantics: `strict` excludes cutoff
ties and `inclusive` includes cutoff ties. `exact_rank` is parked and fails
closed until implemented. MAD=0 is recorded as undefined rather than converted
to an infinite score.

The semantic boundary is frozen:

```text
outlier != error
long branch != anomalous cell
trim != diagnosis
```

Numeric trimming is not implemented. A future diagnostic-only, two-axis
empirical-extremeness record is parked with fields
`within_gene_percentile` and `within_branch_percentile`; QC1-FIX001 does not
compute them and assigns no threshold or admission effect.

Signed matrices use explicitly named absolute-value concentration variants.
Positive-only/log summaries count zero and negative values as outside the
metric domain; they never transform the source values.

`lt_apply_qc()` requires a stable decision ID, authority, and explicit selected
flag IDs, rule IDs, or mask. It only constructs a new Layer-2 declaration.
Selected eligible coordinates receive
`ineligible / measurement_qc_excluded`. Selected unresolved coordinates may be
resolved to that same state while their prior state is retained in application
provenance. A selected already-ineligible coordinate preserves its exact status
and reason. Selecting `not_applicable` fails closed. Application provenance
records the mutually exclusive eligible, already-ineligible, and unresolved
transition counts plus selected identities and prior-state ledgers.

## Parked visualization API

The future interface is recorded but not implemented in QC1:

```r
lt_plot_gene_rate(x, gene, data, tree, ...)
```

`data` may later be a gene-by-branch numeric matrix or compatible La Terra
object; `tree` will be a required first-class phylogenetic coordinate display
object. Values must join by explicit branch identity, never node/row order.
Trait overlays are optional. No plotting code is present in QC1.
