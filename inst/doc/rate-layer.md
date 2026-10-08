# La Terra V1 rate layer

The production baseline is the certified C3 nonnegative raw-scale working
mean. `ADD_LT`, `GBI`, and the historical positive-only `logGBI_strict` are
deterministic views of that one fit; rate construction does not refit or
rescale C3. The primary current log-relative object is constructed explicitly
from GBI with `lt_log_relative(gbi)`.

## Runtime validity and frozen identity

La Terra separates two questions:

1. **Structural/scientific validity.** Ordinary `validate_*()` calls enforce
   dimensions, stable scientific keys, coordinate and Layer-2 state machines,
   baseline/representation consistency, and dependency identities.
2. **Frozen-snapshot identity.** `lt_audit_fingerprint()` and
   `validate_lt_frozen_snapshot()` are opt-in certification/replay tools. A
   mismatch means only `NOT_THE_SAME_FROZEN_SNAPSHOT`; it does not by itself
   make an otherwise valid user object unusable.

An ordinary user-supplied `lt_matrix` does not need a pre-existing SHA-256.
When coordinate provenance is omitted, the constructor records explicitly
unfrozen runtime provenance and derives the ordered branch-ledger identity from
the supplied scientific keys. Display labels and non-authoritative annotation
metadata may be edited without invalidating a fitted rate. Scientific-key,
numeric, coordinate-state, Layer-2, or baseline-fit changes make dependent
results stale and require recomputation (or a future explicit compatible
rekey); they do not make the edited input matrix globally invalid.

```text
ADD_LT = Y - mu_C3
GBI    = Y / mu_C3
logGBI_strict = log(Y / mu_C3), conditional on Y > 0 and mu_C3 > 0
lt_log_relative = boundary-aware representation of GBI
```

Exact zero is valid C3 input. No pseudocount is used. With positive baseline,
an admitted observed numerator zero is a real ratio boundary, not a finite log
coordinate. The primary boundary object stores positive natural-log values,
the three-state ratio layer, and a distinct unavailable-cause ledger. A
post-division zero without certified numerator provenance fails closed as
numeric indeterminacy. The finite global-k10 coordinate survives only as the
explicit historical `logGBI_encoded()` view and is never the default
inferential meaning.

The new GBI is dimensionless. It is not historical legacy GBI, which used its
own mean-normalization construction and historical units. Equality could occur
accidentally in a special case but is not an identity or compatibility claim.

`lt_rate` preserves coordinate truth, measurement eligibility, baseline
support, representation availability, and representation reason as distinct
layers. Values and compact integer state matrices remain aligned to exact
ordered gene and branch ledgers. `lt_rate_ledger()` can expose one keyed row per
cell, but the compact object avoids a large character data frame in memory.

## Required C2 sensitivity

C2 is the independently defined positive-cell two-way log-OLS estimand. It
fits only finite, Layer-2-eligible, observed `Y > 0` cells and uses no
pseudocount. Connected positive components are solved and gauged separately.
Zero cells may receive a C2 baseline only if their gene and branch are supported
in the same positive component. No cross-component outer product is formed.
C2 outputs are always labelled `C2_positive_cell_log_OLS` sensitivity, never
production.

Exact keyed common-domain constructors are required for representation and
C3-vs-C2 comparisons. Unequal full domains are descriptive only and are not
estimator-performance evidence. The sensitivity layer selects no winner.

Outlier QC and trimming are not part of rate construction. Traits are not
inspected. RER uses a separate operator and is not implemented here.

## Trait-free geometry diagnostics

`lt_diagnose_rate()` produces a threshold-free `lt_rate_diagnostic` from one
already materialized rate object. It reports native-domain distribution,
neutral departure, exact finite-extreme repetition, representation loss,
gene/branch summaries, and a bounded keyed extreme ledger. It never changes a
rate value, availability mask, coordinate state, or Layer-2 decision.
Historical finite-log diagnostics remain replay-only. Current boundary
diagnosis derives exact-zero and positive supports from `ratio_state`; it never
infers state from equality to a sentinel.

Optional Y/baseline decomposition requires the caller to provide both the
source matrix and compatible certified fit. A stale dependency stops at the
operation boundary with `STALE_RECOMPUTE_REQUIRED`; the diagnostic never
rekeys, repairs, or refits an object. Repeated extrema are reported with the
neutral interpretation `exact_extreme_repetition_no_upstream_cause_inferred`.
Unusual is not error, extreme is not artifact, and long branch is not an
anomalous cell.
