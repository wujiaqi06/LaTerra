# Boundary-aware log-relative representation and historical migration

## Current V1 contract

`lt_log_relative(gbi)` constructs the primary `lt_log_relative` object with
machine representation ID `log_relative_boundary`. Its authoritative cell
payload is `positive_log_values`, `ratio_state`, and
`unavailability_reason`. The mathematical ratio state has exactly three
values: `POSITIVE`, `EXACT_ZERO_BOUNDARY`, and `UNAVAILABLE`. Layer-2
eligibility remains a separate keyed dependency.

“Exact zero” means a numeric exact zero preserved in the admitted observed
payload under the frozen input/eligibility contract. It does not prove an
ontological biological evolutionary rate of exactly zero, that no tiny
positive quantity existed before admitted serialization, or an
optimizer-estimated strict zero.

An exact-zero boundary is admitted only when authoritative source evidence
shows an observed numerator exact zero, a finite strictly positive denominator,
valid coordinates, and passing measurement/Layer-2 gates. Post-division
`GBI == 0` alone is insufficient: an uncertain or underflow zero becomes
`UNAVAILABLE / numeric_indeterminate` or fails closed. `observed_payload_unavailable`
is never relabelled as `baseline_unavailable`.

Positive cells store natural `log(GBI)`. Exact-zero boundaries have no finite
mathematical log coordinate. Unavailable cells have no ratio/log-relative
coordinate. Default binary analysis uses boundary-aware two-component
estimands: ZERO_MASS (`delta_z`) and POSITIVE_LOG (`delta_L_positive`). Each
has its own estimability and Null A/Null B calibration. No joint p-value is
defined. La Terra V1 does not claim to implement a full generative hurdle model.

## Explicit historical replay

`logGBI_encoded(object, mode = "legacy_global_k10")` is the sole finite
historical view. It preserves the original source-GBI `m_ref`, assigns
`log(m_ref) - 10` only in the derived view, and carries
`default_continuous_inference_authorized = FALSE`. It may be materialized for
diagnosis and frozen replay comparison, but it is not a La Terra V1 inferential
representation. No gene-local or optimized sentinel is a production mode.

Exact-key `gene_ids` and `branch_ids` may restrict this derived view for a
historical replay. Such restriction never recomputes `m_ref`, `L_min_plus`, or
the sentinel: they remain bound to the complete source authority and its
identity ledger. The derived view also carries the compact source ratio state.

`validate_lt_logGBI_encoded()` checks local structural consistency only. Its
self-hashes detect accidental corruption; they do not authenticate provenance.
`validate_lt_logGBI_encoded_source(encoded, source_log_relative)` performs the
scientific source-bound check by independently deriving axes, ratio states,
positive logs, full-source `m_ref`, sentinel, and source identities from the
separately supplied boundary object. A coordinated payload change plus local
rehashing therefore fails against the fixed source.

`validate_legacy_logGBI_replay()` additionally compares the source-bound view
to a separately frozen external replay authority. It computes no association,
calibration, p-value, or adjusted p-value. If both the encoded object and its
purported source are changed and rehashed, they describe a different source;
package-local self-hashes cannot validate that pair against the original frozen
authority. La Terra uses no hidden token, namespace seal, or secret for this
purpose.

`lt_rate_spec(method = "logGBI")` and ordinary finite-log construction return
a classed guidance error. Existing historical RDS objects remain readable as
historical objects, but only `lt_migrate_logGBI()` may construct the new
scientific object. Migration requires compatible source GBI authority, never
infers boundary state from sentinel equality, and never overwrites the old
object.

Historical finite `lt_rate(representation = "logGBI")` objects and
`lt_logGBI_encoded` views fail closed at every ordinary domain, association,
calibration, adjustment, and internal score entry point. Acknowledgement flags,
edited provenance/configuration, or direct `:::` calls cannot unlock inference.
Boundary-aware inference is parked outside La Terra V1; use GBI for V1 ordinary
inference.

New boundary `semantic_identity` and `audit_identity` values execute the named
canonicalization literally: their normalized payloads are serialized with R
serialization version 3 before SHA-256 hashing. Ordered axes, compact raw-code
order, NA placement, and reason/state vocabularies are authoritative. This
does not change any historical package-wide hash identity.

Display labels and non-authoritative metadata are excluded from runtime
scientific identity. Scientific dependency changes make dependent results
stale and require recomputation or explicit compatible rekeying.
