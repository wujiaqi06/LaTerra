# La Terra V1 branch-state and association layer

`lt_trait` remains a taxon-level container. `lt_branch_state` is a distinct
first-class object whose identifiers are exact branch scientific keys. T2A
does not convert taxon traits into branch states and does not perform ancestral
state reconstruction.

## Layer 5

The association domain is downstream of coordinate truth, payload
availability, measurement eligibility, and rate representation availability.
For a gene-branch cell it records exactly one of:

```text
eligible
rate_unavailable
trait_unavailable
branch_user_excluded
```

Native domains use one rate representation. Common domains intersect every
requested representation before association, so cross-representation
comparisons never confuse unequal per-gene branch domains with rate geometry.
Scientific keys are matched exactly. The default warns and proceeds on the
exact intersection; strict mode errors. Position-only, fuzzy, case-folded, and
display-label matching are forbidden.

## Point estimands

For continuous state `x` and rate response `r`, the per-gene point estimate is
the ordinary unweighted slope with intercept on the exact Layer-5 domain:

```text
beta = sum((x - mean(x)) * (r - mean(r))) / sum((x - mean(x))^2)
```

For explicit binary reference/focal coding:

```text
intercept = mean(rate | reference)
beta      = mean(rate | focal) - mean(rate | reference)
```

The sign is always focal minus reference. No lexical level ordering is used.
No automatic centering, standardization, weighting, trimming, winsorization,
pseudocount, imputation, or rate reconstruction occurs.

Every gene remains in the output. Mathematical non-estimability is represented
by an explicit status rather than row deletion. Tiny but mathematically
estimable domains return their estimate with `small_domain_warning`.

## Inference boundary

The response is rate and the predictor is branch state. This is association,
not prediction or causality. T2A releases point estimates and descriptive
correlations only:

```text
inference_status = POINT_ESTIMATE_ONLY
significance_calibration = NOT_PERFORMED
```

No ordinary branch-iid p-value, confidence interval, q-value, or FDR is exposed.
Phylogenetic/null calibration is a later scientific stage.

Historical finite-sentinel `lt_rate(representation = "logGBI")` objects are
read/migrate-only and cannot enter native/common association-domain or
ordinary univariate-association routes. Use `lt_migrate_logGBI()` for the
current boundary object. The `legacy_global_k10` view is diagnostic/replay-only
and is refused by every ordinary association route, including internal helpers.
Source-bound replay validation computes no association statistic. Use GBI for
La Terra V1 ordinary inference.
