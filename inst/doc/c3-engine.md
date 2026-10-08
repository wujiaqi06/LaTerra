# C3 domain, existence, and baseline engine

## Scope

The T1A engine accepts an `lt_matrix` and a separate layer-2 measurement-fit
eligibility declaration. Its fit domain is exactly

```
coord_state == "observed" &
measurement_fit_eligibility == "eligible" &
is.finite(Y) & Y >= 0
```

Exact numeric zero is admitted. There is no pseudocount, semantic-zero
tolerance, trimming, winsorization, magnitude rule, trait field, or
dataset-specific QC operator. When no exclusion ledger is supplied,
`eligible_declared` is explicitly recorded as a caller-domain declaration and
`upstream_qc_certified` remains false. Diagnostics are audit-only inputs and
cannot mutate this layer.

Layer 2 is scoped only to observed coordinates. Every coordinate whose truth
state is not `observed` has status `not_applicable` and reason
`coordinate_not_observed`; it is neither unresolved nor excluded measurement
evidence. Conversely, an observed coordinate cannot use `not_applicable`.

## Graph and boundary sequence

The implementation builds the bipartite gene-branch graph on admitted cells,
ledgers degree-zero vertices, assigns deterministic connected-component IDs,
and computes branch-outer/gene-inner compensated pairwise margins. A component
whose exact total is zero is `ZERO_TOTAL_COMPONENT`. A positive-total vertex
whose exact margin is zero is `ZERO_MARGIN_BOUNDARY`; its incident fitted means
are exact zero and its scope is removed before interior certification.
Components are recomputed after every support reduction.

## Exact minimal-face certificate

The marginal-cone certificate uses the transportation-flow structure directly
and does not treat an LP solver status as authority.

The admitted nonnegative payload vector is an exact rational feasible flow
because every binary64 input is a dyadic rational and the target margins are
defined from the same edge ledger. On a positive-margin candidate component:

1. every allowed edge has a forward residual arc `gene -> branch`;
2. every edge with `Y > 0` also has a reverse residual arc
   `branch -> gene`;
3. positive-edge undirected components are collapsed;
4. SCCs of the resulting residual graph are computed exactly from integer
   edge identities; and
5. an edge is retained precisely when its endpoints are in the same residual
   SCC.

For each retained edge, a residual directed cycle provides a feasible cycle
augmentation with that edge positive. Averaging the original feasible flow
and these cycle-augmented flows is a symbolic exact strictly-positive primal
witness on all retained edges. Thus one residual SCC certifies relative
interior without a floating feasibility tolerance.

For a forced-zero edge, the SCC condensation is a DAG. Longest-path integer
ranks give

```
h_gene   = -rank(SCC_gene) / S
h_branch =  rank(SCC_branch) / S
```

where `S` is the sum of all forced-edge rank differences. Every candidate
edge has nonnegative `a_e^T h`, every retained edge has zero, every forced edge
has a strictly positive integer numerator, their normalized sum is one, and
`t^T h` is exactly zero: positive-payload edges have coefficient zero and all
positive-coefficient edges have payload zero. This is the exact dual exposing
certificate for the minimal face. A positive payload on a forced edge is a
fatal contradiction.

On the graph representation, the residual-SCC and condensation facial-
certificate core is `O(V + E)` time and `O(V + E)` working memory. This bound
does not describe the complete current dense-R pipeline: construction and
validation of dense matrix/mask storage may scan `G x B`. A declared
deterministic resource ceiling fails closed as `NUMERICALLY_INDETERMINATE`; it
is never relaxed until a convenient state appears.

## C3 numerical fit

Only final certified positive-support components enter deterministic IPF.
Updates solve the actual incomplete-mask margin equations. Fitted means are
materialized only on ordered admitted edges. Boundary edges receive exact zero;
matrix holes and cross-component products are never constructed as baseline
cells.

Final margins are recomputed with the frozen compensated pairwise merge tree
and checked using `fit_rtol = 1e-10` and
`fit_atol = 1e-12 * max(1, T / n_edges)`. The stable zero-aware objective is
the generalized-KL/Bregman working measure; it is not a Poisson probability,
variance, likelihood-ratio, or inferential model.

Each final positive component uses the gauge
`sum(gene_log_factor) = 0`. Factors are secondary and component-specific;
supported `mu` is the gauge-invariant scientific baseline identity. Full
same-environment solver replay must reproduce state, iteration count, factor
and `mu` hashes. Any margin, objective, gauge, replay, overflow, underflow,
denominator, support, or positive-Y/zero-mu failure holds the affected run and
releases no partial production baseline.
