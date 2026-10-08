# LT001-A1A architecture hardening

Status: architecture hardened; scientific algorithms not implemented.

Base scientific-object commit:
`e0cff101a88723fa2793bb930f9d1cd60f854312`.

## Frozen decisions

- La Terra V1 is R-native and MIT licensed. Third-party source code must not be
  copied into the package.
- The first-class input is a typed, coordinate-valid, branch-associated
  matrix. Branch length is the main V1 raw payload, not the only valid payload.
- Imported external rate representations are permitted when their payload
  origin and provenance are explicit; La Terra must not imply that it derived
  them.
- Dense base-R matrices are the canonical V1 storage representation. There is
  no backend abstraction, DelayedArray, Arrow, Python, Rcpp, or benchmark.
- Lightweight S3 objects and function-style transformations remain frozen. R6
  mutable state is excluded.
- Alignment, tree inference, branch-length estimation, and SplitAligner
  coordinate construction remain upstream of La Terra core.

## Layering

```text
scientific objects
  -> pure operators
  -> configuration recipe
  -> audited run object
```

A1A hardens scientific and audit containers. Pure scientific operators remain
explicit `not implemented` placeholders.

## Matrix contract

An `lt_matrix` contains:

- `values`: dense numeric gene-by-branch matrix;
- `coord_state`: coordinate-truth mask;
- `value_reason`: reason an observed coordinate lacks a value for the current
  payload;
- `payload`: typed identity of the current branch-associated quantity;
- ordered `gene_ids` and `branch_ids`;
- exact `coordinate_provenance`;
- optional metadata.

The coordinate-state vocabulary is exactly `observed`, `NA_struct`, `NA_fuse`,
`NA_topo`, and `residual_NA`. It is not extended with calculation-induced
missingness. Downstream operators must not rewrite coordinate truth merely
because they filter a value or cannot define the current payload.

The value invariants are:

```text
coord_state != observed
  => values is NA and value_reason is NA

coord_state == observed and values is finite
  => value_reason is NA

coord_state == observed and values is NA
  => value_reason is a non-empty method-specific reason
```

Numeric zero is a valid finite value. NaN and infinite values are prohibited.
No global downstream-reason vocabulary is frozen; a future scientific method
must version and validate its own codes.

## Typed payload

Every matrix payload records `type`, `name`, `units`, `scale`, `origin`, and
`spec_id`. Origin is `imported` or `derived`. A derived payload requires a
non-empty specification identity. An imported payload may have a null spec ID,
or may identify the external specification that produced it.

Rate operators will later declare accepted input payload types. No such
operator is implemented in A1A.

## Coordinate provenance

The minimum coordinate-axis identity records:

- coordinate system;
- reference-tree SHA-256;
- ordered branch-ledger SHA-256;
- ordered taxon-ledger SHA-256;
- branch-label contract;
- source or mapping method identity.

Software-derived mappings may additionally record source software, version,
and commit. External matrices are not required to originate from SplitAligner.
Coordinates are stored and validated, never recomputed by `lt_matrix`.

For the legacy recipes, ordered-ledger hashes use UTF-8 identifiers, one per
line, LF endings, and a final LF. Branches are B1-B601 in matrix/branch-map
order; taxa are the 302 terminal `sub_tree` identifiers in branch-map order.

## Input and configuration contract

Generic matrix inputs support TSV, CSV, RDS, and archive members. An archive is
not mandatory. A complete external matrix may explicitly declare
`all_observed`; otherwise a matrix with coordinate absence must declare a
separate coordinate-state input. Current-payload reason input is independently
declared.

Scientific stages use the generic module shape:

```text
enabled
method
contract_version
parameters
```

The generic schema assumes no particular model, engine, family, penalty, or
null analysis. Null models are an optional list of modules. Method-specific
semantic validators may be added only when those scientific contracts are
separately frozen.

R-native semantic validation is authoritative in V1. The YAML/JSON Schema is a
machine-readable interoperability contract; no JSON-Schema runtime dependency
is required.

## Domain invariants

Ineligible taxa require a non-empty exclusion reason. Eligible taxa must have
an `NA` exclusion reason. Constructors do not infer reasons.

## Audited run and provenance

The provenance container records input and output hashes, configuration hash,
software and environment identity, ordered gene/branch/taxon ledgers,
eligibility masks, specification and filtering ledgers, RNG identity and seed,
taxon-level fold membership, and commands.

`specification_ledger` preserves exact formula/specification identity and
version. `filter_ledger` preserves ordered decisions and reasons.
`fold_ledger` records `fold_id`, `taxon_id`, `role`, and `group` so membership is
machine-auditable. `seed_ledger` includes RNG kind fields and R version.

Run status is computational state only. A completed `lt_run` never implies
scientific certification, and La Terra does not self-certify PASS.

## Marine legacy authority by role

The legacy authority set is not collapsed into a single competing artifact:

- scientific specification: `marine_mammal_paper.NC_revise.pdf`;
- executable reference: reviewer/editor software package v1.2;
- data/result authority: frozen Dryad archive;
- source/code identity: exact Git tag object and peeled commit.

Both legacy recipes record the exact identities and hashes. They remain golden
external references; their source code is not copied into La Terra.

## Explicitly absent from A1A

- SplitAligner mapping or branch-coordinate reconstruction;
- GBI, additive GBI, or log-GBI calculations;
- ancestral-state reconstruction;
- single-gene statistics or FDR procedures;
- imputation, scaling, LASSO, cross-validation, prediction, or null models;
- numerical golden-result tests;
- backend abstraction, Python, Rcpp/C++, benchmarking, UI, or manuscript prose.

## Deferred, non-blocking method contracts

- Method-specific `value_reason` code vocabularies and validators are deferred
  to the work order that freezes each scientific operator.
- Accepted payload types are deferred to each rate/operator specification.
- Scientific certification remains an external review decision rather than a
  package-generated status.

There are no unresolved blocking A1A architecture questions.
