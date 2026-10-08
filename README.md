# La Terra: A Matrix-First Framework for Branch-Wise Comparative Genomics

La Terra is a matrix-first framework and evolving R-native platform for
auditable branch-wise comparative genomics. Analyses start from a typed,
coordinate-valid matrix; scientific keys, value/state layers and provenance
remain explicit.

See [related work and recommended citations](#related-work-and-recommended-citations) for the methodological foundation, upstream software and related theory.

## Start with your own binary trait

The usable development route imports a supported SplitAlignerR exchange,
reads explicit trait/group/fold tables, computes a selected S2 recipe, and
reads its saved results. Neither recipe nor validation design is selected
for you. Inspect the available contracts with `lt_binary_recipes()` first.

```r
library(LaTerra)
input <- lt_import("/path/to/exchange.rds")
trait <- lt_read_trait("/path/to/trait.tsv", trait_name = "my_trait")
run <- lt_run_binary(input, trait, "/path/to/groups.tsv", "/path/to/folds.tsv",
                     recipe = "submitted_S2", validation_recipe = "nested_phase13",
                     output_dir = "/path/to/new_binary_run")
```

The trait table contains exact species keys and a named numeric 0/1 column;
optional 0.5 means explicitly excluded. Both endpoint classes must be present.
The grouping table supplies `species` and `genus`; the ordered fold table
supplies `fold_id`, `genus` and `seed`, including excluded-only genera.
No species-key normalization, trait recoding, grouping inference or fold
generation occurs. The import keeps the complete matrix/tree/coordinate/
upstream-evidence bundle, not only its numeric matrix.

```r
result <- lt_read_binary_run(run)  # or the saved run directory; verifies hashes
summary(result)
report <- lt_report_binary(result, "/path/to/new_binary_report")
lt_plot_binary(result, view = "oof")  # also "roc" or "trim"; requires ggplot2
```

The report and plots read saved outputs only: no ASR, screens or models are
refitted, and the run is not modified. Reports separate D1 input/QC counts,
D2 trimming and domains, D3 branch/screen support, and D4 saved validation,
folds, feature use and warnings. These are bounded diagnostics, not new QC
decisions or scientific certification. Failed/partial runs expose unavailable
dependencies explicitly as `not_run_or_incomplete`, never as successful results.

For an installed synthetic example, locate
`system.file("examples", "binary_small", package = "LaTerra")` and use its
`exchange.rds`, `trait.tsv`, `groups.tsv` and `folds.tsv`. The runnable
`examples/run_binary_small.R` script accompanies it; run it with
`source(system.file("examples", "run_binary_small.R", package = "LaTerra"))`.
It chooses new temporary output directories and prints their paths. Its 24 tips and 32 genes
are synthetic algebraic inputs exported by the real producer; they are not
biological evidence or a performance benchmark.

Imports currently require the exact SplitAlignerR **0.1.0.9002** backend for
schema 0.2.0-development. A default library may contain incompatible 0.1.0;
use an isolated R library and check the loaded version before importing.
Computation requires `ape`, `castor` and `glmnet`; plotting adds `ggplot2`.
See [the binary workflow guide](inst/workflow/BINARY.md) for installation
preflight, table formats, recipes and result boundaries. This is development
software; completed computation does not imply release certification.

## Historical Marine/S2 exact-replay workflow (M1 development)

The study-specific entry point now computes the frozen S2 arithmetic baseline
and GBI, deterministic branch annotations, and the marine/aquatic baseline
Welch screens from raw inputs. It writes outputs and an audited `lt_run`:

```r
library(LaTerra)
run <- lt_run_marine(data_root = "/path/to/marine_inputs",
                     output_dir = "/path/to/new_marine_run")
run$results$screens
```

Supply these three files beneath `data_root` (exact paths and snapshot hashes
are available through `lt_marine_profile()`):

```text
marine_inputs/
  00_traits_and_species/trait_table.mammal302.active_TY_NK_final_18pt.tsv
  17_large_matrix_archives/branch_coordinate_matrices.tar.gz
  supplemental_authority/marine_2026_fixed_coord_state_17432_v1.tsv.gz
```

The first two are the frozen Dryad inputs. The third is the pre-existing A2U
artifact admitted as TARGET001-AUTHORITY-ADDENDUM001. It supplements
state provenance only; the original raw values remain numerical authority,
and the historical 17,428-gene classified member retains its own partial axis.
No precomputed GBI or screen results are needed or read by this entry point.

ADDENDUM006 restores the historical execution path automatically: use the
supplied crosswalk as an inverse computational-order view, compute the unchanged
S2 recipe, write with the original R writer, reorder only the text tokens back
to the public axes, then read this newly generated GBI. Raw values and public
branch identities do not change. No manual conversion, Python runtime, archived
GBI input or precision option is involved. This historical order/serialization
rule does not apply to the ordinary other-data S2 entry point.

The output parent must exist and the run directory must be new and outside
the input root. No files are overwritten and no historical absolute-path
fallback is used. The default frozen environment is R 4.4.2, ape 5.8-1 and
castor 1.8.4 on aarch64-apple-darwin20. An explicit
`environment = "compatibility"` run is warned and non-certifying.

Major outputs: `raw_matrix.rds`, `upstream_tokens.rds`,
`S2_baseline_GBI.rds` (its `gbi` is the generated/read-back representation),
`GBI.generated_native.tsv`, `GBI.generated_consumed.oldlabels.tsv`,
`GBI_computational_order.tsv`, `GBI_consumption_boundary.rds`,
per-trait branch states/tested/significant/gene-ledger
tables, `run.rds`, `input_hashes.tsv`, `frozen_profile.yaml`, `sessionInfo.txt`,
`RUN_STATUS.yaml` and `SHA256SUMS`. Failures after output reservation retain
their actual stage in `RUN_STATUS.yaml`; a completed run is not automatically
a certified replay. The separate development validator compares results to
the frozen TARGET001 outputs after execution.

This historical study profile includes its frozen inclusive upper-2.5% trim
and strict-prefix FDR rule; neither is a recommended generic default. SHA
checks belong to this explicit frozen replay route, not normal user objects.
M1 alone does not execute the later models or confer full 22-target V1 certification.

## Complete Marine computation and reporting (development)

The M2 entry extends the same raw-input workflow through permutation replay,
grouped validation, full-data fits, turnover and descriptive projections. It
requires the thirteen-file admitted bundle, whose layout is documented in
[the installed Marine workflow guide](inst/workflow/MARINE_M2.md).

```r
run <- lt_run_marine_m2(
  data_root = "/path/to/marine_inputs",
  output_dir = "/path/to/results/marine_run_001",
  cores = 4L
)
```

`lt_export_marine_tables()` is a separate terminal reporting entry for a
completed run. It does not rerun the expensive scientific stages or import
old result workbooks. The five table files preserve their historical worksheet
schemas; separate companions retain authority, transformations, historical
references and the 22-target evidence ledger. Reporting-only dependencies are
the R packages `openxlsx`, `xml2` and `zip`.

```r
tables <- lt_export_marine_tables(
  data_root = "/path/to/marine_inputs",
  run_dir = "/path/to/results/marine_run_001",
  output_dir = "/path/to/results/marine_tables_001"
)
tables$run$results$target_status
```

Only the final S5 coefficient display uses the explicitly admitted six-
significant-digit numeric view. Fits, signs, markers and ordering retain full
precision. Original-generator-unlocated presentation rules are disclosed;
no old numeric answers are copied. TGT021 table reproduction remains under
independent development review. Neither a completed model run nor assembled
tables certify the development package.

## Retained development capabilities

This repository also retains the earlier lightweight S3
scientific objects, authoritative R-native structural validators,
machine-readable schemas, provenance containers, certified C3-derived rate
views, branch-state association/calibration operators, and a boundary-aware
log-relative representation.

`lt_matrix` stores dense base-R values, coordinate truth, current-payload
absence reasons, and typed payload identity as separate fields. Imported and
derived branch-associated quantities are explicit and coordinate identities
are stored rather than recomputed.

The current log-relative API uses `lt_log_relative(gbi)` and separates
ZERO_MASS from POSITIVE_LOG without a joint p-value. Historical finite
`logGBI` rates remain readable and migratable but cannot enter ordinary
inference. Their global-k10 numerical representation is a diagnostic/replay
view only: ordinary inference refuses it, structural validation does not
authenticate its provenance, and frozen replay comparison requires separately
supplied source and replay authorities. La Terra does not perform upstream
alignment, tree inference, SplitAligner coordinate construction, or new
scientific null design. The explicit Marine/S2 route above consumes supplied
coordinate identities and performs only the frozen historical ASR recipe.
Retained C3 and zero-inference development APIs are not substitutions for the
current Marine/S2 V1 contract.

The frozen architectural note is installed at
`inst/doc/architecture-freeze.md`.

## Related work and recommended citations

La Terra builds on the gene-branch interaction approach of Wu, Yonezawa and
Kishino (2017) and consumes branch-coordinate matrices produced by
SplitAligner. The related graph-theoretic work below provides a formal
account of cross-gene branch identity under fixed labelled inputs and
conventions.

### Primary recommended citations

For analyses using La Terra with SplitAligner-derived matrices, we recommend
citing both the foundational method and the upstream branch-mapping work:

1. **Wu, J., Yonezawa, T., and Kishino, H. (2017).**
   Rates of Molecular Evolution Suggest Natural History of Life History
   Traits and a Post-K-Pg Nocturnal Bottleneck of Placentals.
   *Current Biology*, **27**(19), 3025-3033.e5.
   [https://doi.org/10.1016/j.cub.2017.08.043](https://doi.org/10.1016/j.cub.2017.08.043).
   This is the methodological foundation for the gene-branch interaction
   approach and rate-based trait reconstruction.

2. **Wu, J. (2026).**
   SplitAligner: A Gene-Species Tree Reconciliation Framework Using
   Split-Based Branch Mapping.
   *bioRxiv* preprint.
   [https://doi.org/10.64898/2026.02.24.707838](https://doi.org/10.64898/2026.02.24.707838).
   This describes the upstream branch-mapping framework used to construct
   coordinate-aware matrices.
   [Software repository](https://github.com/wujiaqi06/SplitAligner).

### Related theoretical work

3. **Wu, J. (2026).**
   A Unique Graph-Theoretic Truth Table for Cross-Gene Branch Identity.
   *bioRxiv* preprint.
   [https://doi.org/10.64898/2026.07.22.740066](https://doi.org/10.64898/2026.07.22.740066).
   This develops the graph-theoretic basis for a unique coordinate-and-state
   ledger under fixed labelled inputs and conventions. Cite it when
   discussing that theoretical foundation.

For reproducibility, also report the exact La Terra version or commit and
the upstream software versions used.
