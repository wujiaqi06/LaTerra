# La Terra: A Matrix-First Framework for Branch-Wise Comparative Genomics

La Terra explores how molecular evolution varies across genes and branches,
and how those patterns relate to organismal traits. Built for branch-length
matrices from [SplitAligner](https://github.com/wujiaqi06/SplitAligner), it
brings gene-level association screens, sparse predictive models and
diagnostic plots into a single R workflow.

The analysis starts from a complete SplitAlignerR exchange: a gene × branch
matrix, its species tree, branch-coordinate keys, cell states and provenance.
These keys preserve each reference branch's identity across genes, while
explicit states distinguish unavailable coordinates from observed values.
You supply the trait table and, for grouped validation, the group and fold
tables.

La Terra is under development. The workflow described here supports binary
traits and saves the input snapshots, analysis settings and software
information needed to trace its results.

## Installation

Use R 4.2 or later and start a fresh R session. The current importer requires
**SplitAlignerR 0.1.0.9002** with exchange schema **0.2.0-development**.
The [SplitAlignerR repository](https://github.com/wujiaqi06/SplitAlignerR)
currently has a different version on its default branch. A public source for
the required backend has not yet been established for these instructions;
installation therefore requires that you already have the matching
SplitAlignerR and LaTerra source archives.

Install into a separate library, replacing the archive paths with the files
you obtained:

```r
lib <- file.path(path.expand("~"), "R", "LaTerra-dev")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

install.packages(c("digest", "yaml", "ape", "castor", "glmnet", "ggplot2"),
                 lib = lib)
install.packages("/path/to/SplitAlignerR_0.1.0.9002.tar.gz",
                 repos = NULL, type = "source", lib = lib)
install.packages("/path/to/LaTerra_source.tar.gz",
                 repos = NULL, type = "source", lib = lib)

stopifnot(as.character(packageVersion("SplitAlignerR")) == "0.1.0.9002")
find.package("SplitAlignerR")
library(LaTerra)
```

`parallel` is included with R. Computation uses `ape`, `castor` and `glmnet`;
`ggplot2` is needed for plots. Marine spreadsheet export additionally requires
`openxlsx`, `xml2` and `zip`. Building the upstream backend from source may
require the compiler toolchain specified by that backend.

Install LaTerra before running an analysis; source sessions created with
`pkgload::load_all()` are not supported. See the
[binary workflow guide](inst/workflow/BINARY.md#install-and-check-the-actual-runtime)
for runtime checks and library selection.

## Quick start: your own binary trait

Begin with a SplitAlignerR exchange and three tab-separated tables:

| Input | Contents |
| --- | --- |
| `exchange.rds` | Matrix, species tree, coordinates, cell states and provenance, exported together by SplitAlignerR. |
| `trait.tsv` | Species names and a named numeric trait column: 0 and 1 for the two classes; optional 0.5 marks excluded species. |
| `groups.tsv` | Species names and the groups used for validation, in columns `species` and `genus`. |
| `folds.tsv` | The ordered validation folds, with columns `fold_id`, `genus` and `seed`. |

Use species names exactly as they appear in the exchange. Both trait classes
must be represented, and the fold table must cover every group, including
those containing only excluded species. Prepare the coding, groups and fold
order before running; La Terra uses these choices as supplied.

Choose an analysis recipe and a validation design with `lt_binary_recipes()`.
The example below uses `submitted_S2` with `nested_phase13`; their meanings
are explained in [Analysis recipes](#analysis-recipes).

```r
library(LaTerra)
input <- lt_import("/path/to/exchange.rds")
trait <- lt_read_trait("/path/to/trait.tsv", trait_name = "my_trait")
run <- lt_run_binary(input, trait, "/path/to/groups.tsv", "/path/to/folds.tsv",
                     recipe = "submitted_S2", validation_recipe = "nested_phase13",
                     output_dir = "/path/to/new_binary_run")
```

Keep the complete object returned by `lt_import()`: the tree, coordinates
and provenance are part of the input. The parent output directory must
already exist, and each run must use a new directory.

The analysis represents branch-length variation as **GBI** (gene–branch
interaction), using the baseline defined by the selected recipe. Gene
screens and predictive models then assess its relationship with the trait;
a GBI value alone is not a test of association.

## Read, inspect and plot the results

```r
result <- lt_read_binary_run(run)  # or the saved run directory; verifies hashes
summary(result)
report <- lt_report_binary(result, "/path/to/new_binary_report")
lt_plot_binary(result, view = "oof")  # also "roc" or "trim"; requires ggplot2
```

Reports and plots read the saved results without refitting models or
changing the analysis. They cover inputs and QC summaries, trimming and
computable values, branch states and gene screens, and validation results
with fold-level feature use and warnings. These sections are labelled
D1–D4 in the report. Diagnostics help you inspect the analysis; they do not
automatically change which observations are included. For a failed or
partial run, results that could not be produced are shown as
`not_run_or_incomplete`.

To try the workflow with a small synthetic dataset, run:

```r
source(system.file("examples", "run_binary_small.R", package = "LaTerra"))
```

The example contains 24 tips and 32 genes, exported through SplitAlignerR.
Its four input files are installed under
`system.file("examples", "binary_small", package = "LaTerra")`. The script
creates new temporary run and report directories and prints their paths.
This is a worked software example, not biological evidence or a performance
benchmark.

See [the binary workflow guide](inst/workflow/BINARY.md) for installation
preflight, table formats, recipes and result boundaries.

## Analysis recipes

A recipe records the scientific choices used in an analysis. Select it
explicitly so that the same choices can be inspected and repeated.

`submitted_S2` preserves the analysis recipe of the 2026 marine-mammal study.
Within each gene, values at or above the 97.5% quantile (R type 7) are excluded.
For the remaining values, GBI is calculated by dividing each branch length
first by its branch arithmetic mean and then by its gene arithmetic mean.
Both means are calculated after trimming; no grand-mean normalization is
applied. The recipe also specifies the historical ancestral-state
reconstruction and Welch screen with its strict-prefix FDR rule.

GBI describes variation relative to this recipe's baseline. Values above or
below 1 do not, by themselves, establish a biological rate shift. Observed
zeros remain distinct from unavailable input cells; nonfinite calculated
ratios are stored as `NA`.

Choose the validation design separately:

- **`nested_phase12B` and `nested_phase13`** repeat feature screening within
  each training fold. They preserve two historical arithmetic variants and
  remain separate choices.
- **`global_screen_foldwise`** fits the globally screened gene set fold by
  fold. Because screening has already used the full dataset, its held-out
  predictions do not provide nested feature-selection validation.

These recipes reproduce specified analyses. Their trimming and screening
rules should be assessed for the dataset and scientific question at hand.

## Historical Marine/S2 replay (worked example)

The Marine workflow reproduces the study-specific downstream analyses from
frozen inputs. Use `lt_run_marine()` for the baseline and gene screens,
`lt_run_marine_m2()` for the full computation, and
`lt_export_marine_tables()` to export tables from a completed run without
refitting. The [Marine replay guide](inst/workflow/MARINE_REPLAY.md) describes
the baseline inputs and output files; the
[Marine full-workflow guide](inst/workflow/MARINE_M2.md) covers the complete
input bundle and computation. `?lt_run_marine` and `lt_marine_profile()`
also describe the baseline entry and its required inputs.

Frozen input snapshot checks and historical computation-order/serialization
rules are specific to this Marine replay route. The trimming, GBI, ancestral
annotation and screening rules of `submitted_S2` also apply when that recipe
is explicitly selected for other datasets. Completion is recorded separately
from comparisons with the historical reference results.

## Scope and supported environments

Sequence alignment and tree inference precede SplitAligner. SplitAligner
constructs branch-coordinate matrices from those trees; La Terra analyzes
the exported measurements alongside the trait and tree information.

Development checks have been performed on macOS arm64 with R 4.4.2.
The Marine exact-replay reference also specifies ape 5.8-1 and castor 1.8.4
on aarch64-apple-darwin20. For the Marine entry points,
`environment = "compatibility"` permits an attempt in another environment
and records a warning. Other platform/version combinations still require
validation; this mode does not guarantee successful execution or equality
with the historical outputs.

Completed runs retain input-hash manifests, `SHA256SUMS` and `sessionInfo.txt`.
If an analysis fails after creating its run directory, `RUN_STATUS.yaml`
records the stage and the files already produced remain available for
inspection.

## Additional interfaces

The package also includes typed matrix objects, structural validators and
provenance utilities, together with separate development interfaces for a
C3 baseline and log-relative analysis. These interfaces are documented in
`inst/doc/` and are not substituted into the binary-trait workflow.
The log-relative interface treats exact zero mass and positive log variation
separately, without a joint p-value. Historical finite-sentinel logGBI objects
can be read or migrated, but cannot enter ordinary inference. See the
[architecture guide](inst/doc/architecture-freeze.md) for details.

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
