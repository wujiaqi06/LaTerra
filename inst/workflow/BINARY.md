# Binary-trait workflow: explicit inputs, recipes and saved results

La Terra: A Matrix-First Framework for Branch-Wise Comparative Genomics is
an evolving R-native framework/platform. The present usable route wraps the
existing historical S2 operators; it does not turn their scientific choices
into universal defaults. The Marine/S2 replay below the generic route in the
README is a historical example, not a future scope limit.

## Install and check the actual runtime

Install the matching development packages using the
[public installation instructions](../../README.md#installation). They download
SplitAlignerR 0.1.0.9002 from its development prerelease, verify its SHA256,
and install a fixed La Terra revision into a separate library. The backend
source archive is available at the
[development prerelease](https://github.com/wujiaqi06/SplitAlignerR/releases/tag/v0.1.0.9002-laterra-dev).

After installation, check the selected runtime in a fresh R session:

```r
lib <- file.path(path.expand("~"), "R", "LaTerra-dev")
.libPaths(c(lib, .libPaths()))
stopifnot(as.character(packageVersion("SplitAlignerR")) == "0.1.0.9002")
stopifnot(identical(SplitAlignerR::splitaligner_exchange_schema()$schema_version,
                    "0.2.0-development"))
find.package("SplitAlignerR")
library(LaTerra)
```

The import contract is schema 0.2.0-development with the matching public
producer validator. An installed 0.1.0 backend is not interchangeable. Start
a fresh R session after changing libraries if the wrong namespace was already
loaded. `ape`, `castor` and `glmnet` are computation dependencies; `ggplot2`
is needed only for the saved-result plots. Installation success is not
release certification or cross-platform validation.

Audited S2 runs require the installed package's `R/LaTerra.rdb` runtime.
`pkgload::load_all()` source sessions are not supported for runs: the run
preflight rejects them before annotation, fitting or output-directory creation
and explains how to install the source package. It does not borrow a runtime
hash from a different installed copy or invent one for developer mode.

## Supply the complete keyed inputs

`lt_import("exchange.rds")` is the generic name for the accepted
`lt_read_splitaligner_exchange()` import. Keep the complete return value:
matrix, tree, coordinates, original exchange and provenance. Original arrays,
states/reasons and label evidence are not replaced by a guessed crosswalk.
La Terra does not align sequences, infer trees or construct upstream branch
coordinates. Bare B-label matrices and unsupported exchange schemas are not
silently repaired.

The TSV readers use headers and explicit columns:

| File | Required columns | Contract |
| --- | --- | --- |
| `trait.tsv` | `species` and your named trait column | Numeric 0 and 1; optional 0.5 is excluded; no missing values; both endpoints present. |
| `groups.tsv` | `species`, `genus` | Exact supplied species keys and grouping keys; no parsing of species names. |
| `folds.tsv` | `fold_id`, `genus`, `seed` | Complete ordered fold ledger, including excluded-only genera; explicit integer-valued seeds. |

Read a trait with `lt_read_trait(file, trait_name = "my_trait")`. If the key
column has another name, specify `taxon_column`; keys themselves are never
inferred or normalized. Groups and folds can be read explicitly with
`lt_read_groups()` and `lt_read_folds()`, or supplied as TSV paths to the run
wrapper. Existing `lt_trait` and data-frame objects are also accepted.
The run validates the exact keyed trait/tree/group domains and fold design
before computation. No automatic discretization, fold generator or
measurement-eligibility decision is introduced.

## Inspect and choose recipes, then run

```r
lt_binary_recipes()
input <- lt_import("/path/to/exchange.rds")
trait <- lt_read_trait("/path/to/trait.tsv", trait_name = "my_trait")
run <- lt_run_binary(input, trait, "/path/to/groups.tsv", "/path/to/folds.tsv",
                     recipe = "submitted_S2", validation_recipe = "nested_phase13",
                     output_dir = "/path/to/new_binary_run", cores = 1L)
```

`input` may instead be the supported exchange RDS path. If `trait` is a TSV
path rather than an `lt_trait`, select its column explicitly with
`trait_name`. The output parent must already exist and the run directory must
be new; no existing run is overwritten. `cores` is a positive integer;
parallel execution remains limited to platforms supported by the underlying
S2 operator.

The required `recipe = "submitted_S2"` retains q0.975/type7 inclusive
trimming, arithmetic marginals, branch-mean then gene-mean division,
numeric-trait-plus-one parsimony with first ties, descendant-node annotation
and the strict-prefix 0.01 Welch screen. This is an explicit historical
recipe, not a recommended generic QC policy. Zero values, unavailable cells,
coordinate states and derived trim/ratio availability remain distinct.

Choose one validation recipe explicitly:

| `validation_recipe` | Saved computation contract |
| --- | --- |
| `nested_phase12B` | Foldwise screening retaining the phase12B historical arithmetic expression. |
| `nested_phase13` | Foldwise screening retaining the phase13 historical arithmetic expression. |
| `global_screen_foldwise` | Foldwise fitting using the globally screened feature order; does not repeat supervised feature selection inside each fold. |

The two nested recipes are retained separately rather than silently unified.
The global-screen recipe has a different selection boundary; its saved OOF
predictions do not establish nested feature-selection validation. Internal
branches may participate in screening, never terminal model fitting. All
groups, fold order and seeds are caller-supplied. Computational completion is
not a scientific decision, frozen Marine replay certification or release gate.

## Read and inspect a saved run

```r
result <- lt_read_binary_run(run)  # also accepts "/path/to/new_binary_run"
print(result)
summary(result)
names(result$tables)
result$diagnostics
result$configuration
```

The `lt_binary_result` contains `status`, `run_dir`, `completed`, `integrity`,
`tables`, `diagnostics` and `configuration`. `verify = TRUE` checks required
files and recorded hashes
for completed runs. These hashes detect changed bytes; they do not authenticate
scientific authority. `verify = FALSE` skips hash checks but still requires the
completed-run file set; it is explicitly unverified inspection. Failed runs keep their
partial outputs and recorded failure stage; diagnostics whose saved
dependencies are absent remain explicitly `not_run_or_incomplete` in product
views; the dependency ledger distinguishes `not_run` from
`incomplete_failed_stage`. Screen diagnostics distinguish `not_estimable`
from tested `not_selected` genes.

The run remains the existing flat output directory, including
`supplied_inputs.rds`, `configuration.rds`, `S2_baseline_GBI.rds`, branch and
screen TSVs, `validation.rds`, `OOF_predictions.tsv`, fold/seed ledgers,
`RUN_STATUS.yaml`, `run.rds` and `SHA256SUMS`. The reader exposes those saved
artifacts without moving or replacing them. It does not create fits or infer
successful dependencies from a partial directory.

Saved scientific keys remain literal character strings, including group keys
`001` and `1` (distinct), `NA`, `TRUE` and quoted labels. Numeric result columns
remain numeric. An internal branch may have the writer's missing terminal-taxon
sentinel, but an unexpected nonempty taxon on an internal branch is rejected,
not silently erased; a terminal literally named `NA` remains a taxon key.

`input_hashes.tsv` indexes the actual-input `supplied_inputs.rds` snapshot;
it is not an index of the original four input files. For file-backed inputs:

```r
supplied <- readRDS(file.path(result$run_dir, "supplied_inputs.rds"))
supplied$matrix$metadata$splitaligner_import[c("input_file", "input_file_sha256")]
supplied$trait$metadata$input_source
attr(supplied$terminal_groups, "input_source")
attr(supplied$folds, "input_source")
```

Directly supplied in-memory objects need not have original file metadata.
No original file identity is invented for them; the serialized supplied-input
snapshot is still recorded for the run.

## Report and plot without rerunning

```r
report <- lt_report_binary(result, "/path/to/new_binary_report")
p <- lt_plot_binary(result, view = "oof")
print(p)
lt_plot_binary(result, view = "roc")
lt_plot_binary(result, view = "trim")
```

The report directory must be new and outside the run. It contains `report.md`
and exported saved/diagnostic tables. It is not a replacement for the original
run bundle. D1-D4 organize only bounded, source-backed inspection:

| Layer | Meaning and boundary |
| --- | --- |
| D1 | Input/QC counts, keyed domains and saved state/reason summaries; no new QC decisions. |
| D2 | Saved trimming, denominator/ratio availability and domain summaries; no new threshold selection. |
| D3 | Saved branch annotation and screen support; no ASR or screen refit. |
| D4 | Saved OOF predictions, folds, feature use and warnings; no model fitting, tuning or regenerated folds. |

Each plot returns a `ggplot` built from saved data. `oof` displays saved
predictions, `roc` summarizes their ranking and `trim` displays saved trim
information. The OOF view draws every finite eligible prediction row with
one point shape, uses supplied groups for horizontal offsets and reports the
row/group counts. Group identities stay in the plot data for customization;
the number of groups does not require a shape palette or a large legend.
Missing dependencies cannot be supplied by refitting. Plotting
and reporting never modify the run, certify it or fill an unavailable result
with an inferred answer.
Every available plot from a failed run is explicitly titled and captioned
`PARTIAL / FAILED (stage: ...)`, even when its own early-stage data exist.
Missing plot dependencies still fail explicitly. The marker is read from
`RUN_STATUS.yaml`, never inferred from the user-selected directory name.

Interpret diagnostic denominators literally. D1 distinguishes all matrix
cells from available values; point masses and repeated exact extrema do not
establish their upstream causes. Inclusive trimming does not guarantee an
observed removal fraction of exactly 2.5%. Finite trimmed values whose GBI is
unavailable indicate denominator-domain loss, not a newly invented exclusion.
The historical Welch statistic is reference-minus-focal; an untested gene is
distinct from a tested but FDR-nonselected gene.

D4 preserves special-fold `NA` fields and saved null-intercept statuses.
Fold `n_removed` and `n_screened` count branches, not genes. Screen-selected
features, design features after all-missing/zero-variance drops, and nonzero
coefficients are distinct. Feature use across folds is descriptive, not
independent replication. Per-cell imputation and per-gene design-drop reason
ledgers were not persisted and are not manufactured by the reader. Pooled
OOF ROC/AUC uses finite eligible saved prediction rows; it does not establish
biological performance or independent fold replication.

## Installed small example

```r
example_dir <- system.file("examples", "binary_small", package = "LaTerra",
                           mustWork = TRUE)
input <- lt_import(file.path(example_dir, "exchange.rds"))
trait <- lt_read_trait(file.path(example_dir, "trait.tsv"), trait_name = "binary_trait")
run <- lt_run_binary(input, trait, file.path(example_dir, "groups.tsv"),
                     file.path(example_dir, "folds.tsv"), recipe = "submitted_S2",
                     validation_recipe = "nested_phase13",
                     output_dir = "/path/to/new_example_run")
```

The accompanying script can run directly from the installed package:

```r
source(system.file("examples", "run_binary_small.R", package = "LaTerra",
                   mustWork = TRUE))
```

It chooses new temporary run/report directories and prints their paths; no
command-line arguments are needed. The 24-tip, 32-gene exchange was generated from synthetic algebraic
inputs through the real producer; it is not biological evidence. Its purpose
is to demonstrate installation, explicit inputs, run/result/report/plot flow
and provenance. It is not a representative efficacy test, performance promise
or substitute for the separately authorized historical Marine replay.
