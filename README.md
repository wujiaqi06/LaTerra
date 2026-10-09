# La Terra: A Matrix-First Framework for Branch-Wise Comparative Genomics

**La Terra is a matrix-first framework for branch-wise comparative genomics: it performs downstream analysis on matrices whose columns are scientifically identified branch coordinates, rather than arbitrary numerical features.**

**Branch coordinates are established and audited upstream by SplitAligner; La Terra begins once a coordinate-valid gene × branch matrix is available.**

Start with a matrix from [SplitAligner](https://github.com/wujiaqi06/SplitAligner)
and a trait you want to study. La Terra connects branch-length variation to
gene associations, predictive models and diagnostic plots, keeping each
measurement linked to the branch it represents.

La Terra is an R package under active development. Its current reusable
workflow supports binary traits and saves the inputs, analysis settings and
results together, so you can return to an analysis without reconstructing it
from a collection of scripts.

## Installation

Use R 4.2 or later and start a fresh R session:

```r
install.packages(c("remotes", "digest", "yaml", "ape", "castor",
                   "glmnet", "ggplot2"))
remotes::install_github("wujiaqi06/LaTerra",
                        dependencies = FALSE, upgrade = "never")
library(LaTerra)
```

SplitAlignerR is not required to install La Terra or to work with explicitly
prepared matrix objects. It is needed only for the optional exchange-file
interface described below. Development checks have been performed on macOS
arm64 with R 4.4.2; other environments still need validation.

## From SplitAligner to La Terra

The primary source of input is the **[Perl implementation of SplitAligner](https://github.com/wujiaqi06/SplitAligner)**
described in the SplitAligner paper. Keep its matrix together with the
reference species tree and the accompanying branch-coordinate and cell-state
information. These files tell La Terra what each column represents and which
measurements are available.

| Input | What it contributes |
| --- | --- |
| Gene × branch matrix | One row per gene and one column per reference branch, containing branch-length measurements. |
| Reference species tree | The topology and species names used for branch annotation and tree views. |
| Branch-coordinate information | The split defining each reference branch, so values follow branch identity rather than row order or a display label. |
| Cell-state information | Which cells contain observations and which are unavailable, retaining the upstream state labels. An observed zero remains a measurement. |

A one-step importer for Perl SplitAligner output files is in development.
For now, the runnable file-to-analysis route is the
[optional SplitAlignerR interface below](#try-the-optional-splitalignerr-interface);
`lt_import()` currently reads that exchange format.

## Add your trait and choose the validation groups

Your trait table can contain several traits. Select the one to analyze by
name, using species names that match the reference tree exactly.

```r
trait <- lt_read_trait("trait.tsv", trait_name = "my_trait")
```

For the current binary workflow, code the two classes as 0 and 1. Both
classes must be present. An optional 0.5 label excludes a species from
terminal prediction and its terminal branch from the gene screen; that
label is still supplied to ancestral-state reconstruction. Alongside the
trait, supply the groups to hold out and the order in which to test them:

| File | Columns |
| --- | --- |
| `trait.tsv` | `species` and a named numeric trait column, such as `my_trait`. |
| `groups.tsv` | `species` and `genus`, the group key used for validation. |
| `folds.tsv` | `fold_id`, `genus` and an integer `seed`. |

Each fold holds out one supplied group. Include every group in the fold
table, including groups containing only excluded species. Choosing these
groups is part of the study design: it determines which predictions will
count as held out.

## Follow the analysis through its results

An analysis recipe records the choices that turn the input matrix into
results. Inspect the available settings with `lt_binary_recipes()` and
choose the analysis and validation recipes explicitly.

The current workflow calculates **GBI** (gene–branch interaction) using a
trimmed arithmetic baseline, then connects ancestral-state annotation and
gene screening to sparse predictive models. The
[binary workflow guide](inst/workflow/BINARY.md#inspect-and-choose-recipes-then-run)
describes the calculations and the available validation designs. In the
nested recipes, screening is repeated within training folds; ancestral-state
annotation is computed once on the full tree. A separate global-screen recipe
reuses the full-data screen and therefore has a different validation meaning.

A run produces more than a gene list. Its saved results let you inspect:

- **Gene associations:** test results, effect directions and selected genes.
- **Predictive models:** sparse predictors, held-out predictions and ROC/AUC
  summaries, with the corresponding folds and seeds.
- **Branch annotation:** the trait states assigned to reference branches.
- **Diagnostics:** input availability, trimming, screen support, and
  fold-level feature use and warnings.

Read these together: they show both the result and the observations behind
it. A saved run can be reopened, reported and plotted without refitting:

```r
result <- lt_read_binary_run("/path/to/completed_run")
summary(result)
report <- lt_report_binary(result, "/path/to/new_report")
lt_plot_binary(result, view = "oof")  # also "roc" or "trim"
```

Use a new directory for each run and report, with an existing parent
directory. Saved runs include input snapshots, file hashes and software
information. If a run stops early, its partial outputs and failure stage
remain available for inspection.

## Try the optional SplitAlignerR interface

[SplitAlignerR](https://github.com/wujiaqi06/SplitAlignerR) is an optional
R interface that is still under development. Its exchange files bundle a
matrix with the species tree, branch coordinates, cell states and provenance.
The current La Terra importer supports **SplitAlignerR 0.1.0.9002** with
exchange schema **0.2.0-development**.

<details>
<summary>Install the matching development backend</summary>

Use the public [La Terra backend prerelease](https://github.com/wujiaqi06/SplitAlignerR/releases/tag/v0.1.0.9002-laterra-dev).
In a fresh R session, download and check the archive before installing it:

```r
backend <- tempfile(fileext = ".tar.gz")
download.file(
  paste0("https://github.com/wujiaqi06/SplitAlignerR/releases/download/",
         "v0.1.0.9002-laterra-dev/SplitAlignerR_0.1.0.9002.tar.gz"),
  backend, mode = "wb"
)
stopifnot(identical(
  digest::digest(file = backend, algo = "sha256"),
  "c4fcc16cce3aaf90189a686b57dd0eadba2d20fccfd246d42deeef353f83522f"
))
install.packages(backend, repos = NULL, type = "source")
stopifnot(as.character(packageVersion("SplitAlignerR")) == "0.1.0.9002")
```

Building this optional backend requires its source-compilation toolchain.

</details>

With the matching backend installed, try the bundled synthetic example:

```r
source(system.file("examples", "run_binary_small.R", package = "LaTerra"))
```

It contains 24 tips and 32 genes. The script runs the analysis, writes a
report and prints the output paths. It is a small software demonstration,
not a biological result or a performance benchmark.

For your own supported exchange file, trait and validation tables:

```r
settings <- subset(lt_binary_recipes(), validation_recipe == "nested_phase13")
stopifnot(nrow(settings) == 1L)
input <- lt_import("exchange.rds")
trait <- lt_read_trait("trait.tsv", trait_name = "my_trait")
run <- lt_run_binary(input, trait, "groups.tsv", "folds.tsv",
                     recipe = settings$recipe,
                     validation_recipe = settings$validation_recipe,
                     output_dir = "/path/to/new_run")
```

Keep the complete object returned by `lt_import()`; its tree, coordinates
and provenance belong with the matrix. See the
[input-table guide](inst/workflow/BINARY.md#supply-the-complete-keyed-inputs)
for details.

The [tree annotation and plotting guide](inst/workflow/SPLITALIGNER_PLOTS.md)
also shows how to display imported gene measurements on the reference tree
through this exchange interface.

## Related work and recommended citations

La Terra builds on the gene–branch interaction approach of Wu, Yonezawa and
Kishino (2017). SplitAligner provides the branch-coordinate foundation that
keeps measurements aligned across genes.

### Primary recommended citations

For analyses using La Terra with SplitAligner-derived matrices, please cite
both the foundational method and the upstream branch-mapping work:

1. **Wu, J., Yonezawa, T., and Kishino, H. (2017).**
   Rates of Molecular Evolution Suggest Natural History of Life History
   Traits and a Post-K-Pg Nocturnal Bottleneck of Placentals.
   *Current Biology*, **27**(19), 3025–3033.e5.
   [doi:10.1016/j.cub.2017.08.043](https://doi.org/10.1016/j.cub.2017.08.043).

2. **Wu, J. (2026).**
   SplitAligner: A Gene-Species Tree Reconciliation Framework Using
   Split-Based Branch Mapping.
   *bioRxiv* preprint.
   [doi:10.64898/2026.02.24.707838](https://doi.org/10.64898/2026.02.24.707838).
   [Software repository](https://github.com/wujiaqi06/SplitAligner).

### Related theoretical work

3. **Wu, J. (2026).**
   A Unique Graph-Theoretic Truth Table for Cross-Gene Branch Identity.
   *bioRxiv* preprint.
   [doi:10.64898/2026.07.22.740066](https://doi.org/10.64898/2026.07.22.740066).
   Cite this work when discussing the formal basis of cross-gene branch
   identity.

For reproducibility, also report the La Terra version or commit and the
upstream software versions used in your analysis.
