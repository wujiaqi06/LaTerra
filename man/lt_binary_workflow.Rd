\name{lt_binary_workflow}
\alias{lt_import}
\alias{lt_read_trait}
\alias{lt_read_groups}
\alias{lt_read_folds}
\alias{lt_run_binary}
\alias{lt_binary_recipes}
\alias{lt_read_binary_run}
\alias{lt_report_binary}
\alias{lt_plot_binary}
\alias{print.lt_binary_result}
\alias{summary.lt_binary_result}
\title{Run an explicit binary-trait workflow and inspect saved results}
\usage{
lt_import(file)

lt_read_trait(file, trait_name, taxon_column = "species")

lt_read_groups(file)

lt_read_folds(file)

lt_run_binary(input, trait, terminal_groups, folds, recipe,
  validation_recipe, output_dir, cores = 1L, trait_name = NULL)

lt_binary_recipes()

lt_read_binary_run(run, verify = TRUE)

lt_report_binary(run, output_dir)

lt_plot_binary(run, view = c("oof", "roc", "trim"))

\method{print}{lt_binary_result}(x, ...)

\method{summary}{lt_binary_result}(object, ...)
}
\arguments{
\item{file}{An existing trusted exchange RDS for \code{lt_import()}, or a
header-bearing TSV for the trait, grouping or fold reader.}
\item{trait_name}{Explicit trait-column name. Required by \code{lt_read_trait()}
and when \code{lt_run_binary()} receives a trait TSV path.}
\item{taxon_column}{Explicit species-key column in the trait TSV. Keys are not
inferred, normalized or parsed from other fields.}
\item{input}{Complete accepted import returned by \code{lt_import()} or
\code{\link{lt_read_splitaligner_exchange}}, or its supported exchange RDS path.}
\item{trait}{Explicit binary \code{lt_trait}, or a TSV path together with
\code{trait_name}. Numeric 0 and 1 are required, optional 0.5 is explicitly
excluded, missing values are refused and both endpoint classes must occur.}
\item{terminal_groups}{Data frame or TSV path with explicit \code{species}
and \code{genus} keys. Species grouping is never inferred from spelling.}
\item{folds}{Complete ordered data frame or TSV path with \code{fold_id},
\code{genus} and integer-valued \code{seed}, including excluded-only genera.
No folds or seeds are generated.}
\item{recipe}{Required explicit \code{"submitted_S2"}. The wrapper does not
choose a scientific recipe for the caller.}
\item{validation_recipe}{Required explicit \code{"nested_phase12B"},
\code{"nested_phase13"} or \code{"global_screen_foldwise"}.}
\item{output_dir}{New output directory whose parent exists. A report directory
must also be outside the run. Existing destinations are not overwritten.}
\item{cores}{Positive integer worker count passed to the unchanged S2
grouped-validation operator; parallelism is platform dependent.}
\item{run}{Saved run directory or \code{lt_run}; reporting and plotting also
accept a previously read \code{lt_binary_result}.}
\item{verify}{Check recorded hashes for completed runs; required saved files
are checked even when \code{FALSE}. Disabling hash checks is explicitly
unverified inspection, not integrity evidence.}
\item{view}{Saved-result plot: \code{"oof"}, \code{"roc"} or \code{"trim"}.}
\item{x, object}{An \code{lt_binary_result} to print or summarize.}
\item{...}{Further method arguments.}
}
\description{
Provides an input-friendly route to the existing explicitly selected S2
binary-trait computation, and source-backed inspection of saved outputs.
La Terra is a matrix-first framework and evolving R-native platform. Marine/S2
is its first exact-replay historical example, not a future scope ceiling.
These wrappers do not infer traits, groupings, coordinates or validation design.
}
\details{
\code{lt_import()} retains the complete accepted SplitAlignerR matrix, tree,
coordinate, original-exchange and provenance bundle. It uses the same strict
boundary as \code{\link{lt_read_splitaligner_exchange}}; it does not guess
legacy B-label correspondence or perform upstream alignment, tree inference or
coordinate construction. The current import requires the matching SplitAlignerR
0.1.0.9002 public validator for schema 0.2.0-development. Check the package
actually loaded in an isolated library; 0.1.0 is not interchangeable.
S2 runs require the installed LaTerra runtime. Unsupported
\code{pkgload::load_all()} source sessions fail before computation or run
directory creation with source-package installation guidance; no alternate
runtime hash or developer scientific fallback is used.

The TSV readers preserve explicit keyed inputs. Trait reading requires numeric
0/0.5/1 coding, no missing values and both endpoints. Group and fold readers
require the declared columns. Full trait/tree/group domains and the supplied
fold design are checked by the run route. Extra measurement eligibility,
automatic trait recoding and fold generation are not introduced.

\code{lt_run_binary()} delegates computation to \code{\link{lt_run_s2}}.
The historical \code{submitted_S2} contract retains q0.975/type7 inclusive
trimming, arithmetic marginals, branch-mean then gene-mean division,
numeric-trait-plus-one parsimony with first ties, descendant-node annotation
and the strict-prefix 0.01 Welch screen. It is not a recommended generic
trimming or QC policy. The two nested validation recipes preserve their
distinct historical arithmetic expressions. \code{global_screen_foldwise}
uses globally screened features, not supervised feature selection repeated
inside each fold. Its saved OOF predictions do not establish nested
feature-selection validation. Internal branches can enter screens but never
terminal model fitting. See \code{\link{lt_s2_gloocv}}.

\code{lt_read_binary_run()} reads the existing flat run outputs without moving
or replacing them. Completed-run verification checks required files and
recorded hashes. Hashes detect changed bytes, not scientific authority.
Failed runs retain their partial files and failure stage; absent diagnostic
dependencies remain explicitly \code{not_run_or_incomplete}, never imputed successes.
Saved scientific identity columns, including \code{group}, retain literal
character keys such as distinct \code{001} and \code{1} and literal
\code{NA}/\code{TRUE}. Numeric result columns remain numeric. Unexpected
internal-branch terminal identities are rejected rather than erased.

\code{input_hashes.tsv} indexes the \code{supplied_inputs.rds} snapshot,
not original files. File-backed identities are retained inside that snapshot
at \code{matrix$metadata$splitaligner_import} (\code{input_file} and
\code{input_file_sha256}), \code{trait$metadata$input_source}, and
\code{input_source} attributes of \code{terminal_groups} and \code{folds}.
Direct in-memory inputs need not have original file identities; none are invented.

\code{lt_report_binary()} writes \code{report.md} and tables to a new external
report directory. Bounded diagnostics separate D1 input/QC counts and
state/reason summaries, D2 saved trimming and domains, D3 branch/screen
support, and D4 saved OOF predictions, folds, feature use and warnings.
Reporting does not refit, tune, regenerate folds or modify the run.
\code{lt_plot_binary()} returns plots from saved OOF, ROC or trimming data
and requires \pkg{ggplot2}; it does not recompute scientific stages.
The OOF view uses a single point shape and supplied groups for horizontal
offsets, draws all finite eligible prediction rows, and reports row/group
counts. Exact group keys remain in the editable plot data.
Available plots from failed runs explicitly display
\code{PARTIAL / FAILED (stage: ...)} in title and caption. Absent dependencies
still fail without refitting; completed-run plot semantics are unchanged.

The installed \code{examples/binary_small} directory contains a 24-tip,
32-gene synthetic algebraic example exported through the real producer, with
explicit trait/group/fold TSVs. \code{examples/run_binary_small.R} is its
runnable companion. This is a software demonstration, not biological evidence,
a performance promise, scientific acceptance or release certification.
Computational completion alone confers none of those conclusions.
}
\value{
\code{lt_import()} returns the complete accepted import list.
\code{lt_read_trait()} returns an explicit binary \code{lt_trait}; the group
and fold readers return keyed data frames. \code{lt_run_binary()} returns a
computed, uncertified \code{lt_run}. \code{lt_binary_recipes()} returns a
data frame describing the explicit computation and validation choices.
\code{lt_read_binary_run()} returns an \code{lt_binary_result} containing
\code{status}, \code{run_dir}, \code{completed}, \code{integrity},
\code{tables}, \code{diagnostics} and \code{configuration}.
The print and summary methods provide concise saved-run
inspection. \code{lt_report_binary()} returns report-output information;
\code{lt_plot_binary()} returns a \code{ggplot}.
}
\seealso{
\code{\link{lt_read_splitaligner_exchange}}, \code{\link{lt_run_s2}},
\code{\link{lt_s2_gloocv}}, \code{\link{lt_trait}},
\code{\link{lt_plot_splitaligner}}
}
\examples{
lt_binary_recipes()
\dontrun{
input <- lt_import("/path/to/exchange.rds")
trait <- lt_read_trait("/path/to/trait.tsv", trait_name = "my_trait")
run <- lt_run_binary(input, trait, "/path/to/groups.tsv", "/path/to/folds.tsv",
                     recipe = "submitted_S2", validation_recipe = "nested_phase13",
                     output_dir = "/path/to/new_binary_run")
result <- lt_read_binary_run(run)
summary(result)
lt_report_binary(result, "/path/to/new_binary_report")
lt_plot_binary(result, view = "oof")
}
}
