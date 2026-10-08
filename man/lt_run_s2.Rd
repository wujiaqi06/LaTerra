\name{lt_run_s2}
\alias{lt_run_s2}
\title{Run an explicitly selected S2 workflow on user-owned data}
\usage{
lt_run_s2(x, tree, coordinates, trait, terminal_groups, folds,
  recipe, validation_recipe, output_dir, cores = 1L)
}
\arguments{
\item{x}{Coordinate-valid nonnegative branch-length \code{lt_matrix}.}
\item{tree}{Supplied compatible \code{phylo} tree; never inferred here.}
\item{coordinates}{Existing ledger: branch_id, canonical_split_key, side_A_taxa,
side_B_taxa (semicolon-delimited), branch_type, terminal_taxon.}
\item{trait}{Binary \code{lt_trait} explicitly encoded 0/1 with optional excluded 0.5.}
\item{terminal_groups}{Data frame with species and explicit genus grouping keys.}
\item{folds}{Full ordered fold_id/genus/seed ledger, including excluded-only groups.}
\item{recipe}{Required literal \code{"submitted_S2"}; no automatically selected recipe.}
\item{validation_recipe}{Required \code{"nested_phase12B"}, \code{"nested_phase13"}, or \code{"global_screen_foldwise"}.}
\item{output_dir}{New output directory; no overwrite.}
\item{cores}{Positive integer grouped-validation worker count.}
}
\description{
Validates user data structurally, then executes the explicitly selected historical
S2 baseline, GBI, deterministic annotation, single-gene screen and grouped model
workflow. No Marine hashes, species, matrix dimensions or annotation tables are
required. The caller's matrix, coordinate states, display labels and metadata are
not modified. Actual-input hashes are recorded after validation, not required as
a pre-existing frozen identity.
}
\details{
Requires an installed LaTerra source package with its own
\code{R/LaTerra.rdb} runtime. The loaded namespace is checked before
annotation, computation or output-directory creation. Unsupported
\code{pkgload::load_all()} sessions receive installation guidance; the
installed-runtime hash is not borrowed from another package copy or invented.

This is an explicit historical recipe, not a recommended generic trimming policy.
It retains q0.975 type7 inclusive trimming, arithmetic marginals, branch-mean then
gene-mean division, numeric-trait-plus-one castor parsimony with first ties,
descendant-node annotation and the strict-prefix 0.01 Welch screen. The original
screen columns named marine/non_marine mean focal/reference counts in this recipe;
the supplied trait determines their actual domain. Internal branches can enter
screening but never terminal model fitting. Grouping, folds and seeds are explicit.

No new coordinate construction, baseline choice, trait discretization, null
generator, model fallback or Marine-specific legacy direction exception is
introduced. Continuous traits and unsupported input structures fail explicitly.
Output completion is computational status only, not scientific certification.
The input object and all state/reason layers are preserved separately from the
derived trimming and ratio availability. This entry does not implement general
Layer-2 QC decisions or infer missing measurement eligibility.
}
\value{A computed \code{lt_run}; source inputs, derived baseline, annotations,
screen, per-fold fits/predictions and provenance are retained in the output directory.}
\seealso{\code{\link{lt_matrix}}, \code{\link{lt_trait}}, \code{\link{lt_s2_gloocv}}}
