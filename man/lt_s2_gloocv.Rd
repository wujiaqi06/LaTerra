\name{lt_s2_gloocv}
\alias{lt_s2_gloocv}
\title{Execute an explicitly selected S2 grouped validation recipe}
\usage{
lt_s2_gloocv(gbi, branch_states, terminals, folds, recipe,
  global_features = NULL, progress = NULL, cores = 1L)
}
\arguments{
\item{gbi}{Nonnegative numeric gene-by-branch GBI matrix with unique scientific
keys in its dimnames; NA denotes unavailable payload.}
\item{branch_states}{Named numeric branch states: 0, 0.5 (excluded), or 1.}
\item{terminals}{Data frame with unique character species and branch keys,
character genus and numeric response. Row order defines sample order.}
\item{folds}{Complete ordered \code{fold_id}/genus/seed data frame, including genera
with no eligible test terminal. Seeds must be integer-valued.}
\item{recipe}{Required choice of \code{nested_phase12B}, \code{nested_phase13}
or \code{global_screen_foldwise}. No automatic scientific default.}
\item{global_features}{Ordered global feature keys; required for the global
screen recipe and optional reference keys for nested reporting.}
\item{progress}{Optional callback taking completed count, total and current
result. Its return value does not replace the computed result.}
\item{cores}{Positive integer worker count. Ordered fork batches are used for
more than one worker on supported platforms; unavailable parallelism fails explicitly.}
}
\value{A list of per-fold results and input/specification provenance snapshots.
Completion does not certify frozen replay, the method or a release.}
\description{Reusable S2 operators for explicitly declared data, binary traits
and folds. There are no Marine hash, dimension, species or gene requirements.
No baseline, ASR, world generator or automatic trait inference is performed.}
\details{The nested recipes preserve their historical floating-point expressions,
remove only held-out terminal branches from the screen, retain internal screening
branches, and select a strict sorted-p prefix at 0.01 (not a general BH step-up).
The global recipe uses the supplied globally screened feature order; it does not
pretend to perform nested feature selection. Models use eligible training
terminals only, training mean imputation, sample-SD scaling, binomial glmnet
with standardization disabled, inner foldid equal to the training sample sequence,
deviance CV and lambda.min. Declared null fallbacks have explicit status labels;
model errors are not silently converted into null fits. Caller RNG state is
preserved. Provenance records actual input hashes but does not compare ordinary
user inputs to audit snapshot identities.}
