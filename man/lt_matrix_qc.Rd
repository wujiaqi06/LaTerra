\name{lt_matrix_qc}
\alias{lt_qc_spec}
\alias{validate_lt_qc_spec}
\alias{lt_qc_matrix}
\alias{validate_lt_diagnostic}
\alias{lt_apply_qc}
\alias{print.lt_qc_spec}
\alias{print.lt_diagnostic}
\title{Generic matrix-QC diagnostics and optional Layer-2 admission}
\usage{
lt_qc_spec(profile = c("standard", "custom"),
  value_domain = c("auto", "signed", "nonnegative", "positive"),
  outlier_rules = NULL, top_k = c(1L, 5L, 10L),
  top_fractions = c(0.01, 0.05), tail_fractions = c(0.01, 0.05),
  point_mass_top_n = 10L)

validate_lt_qc_spec(x)

lt_qc_matrix(x, spec = NULL, profile = "standard")

validate_lt_diagnostic(x)

lt_apply_qc(x, diagnostic, exclude,
  measurement_fit_eligibility = NULL, decision_id, authority)

\method{print}{lt_qc_spec}(x, ...)

\method{print}{lt_diagnostic}(x, ...)
}
\arguments{
  \item{x}{An \code{lt_matrix}, \code{lt_qc_spec}, or
  \code{lt_diagnostic}, as appropriate.}
  \item{profile}{A named generic profile: \code{standard} or \code{custom}.}
  \item{value_domain}{The requested QC metric domain. \code{auto} resolves to
  signed when any available value is negative, otherwise nonnegative.}
  \item{outlier_rules}{A list of named rule lists. Methods are \code{none},
  \code{quantile}, \code{iqr}, \code{mad}, \code{absolute}, and
  \code{custom_mask}; scopes are whole matrix, within gene, and within branch.
  Rules are explicit optional parameters of the \code{custom} profile. A
  quantile rule must declare \code{tie_semantics} as \code{strict} or
  \code{inclusive}; \code{exact_rank} is parked and not implemented.}
  \item{top_k}{Positive integer cutoffs for leverage summaries.}
  \item{top_fractions}{Fractions of units for leverage summaries.}
  \item{tail_fractions}{Fractions for lower and upper tail summaries.}
  \item{point_mass_top_n}{Maximum repeated exact values to retain.}
  \item{spec}{An optional explicit \code{lt_qc_spec}.}
  \item{diagnostic}{An exact-matrix-aligned \code{lt_diagnostic}.}
  \item{exclude}{An explicit named list containing any of \code{flag_ids},
  \code{rule_ids}, or an aligned logical \code{mask}.}
  \item{measurement_fit_eligibility}{Optional existing Layer-2 declaration.}
  \item{decision_id}{Stable explicit application-decision identifier.}
  \item{authority}{Human or machine authority for applying selected exclusions.}
  \item{...}{Unused print arguments.}
}
\description{
\code{lt_qc_matrix()} produces deterministic, trait-independent cell, gene,
branch, and matrix summaries. Diagnosis never changes raw matrix values,
coordinate truth, payload reasons, axes, or Layer-2 eligibility. Exact zero is
valid. Signed concentration is explicitly based on absolute values; positive-
only log summaries report domain loss without a pseudocount. Global finite
minimum and maximum repetition are reported separately from the frequency-
ranked top-N point-mass table, with no upstream cause inferred.
}
\details{
The standard profile is diagnose-only and threshold-free: it computes
continuous distribution, tail, leverage, concentration, and availability
summaries without a binary outlier rule. Optional quantile, IQR, MAD, absolute,
or custom-mask rules require an explicit custom profile. A flag or warning is
not an exclusion. \code{lt_apply_qc()} requires an explicit decision. Selected
eligible cells become \code{ineligible / measurement_qc_excluded}; selected
unresolved cells may resolve to that state with their prior state recorded;
selected already-ineligible cells retain their exact reason. Selecting a
\code{not_applicable} cell fails closed. No trimming, winsorization,
clipping, smoothing, replacement, or other value repair is performed.

HOLD is restricted to a structural or mathematical contradiction established
by La Terra, such as finite negative values under an explicitly declared
nonnegative diagnostic domain. An unexplained finite extreme is a WARNING.

QC freezes the distinctions outlier is not error, long branch is not anomalous
cell, and trim is not diagnosis. Numeric trimming and a final trimming policy
are not implemented. The diagnostic object parks, but does not compute, future
\code{within_gene_percentile} and \code{within_branch_percentile} fields.
}
\value{
\code{lt_qc_spec()} returns an \code{lt_qc_spec}; \code{lt_qc_matrix()}
returns an \code{lt_diagnostic}; \code{lt_apply_qc()} returns an
\code{lt_measurement_eligibility}; validators return their input invisibly.
}
