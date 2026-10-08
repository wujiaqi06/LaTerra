\name{lt_calibrate_logGBI_encoded}
\alias{lt_calibrate_logGBI_encoded}
\title{Calibrate the explicit historical finite-encoded view}
\usage{
lt_calibrate_logGBI_encoded(
  encoded,
  branch_state,
  association,
  null_ensembles,
  encoded_mode = NULL,
  acknowledge_encoded_zero = FALSE,
  alternative = c("two_sided", "greater", "less"),
  calibration_id = paste0(association$association_id, "_calibration")
)
}
\arguments{
\item{encoded}{An \code{lt_logGBI_encoded} view.}
\item{branch_state}{Observed binary \code{lt_branch_state}.}
\item{association}{Matching \code{lt_legacy_encoded_association}.}
\item{null_ensembles}{Named frozen external null ensembles.}
\item{encoded_mode}{Deprecated compatibility argument; it cannot authorize
inference.}
\item{acknowledge_encoded_zero}{Deprecated compatibility argument. Its value
has no effect on the mandatory refusal.}
\item{alternative}{Empirical alternative.}
\item{calibration_id}{Stable bundle identifier.}
}
\value{This function always stops with a classed refusal.}
\description{
This retired entry point always fails closed. Source-bound replay identity
validation is available through \code{validate_legacy_logGBI_replay()}, which
computes no association statistic.
}
