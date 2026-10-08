% Generated for LT000A-T2B-CAL001.
\name{lt_calibration_spec}
\alias{lt_calibration_spec}
\alias{validate_lt_calibration_spec}
\title{Declare or validate an empirical calibration specification}
\usage{
lt_calibration_spec(
  statistic = "beta",
  alternative = c("two_sided", "greater", "less"),
  multiple_testing = c("none", "BH", "BY", "Holm", "Bonferroni"),
  domain_mode = c("NATIVE_ASSOCIATION_DOMAIN", "COMMON_ASSOCIATION_DOMAIN"),
  nonestimable_null_policy = "fail_closed",
  metadata = list()
)

validate_lt_calibration_spec(x)
}
\arguments{
\item{statistic}{Explicit T2A statistic.}
\item{alternative}{\code{greater}, \code{less}, or \code{two_sided}.}
\item{multiple_testing}{Requested later adjustment method. No adjustment is
applied during calibration.}
\item{domain_mode}{Native or common frozen association-domain mode.}
\item{nonestimable_null_policy}{Only \code{fail_closed} is supported.}
\item{metadata}{Named non-authoritative metadata list.}
\item{x}{Object to validate.}
}
\value{An \code{lt_calibration_spec}; the validator returns \code{x}
invisibly.}
\description{
Declares the association statistic, inclusive-tail alternative, domain mode,
and fail-closed null-replicate policy. The association-statistic neutral point
is zero and is distinct from numerical zero in any rate representation.
}
