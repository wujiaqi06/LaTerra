\name{lt_rate}
\alias{lt_rate}
\alias{validate_lt_rate}
\alias{lt_rate_reason}
\alias{lt_rate_ledger}
\alias{print.lt_rate}
\title{Certified C3 production rate representations}
\usage{
lt_rate(x, c3_fit, rate_spec,
  measurement_fit_eligibility = c3_fit$domain$eligibility,
  metadata = list())

validate_lt_rate(x)

lt_rate_reason(x)

lt_rate_ledger(x, available_only = FALSE)

\method{print}{lt_rate}(x, ...)
}
\arguments{
  \item{x}{An exact source \code{lt_matrix} or an \code{lt_rate}.}
  \item{c3_fit}{A certified \code{lt_c3_fit}.}
  \item{rate_spec}{A production specification for \code{ADD_LT}, \code{GBI},
  or \code{logGBI_strict}. Historical \code{logGBI} specifications are
  read/migrate-only and cannot execute this constructor.}
  \item{measurement_fit_eligibility}{The exact Layer-2 declaration used by C3.}
  \item{metadata}{Optional named metadata. Production/sensitivity contract
  flags are authoritative; ordinary annotations and display labels are not.}
  \item{available_only}{Whether to materialize only available ledger rows.}
  \item{...}{Unused print arguments.}
}
\description{
Constructs deterministic rate views of one previously certified C3 fitted
mean. Input matrix, ordered axes, coordinate provenance, Layer-2, and baseline
identities are verified and no baseline is refit.
}
\details{
\code{ADD_LT = Y - mu} and \code{GBI = Y / mu}. Current log-relative analysis
is constructed explicitly with \code{lt_log_relative(GBI)}. It stores finite
natural-log values only where \code{GBI > 0}; an admitted observed \code{Y = 0}
with positive denominator is \code{EXACT_ZERO_BOUNDARY}, not a finite log
coordinate. \code{logGBI_strict} is the positive-only conditional companion.

Historical finite-sentinel \code{lt_rate(representation = "logGBI")} objects
remain readable and migratable but are forbidden in ordinary association or
calibration. Their global-k10 numerical view is available only through
\code{logGBI_encoded(..., mode = "legacy_global_k10")} only for diagnostic or
source-bound replay validation. No sentinel is the default public inferential
meaning. GBI and both log-relative forms are unavailable at zero baseline;
ADD_LT remains defined at a certified zero-boundary baseline. New GBI is
dimensionless and is not legacy GBI. Coordinate truth,
Layer-2 eligibility, baseline support, availability, and representation reason
remain separate state layers. QC trimming, traits, and RER are outside this
operator. Validation recomputes one scientific semantic identity from the
authoritative fields and embedded Layer-2 matrices/provenance. Mutable display
labels and explicitly non-authoritative annotations are excluded; exact frozen
snapshot comparison is available separately through \code{lt_audit_fingerprint()}.
The boundary-aware ZERO_MASS and POSITIVE_LOG components are separate
estimands; no joint p-value is defined by the V1 contract.
}
\value{
\code{lt_rate()} returns an \code{lt_rate}; \code{lt_rate_reason()} returns an
aligned character matrix; \code{lt_rate_ledger()} returns a keyed data frame;
the validator returns its input invisibly.
}
