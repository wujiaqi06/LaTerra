\name{lt_diagnose_rate}
\alias{lt_diagnose_rate}
\alias{validate_lt_rate_diagnostic}
\alias{print.lt_rate_diagnostic}
\title{Trait-free threshold-free rate diagnostics}
\usage{
lt_diagnose_rate(rate, source_matrix = NULL, baseline_fit = NULL,
  measurement_fit_eligibility = NULL, extreme_n = 20L)

validate_lt_rate_diagnostic(x)

\method{print}{lt_rate_diagnostic}(x, ...)
}
\arguments{
  \item{rate}{A valid \code{lt_rate}.}
  \item{source_matrix}{Optional compatible \code{lt_matrix}; supplied only
  together with \code{baseline_fit}.}
  \item{baseline_fit}{Optional compatible \code{lt_c3_fit} or \code{lt_c2_fit}.}
  \item{measurement_fit_eligibility}{Optional exact Layer-2 declaration.}
  \item{extreme_n}{Nonnegative number of cells retained for each bounded
  descriptive extreme ranking.}
  \item{x}{An \code{lt_rate_diagnostic}.}
  \item{...}{Unused print arguments.}
}
\description{
Summarizes the native geometry of one materialized rate without inspecting a
trait or changing values, availability, coordinate state, or Layer-2.
}
\details{
The diagnostic is threshold-free and diagnose-only. It does not trim,
winsorize, classify normal-theory outliers, infer upstream causes, or select a
scientific winner. Optional source/baseline decomposition is permitted only
when both dependencies remain current; stale dependencies fail with
\code{STALE_RECOMPUTE_REQUIRED} and are never repaired automatically.

Validation checks scientific structure and internal semantic identity. It does
not require an ordinary user object to match a registered review-package or
golden-file SHA-256.

Current boundary-aware diagnosis treats \code{GBI > 0} as the sole finite
positive-log domain and exact admitted \code{Y = 0} with positive denominator
as a mathematical boundary. Historical finite-sentinel \code{logGBI}
diagnostics are descriptive/read-migrate-only evidence; they do not authorize
ordinary inference. Zero state is never inferred from equality to a sentinel,
and no upstream cause is inferred.
}
\value{
\code{lt_diagnose_rate()} returns an \code{lt_rate_diagnostic}; the validator
returns its input invisibly; the print method returns the object invisibly.
}
