\name{lt_c3_engine}
\alias{lt_c3_existence_states}
\alias{lt_measurement_eligibility}
\alias{validate_lt_measurement_eligibility}
\alias{lt_c3_domain}
\alias{validate_lt_c3_domain}
\alias{lt_c3_control}
\alias{lt_c3_support}
\alias{lt_c3_fit}
\alias{validate_lt_c3_certificate}
\alias{validate_lt_c3_fit}
\alias{print.lt_c3_domain}
\alias{print.lt_c3_fit}
\title{Certified C3 baseline engine}
\usage{
lt_c3_existence_states()

lt_measurement_eligibility(x, status = NULL, reason = NULL,
  provenance = list())

validate_lt_measurement_eligibility(x, matrix_object)

lt_c3_domain(x, measurement_fit_eligibility = NULL,
  diagnostic_warnings = character())

validate_lt_c3_domain(x)

lt_c3_control(fit_rtol = 1e-10, fit_atol_scale = 1e-12,
  max_iter = 10000L, certificate_max_edges = 2e7, replay = TRUE)

lt_c3_support(domain, control = lt_c3_control())

lt_c3_fit(domain, control = lt_c3_control(), run_id = "c3_run")

validate_lt_c3_certificate(certificate)

validate_lt_c3_fit(x)

\method{print}{lt_c3_domain}(x, ...)

\method{print}{lt_c3_fit}(x, ...)
}
\arguments{
  \item{x}{An appropriate La Terra matrix or C3 object.}
  \item{status}{Optional layer-2 eligibility-state matrix.}
  \item{reason}{Optional layer-2 reason matrix.}
  \item{provenance}{Named layer-2 declaration provenance.}
  \item{measurement_fit_eligibility}{Optional explicit layer-2 object.}
  \item{matrix_object}{The aligned \code{lt_matrix} providing coordinate truth.}
  \item{diagnostic_warnings}{Audit-only warnings that cannot mutate eligibility.}
  \item{fit_rtol}{Relative fitted-margin tolerance.}
  \item{fit_atol_scale}{Absolute fitted-margin tolerance multiplier.}
  \item{max_iter}{Maximum deterministic IPF sweeps.}
  \item{certificate_max_edges}{Fail-closed certificate resource ceiling.}
  \item{replay}{Whether to rerun and hash-compare the numerical fit.}
  \item{domain}{A declared C3 fit domain.}
  \item{control}{A deterministic C3 control list.}
  \item{run_id}{Stable run identifier.}
  \item{certificate}{A machine-readable C3 fit certificate.}
  \item{...}{Unused print arguments.}
}
\description{
Constructs the declared nonnegative fit domain, bipartite support graph,
zero-margin and zero-total boundary ledgers, exact residual-SCC minimal-face
certificate, and deterministic generalized-KL/Bregman margin working mean.
No off-domain outer product is returned. This engine is not a literal Poisson
probability model and does not issue scientific PASS decisions.
Layer 2 is \code{not_applicable} with reason
\code{coordinate_not_observed} outside observed coordinate truth.
}
\value{The state vocabulary, an eligibility/domain/support object, a control
list, a certified \code{lt_c3_fit}, or the validated input invisibly.}
