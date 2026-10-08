\name{lt_c2_sensitivity}
\alias{lt_c2_domain}
\alias{validate_lt_c2_domain}
\alias{lt_c2_control}
\alias{lt_c2_fit}
\alias{validate_lt_c2_fit}
\alias{lt_c2_sensitivity}
\alias{print.lt_c2_domain}
\alias{print.lt_c2_fit}
\title{Required C2 positive-cell log-OLS sensitivity}
\usage{
lt_c2_domain(x, measurement_fit_eligibility = NULL)

validate_lt_c2_domain(x)

lt_c2_control(tolerance = 1e-10, max_iterations = 1000L,
  normal_equation_tolerance = 4e-10, replay = TRUE)

lt_c2_fit(domain, control = lt_c2_control(), run_id = "c2_sensitivity")

validate_lt_c2_fit(x)

lt_c2_sensitivity(x, c2_fit,
  representation = c("ADD_C2", "GBI_C2", "logGBI_C2"),
  measurement_fit_eligibility = c2_fit$domain$eligibility,
  metadata = list())

\method{print}{lt_c2_domain}(x, ...)

\method{print}{lt_c2_fit}(x, ...)
}
\arguments{
  \item{x}{An exact input \code{lt_matrix}, \code{lt_c2_domain}, or
  \code{lt_c2_fit}.}
  \item{measurement_fit_eligibility}{Exact Layer-2 declaration.}
  \item{tolerance}{Maximum alternating-effect change.}
  \item{max_iterations}{Maximum deterministic sweeps.}
  \item{normal_equation_tolerance}{Maximum absolute mean score residual.}
  \item{replay}{Whether to require deterministic solver replay.}
  \item{domain}{An \code{lt_c2_domain}.}
  \item{control}{A normalized \code{lt_c2_control()} list.}
  \item{run_id}{Stable sensitivity-run ID.}
  \item{c2_fit}{A certified \code{lt_c2_fit}.}
  \item{representation}{A sensitivity-only algebraic view.}
  \item{metadata}{Optional named metadata.}
  \item{...}{Unused print arguments.}
}
\description{
Fits the required C2 positive-cell two-way log-OLS alternative estimand and
derives explicitly sensitivity-labelled rate views.
}
\details{
The fit domain contains only observed, Layer-2 eligible, finite positive cells.
Exact zeros are excluded with no pseudocount. Components are solved and gauged
independently; factor values are gauge dependent but fitted means are invariant.
Evaluation at a zero cell requires supported gene and branch vertices in the
same component. Cross-component products are never formed. C2 cannot be selected
as the production baseline. Fit validation recomputes effects, factors,
baseline hashes, component-certificate identities, and the top-level semantic
fit identity without rerunning the solver.
}
\value{
Constructors return \code{lt_c2_domain}, control list, \code{lt_c2_fit}, or a
sensitivity-labelled \code{lt_rate}; validators return their input invisibly.
}
