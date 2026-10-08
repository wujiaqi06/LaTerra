% Generated for LT000A-T2B-CAL001.
\name{lt_calibrate_association}
\alias{lt_calibrate_association}
\alias{validate_lt_calibration}
\title{Calibrate a T2A association against external null worlds}
\usage{
lt_calibrate_association(
  rate,
  branch_state,
  association,
  domain,
  null_ensemble,
  calibration_spec = lt_calibration_spec(domain_mode = domain$domain_mode),
  calibration_id = paste0(association$association_id, "_",
    calibration_spec$statistic, "_", calibration_spec$alternative,
    "_calibration"),
  metadata = list()
)

validate_lt_calibration(x)
}
\arguments{
\item{rate}{Frozen \code{lt_rate} used by the observed association.}
\item{branch_state}{Observed \code{lt_branch_state}.}
\item{association}{Frozen observed \code{lt_association}.}
\item{domain}{Frozen \code{lt_association_domain}.}
\item{null_ensemble}{External \code{lt_branch_state_ensemble}.}
\item{calibration_spec}{An \code{lt_calibration_spec}.}
\item{calibration_id}{Stable result identifier.}
\item{metadata}{Named non-authoritative metadata list.}
\item{x}{Object to validate.}
}
\value{An \code{lt_calibration} retaining every gene; the validator returns
\code{x} invisibly.}
\description{
Recomputes the selected T2A statistic for each coherent external null world on
the frozen gene-specific domain. Empirical p-values use inclusive ties and the
plus-one correction. Non-estimable null replicates fail the affected gene
closed. Scores are streamed by replicate; null generation is out of scope.
}
\details{
The result validates calibration computation only. It does not certify the
scientific validity of the external null generator. Historical finite-sentinel
\code{lt_rate(representation = "logGBI")} objects fail closed in this ordinary
route. Migrate them with \code{lt_migrate_logGBI()} and calibrate ZERO_MASS and
POSITIVE_LOG separately. The \code{legacy_global_k10} view is
diagnostic/replay-only and is refused by ordinary calibration. This route never
authorizes a finite zero sentinel and never defines a joint boundary p-value.
}
