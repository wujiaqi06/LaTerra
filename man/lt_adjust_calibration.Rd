% Generated for LT000A-T2B-CAL001.
\name{lt_adjust_calibration}
\alias{lt_adjust_calibration}
\alias{validate_lt_calibration_adjustment}
\title{Adjust empirical p-values in one explicit calibration family}
\usage{
lt_adjust_calibration(calibration, method = NULL, metadata = list())

validate_lt_calibration_adjustment(x)
}
\arguments{
\item{calibration}{One \code{lt_calibration}.}
\item{method}{One of \code{none}, \code{BH}, \code{BY}, \code{Holm}, or
\code{Bonferroni}. If omitted, use the method declared by the calibration
specification.}
\item{metadata}{Named non-authoritative metadata list.}
\item{x}{Object to validate.}
}
\value{An \code{lt_calibration_adjustment}; the validator returns \code{x}
invisibly.}
\description{
Applies an explicitly requested adjustment to calibrated genes from exactly
one trait by rate representation by domain mode by statistic by alternative
family. Adjustment is never automatic. BH is available without a claim of
universal validity under arbitrary dependence; BY is the conservative
arbitrary-dependence FDR option. Holm and Bonferroni are FWER procedures.
}
