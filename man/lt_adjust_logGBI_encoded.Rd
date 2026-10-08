\name{lt_adjust_logGBI_encoded}
\alias{lt_adjust_logGBI_encoded}
\title{Adjust explicit historical encoded calibration families}
\usage{
lt_adjust_logGBI_encoded(calibration, methods = c("BH", "BY"))
}
\arguments{
\item{calibration}{An \code{lt_legacy_encoded_calibration}.}
\item{methods}{Separate adjustment methods, normally \code{BH} and \code{BY}.}
}
\value{This function always stops with a classed refusal.}
\description{
This retired adjustment entry point always fails closed. Historical encoded
objects are diagnostic/replay-only and cannot enter La Terra V1 inference.
}
