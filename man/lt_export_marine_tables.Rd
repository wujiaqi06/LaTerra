\name{lt_export_marine_tables}
\alias{lt_export_marine_tables}
\title{Export the frozen Marine reporting tables from a completed run}
\usage{lt_export_marine_tables(data_root, run_dir, output_dir)}
\arguments{
\item{data_root}{Root of the same frozen input bundle used by the M2 run.}
\item{run_dir}{Directory produced by \code{lt_run_marine_m2()}.}
\item{output_dir}{New directory outside the input and scientific-run roots.}
}
\description{A terminal reporting step using this-run outputs and explicitly
admitted reference metadata. It does not refit models, regenerate GBI, traits
or random worlds, or read archived result workbooks. Requires the optional R
packages \pkg{openxlsx}, \pkg{xml2}, and \pkg{zip}.}
\details{Preserves the original five workbook/eighteen worksheet logical schemas.
Historical references and later-frozen presentation rules are disclosed in
companion provenance. Computation is not independent TARGET001, scientific or
release certification. This fixed-profile identity boundary does not affect
ordinary LaTerra objects or other-data workflows. The terminal \code{run.rds}
retains all scientific ledgers with explicit source-directory scope, separately
identifies the reporting software/configuration and links all 22 required
targets to hashed evidence. No independent target acceptance is self-assigned.}
\value{Invisibly a list of workbook paths, reporting provenance, and an
\code{lt_run} combining source-scoped scientific and reporting provenance.}
\seealso{\code{\link{lt_run_marine_m2}}}
