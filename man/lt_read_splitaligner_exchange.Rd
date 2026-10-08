\name{lt_read_splitaligner_exchange}
\alias{lt_read_splitaligner_exchange}
\title{Read a versioned SplitAlignerR exchange into explicit S2 inputs}
\usage{
lt_read_splitaligner_exchange(file, allow_synthetic_migration = FALSE)
}
\arguments{
\item{file}{Existing exchange RDS from a trusted producer. Stored historical source paths are not reopened.}
\item{allow_synthetic_migration}{Explicitly permit declared synthetic label-migration fixtures. Never authorizes real historical migration.}
}
\value{A list with \code{matrix} (\code{lt_matrix}), complete \code{tree}, S2 \code{coordinates}, unchanged original \code{exchange}, and import \code{provenance}. Save the complete bundle to retain all upstream evidence.}
\description{Experimental single-fixed primitive branch-length reader. All six original matrices must have identical declared keyed axes before construction. An explicit transpose applies to all six. Scientific keys, full phylo, source reasons and composite ledgers are preserved; composites are not ordinary edges. Reading does not compute coordinates, baseline, ASR or models.}
\details{Supports schema 0.2.0-development with the SplitAlignerR 0.1.0.9002 public validator. Canonical hashes and native label/reference agreement are mandatory upstream checks, not duplicate LT graph algorithms. Old schema 0.1.0-development (INTEROP001-H01) and unknown versions are refused without resealing, automatic conversion, locale changes or bypass. Scoped H01 evidence covers R4.4.2 C/UTF-8, not actual execution under another R version or release certification. The explicit source-bound legacy profile \code{marine302-user-bound-20260908-v1} is admitted with its complete migration history, schema 0.2.0-development and SAR-C14N-1-SHA256 contract. The upstream public validator checks embedded original source bytes, exact tree instance and both label ledgers. Unknown profiles and bare ambiguous historical matrices are not admitted; no crosswalk is guessed or applied twice. This supports only the named source-bound correspondence, not historical generation-chain certification. This strict import boundary does not make golden SHA hashes a requirement on ordinary user matrices. Traits, grouping, folds and scientific recipes remain explicit arguments to \code{\link{lt_run_s2}}. Complete original source description is retained even though the La Terra payload origin is normalized to \code{imported}. Use \code{\link{lt_plot_splitaligner}} to plot without matching branch labels or node positions yourself.}
