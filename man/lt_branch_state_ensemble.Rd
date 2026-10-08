% Generated for LT000A-T2B-CAL001.
\name{lt_branch_state_ensemble}
\alias{lt_branch_state_ensemble}
\alias{validate_lt_branch_state_ensemble}
\title{Construct or validate an external branch-state null ensemble}
\usage{
lt_branch_state_ensemble(
  ensemble_id,
  trait_id,
  branch_ids,
  replicate_ids,
  values,
  type = c("binary", "continuous"),
  coding = list(),
  availability,
  null_hypothesis,
  generator_provenance = list(),
  observed_realization_included = FALSE,
  metadata = list()
)

validate_lt_branch_state_ensemble(x)
}
\arguments{
\item{ensemble_id}{Stable ensemble identifier.}
\item{trait_id}{Trait concept represented by the null worlds.}
\item{branch_ids}{Ordered unique scientific branch keys.}
\item{replicate_ids}{Ordered unique null-world identifiers.}
\item{values}{Atomic branch-by-replicate matrix.}
\item{type}{Either \code{"binary"} or \code{"continuous"}.}
\item{coding}{Named coding declaration. Binary ensembles require explicit
\code{reference_level} and \code{focal_level}.}
\item{availability}{Branch-by-replicate availability states.}
\item{null_hypothesis}{Scientific null declaration supplied by the external
authority.}
\item{generator_provenance}{Named external-generator provenance list.}
\item{observed_realization_included}{Must be \code{FALSE}.}
\item{metadata}{Named non-authoritative metadata list.}
\item{x}{Object to validate.}
}
\value{An \code{lt_branch_state_ensemble}; the validator returns \code{x}
invisibly.}
\description{
Records externally generated coherent branch-state null worlds. Each column
is one complete realization shared across all genes. The constructor does not
generate, shuffle, permute, or simulate branch states.
}
