#' Construct a La Terra taxon domain
#'
#' @param domain_id Stable domain identifier.
#' @param taxon_ids Ordered unique taxon identifiers.
#' @param eligible Logical eligibility mask aligned to `taxon_ids`.
#' @param exclusion_reason Character exclusion ledger aligned to `taxon_ids`.
#' @param group Optional grouping identifiers, such as genus.
#' @param metadata Optional named metadata list.
#'
#' @return An object of class `lt_domain`.
#' @export
lt_domain <- function(domain_id,
                      taxon_ids,
                      eligible = rep(TRUE, length(taxon_ids)),
                      exclusion_reason = rep(NA_character_, length(taxon_ids)),
                      group = rep(NA_character_, length(taxon_ids)),
                      metadata = list()) {
  x <- structure(
    list(
      domain_id = domain_id,
      taxon_ids = taxon_ids,
      eligible = eligible,
      exclusion_reason = exclusion_reason,
      group = group,
      metadata = metadata
    ),
    class = "lt_domain"
  )
  validate_lt_domain(x)
  x
}

#' Validate a La Terra taxon domain
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_domain <- function(x) {
  if (!inherits(x, "lt_domain")) {
    .lt_abort("`x` must inherit from `lt_domain`.")
  }
  .lt_assert_scalar_character(x$domain_id, "domain_id")
  .lt_assert_unique_ids(x$taxon_ids, "taxon_ids")
  n <- length(x$taxon_ids)
  if (!is.logical(x$eligible) || length(x$eligible) != n || anyNA(x$eligible)) {
    .lt_abort("`eligible` must be a non-missing logical vector aligned to `taxon_ids`.")
  }
  if (!is.character(x$exclusion_reason) || length(x$exclusion_reason) != n) {
    .lt_abort("`exclusion_reason` must be a character vector aligned to `taxon_ids`.")
  }
  explicit_reason <- !is.na(x$exclusion_reason) &
    nzchar(trimws(x$exclusion_reason))
  if (any(!x$eligible & !explicit_reason)) {
    .lt_abort("Every ineligible taxon requires a non-empty `exclusion_reason`.")
  }
  if (any(x$eligible & !is.na(x$exclusion_reason))) {
    .lt_abort("Eligible taxa must have an NA `exclusion_reason`.")
  }
  if (!is.character(x$group) || length(x$group) != n) {
    .lt_abort("`group` must be a character vector aligned to `taxon_ids`.")
  }
  .lt_assert_named_list(x$metadata, "metadata")
  invisible(x)
}
