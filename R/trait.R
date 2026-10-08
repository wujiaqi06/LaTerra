#' Construct a La Terra trait
#'
#' @param trait_id Stable trait identifier.
#' @param taxon_ids Ordered unique taxon identifiers.
#' @param values Atomic trait values aligned to `taxon_ids`.
#' @param type One of `binary`, `categorical`, `ordinal`, or `continuous`.
#' @param coding Named list describing positive, negative, intermediate, or
#'   excluded values without interpreting them.
#' @param metadata Optional named metadata list.
#'
#' @return An object of class `lt_trait`.
#' @export
lt_trait <- function(trait_id,
                     taxon_ids,
                     values,
                     type = c("binary", "categorical", "ordinal", "continuous"),
                     coding = list(),
                     metadata = list()) {
  type <- match.arg(type)
  x <- structure(
    list(
      trait_id = trait_id,
      taxon_ids = taxon_ids,
      values = values,
      type = type,
      coding = coding,
      metadata = metadata
    ),
    class = "lt_trait"
  )
  validate_lt_trait(x)
  x
}

#' Validate a La Terra trait
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_trait <- function(x) {
  if (!inherits(x, "lt_trait")) {
    .lt_abort("`x` must inherit from `lt_trait`.")
  }
  .lt_assert_scalar_character(x$trait_id, "trait_id")
  .lt_assert_unique_ids(x$taxon_ids, "taxon_ids")
  if (!is.atomic(x$values) || is.list(x$values)) {
    .lt_abort("`values` must be an atomic vector.")
  }
  if (length(x$values) != length(x$taxon_ids)) {
    .lt_abort("Trait values and taxon IDs must have identical lengths.")
  }
  if (!x$type %in% c("binary", "categorical", "ordinal", "continuous")) {
    .lt_abort("Unknown trait type.")
  }
  .lt_assert_named_list(x$coding, "coding")
  .lt_assert_named_list(x$metadata, "metadata")
  invisible(x)
}

#' Validate trait and domain taxon compatibility
#'
#' @param trait An `lt_trait`.
#' @param domain An `lt_domain`.
#' @param require_exact If `TRUE`, both taxon sets must be identical. If
#'   `FALSE`, every trait taxon must be present in the domain.
#'
#' @return `TRUE`, invisibly, or an informative mismatch error.
#' @export
validate_lt_trait_domain <- function(trait, domain, require_exact = TRUE) {
  validate_lt_trait(trait)
  validate_lt_domain(domain)
  if (!is.logical(require_exact) || length(require_exact) != 1L || is.na(require_exact)) {
    .lt_abort("`require_exact` must be TRUE or FALSE.")
  }

  missing_from_domain <- setdiff(trait$taxon_ids, domain$taxon_ids)
  missing_from_trait <- if (require_exact) {
    setdiff(domain$taxon_ids, trait$taxon_ids)
  } else {
    character()
  }

  if (length(missing_from_domain) > 0L || length(missing_from_trait) > 0L) {
    detail <- c(
      if (length(missing_from_domain) > 0L) {
        paste0("trait-only: ", paste(utils::head(missing_from_domain, 5L), collapse = ", "))
      },
      if (length(missing_from_trait) > 0L) {
        paste0("domain-only: ", paste(utils::head(missing_from_trait, 5L), collapse = ", "))
      }
    )
    .lt_abort(paste0("Trait/taxon mismatch detected (", paste(detail, collapse = "; "), ")."))
  }
  invisible(TRUE)
}
