# Explicit post-calibration multiple-testing adjustment.

.lt_adjustment_payload <- function(x) {
  x[c(
    "calibration_identity", "adjustment_method", "family_identity",
    "n_tests", "gene_ids", "adjusted_p", "procedure_class",
    "provenance"
  )]
}

#' Adjust empirical calibration p-values in one explicit test family
#'
#' @param calibration One `lt_calibration` defining exactly one trait by rate
#'   representation by domain mode by statistic by alternative family.
#' @param method One of `none`, `BH`, `BY`, `Holm`, or `Bonferroni`. When
#'   omitted, the method declared in the calibration spec is used.
#' @param metadata Named non-authoritative metadata list.
#' @return An `lt_calibration_adjustment`. Non-calibrated genes remain present
#'   with `NA` adjusted p-values.
#' @export
lt_adjust_calibration <- function(calibration, method = NULL,
                                  metadata = list()) {
  validate_lt_calibration(calibration)
  if (is.null(method)) method <- calibration$calibration_spec$multiple_testing
  if (!is.character(method) || length(method) != 1L ||
      !method %in% .lt_adjustment_methods()) {
    .lt_abort("Unknown multiple-testing adjustment method.")
  }
  .lt_assert_named_list(metadata, "metadata")
  selected <- calibration$calibration_status == "CALIBRATED"
  adjusted <- rep(NA_real_, length(calibration$gene_ids))
  p_method <- switch(
    method, none = "none", BH = "BH", BY = "BY", Holm = "holm",
    Bonferroni = "bonferroni"
  )
  adjusted[selected] <- if (method == "none")
    calibration$p_empirical[selected] else
    stats::p.adjust(calibration$p_empirical[selected], method = p_method)
  family_payload <- list(
    trait_id = calibration$trait_id,
    rate_representation = calibration$rate_representation,
    domain_mode = calibration$domain_mode,
    statistic = calibration$statistic,
    alternative = calibration$alternative,
    calibrated_gene_ids = calibration$gene_ids[selected]
  )
  procedure_class <- switch(
    method,
    none = "NO_ADJUSTMENT",
    BH = "FDR_PROCEDURE_DEPENDENCE_NOT_UNIVERSALLY_GUARANTEED",
    BY = "FDR_PROCEDURE_ARBITRARY_DEPENDENCE_CONSERVATIVE",
    Holm = "FWER_PROCEDURE",
    Bonferroni = "FWER_PROCEDURE"
  )
  core <- list(
    calibration_identity = calibration$calibration_identity,
    adjustment_method = method,
    family_identity = .lt_hash(family_payload),
    n_tests = as.integer(sum(selected)),
    gene_ids = calibration$gene_ids,
    adjusted_p = adjusted,
    procedure_class = procedure_class,
    provenance = list(
      layer = "post_calibration_multiple_testing",
      automatic = FALSE,
      family_definition = paste(
        "one trait x one rate representation x one domain mode x",
        "one statistic x one alternative"
      ),
      gene_test_independence = "not_assumed",
      bh_arbitrary_dependence_validity = "not_claimed"
    )
  )
  out <- structure(c(core, list(
    adjustment_identity = .lt_hash(core), metadata = metadata
  )), class = "lt_calibration_adjustment")
  validate_lt_calibration_adjustment(out)
  out
}

#' Validate a calibration adjustment
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_calibration_adjustment <- function(x) {
  if (!inherits(x, "lt_calibration_adjustment")) {
    .lt_abort("`x` must inherit from `lt_calibration_adjustment`.")
  }
  .lt_assert_sha256(x$calibration_identity, "calibration_identity")
  if (!x$adjustment_method %in% .lt_adjustment_methods()) {
    .lt_abort("Unknown calibration adjustment method.")
  }
  .lt_assert_sha256(x$family_identity, "family_identity")
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  if (!is.numeric(x$adjusted_p) ||
      length(x$adjusted_p) != length(x$gene_ids) ||
      any(x$adjusted_p < 0 | x$adjusted_p > 1, na.rm = TRUE) ||
      !is.integer(x$n_tests) || length(x$n_tests) != 1L ||
      x$n_tests != sum(!is.na(x$adjusted_p))) {
    .lt_abort("Adjusted-p ledger is malformed.")
  }
  .lt_assert_scalar_character(x$procedure_class, "procedure_class")
  .lt_assert_named_list(x$provenance, "provenance")
  .lt_assert_named_list(x$metadata, "metadata")
  .lt_assert_sha256(x$adjustment_identity, "adjustment_identity")
  if (!identical(x$adjustment_identity,
                 .lt_hash(.lt_adjustment_payload(x)))) {
    .lt_abort("Calibration-adjustment semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

print.lt_calibration_adjustment <- function(x, ...) {
  cat("<lt_calibration_adjustment>", x$adjustment_method, "\n")
  cat("  tests:", x$n_tests,
      " family:", substr(x$family_identity, 1L, 12L), "...\n")
  invisible(x)
}
