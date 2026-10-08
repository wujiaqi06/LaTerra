.lt_supported_object_classes <- function() {
  c(
    "lt_matrix", "lt_trait", "lt_domain", "lt_rate_spec", "lt_analysis_spec",
    "lt_run", "lt_c3_domain", "lt_c3_fit", "lt_qc_spec", "lt_diagnostic",
    "lt_measurement_eligibility", "lt_rate", "lt_c2_domain", "lt_c2_fit",
    "lt_rate_common_domain", "lt_rate_diagnostic", "lt_branch_state",
    "lt_association_domain", "lt_association", "lt_branch_state_ensemble",
    "lt_calibration_spec", "lt_calibration", "lt_calibration_adjustment",
    "lt_log_relative", "lt_logGBI_encoded", "lt_boundary_association",
    "lt_boundary_calibration", "lt_boundary_adjustments",
    "lt_legacy_encoded_association", "lt_legacy_encoded_calibration",
    "lt_legacy_encoded_adjustments"
  )
}

.lt_validate_supported_object <- function(x) {
  class_name <- intersect(class(x), .lt_supported_object_classes())
  if (length(class_name) != 1L) {
    .lt_abort("Object must be one supported La Terra S3 class.")
  }
  switch(
    class_name,
    lt_matrix = validate_lt_matrix(x),
    lt_trait = validate_lt_trait(x),
    lt_domain = validate_lt_domain(x),
    lt_rate_spec = validate_lt_rate_spec(x),
    lt_analysis_spec = validate_lt_analysis_spec(x),
    lt_run = validate_lt_run(x),
    lt_c3_domain = validate_lt_c3_domain(x),
    lt_c3_fit = validate_lt_c3_fit(x),
    lt_c2_domain = validate_lt_c2_domain(x),
    lt_c2_fit = validate_lt_c2_fit(x),
    lt_rate = validate_lt_rate(x),
    lt_rate_diagnostic = validate_lt_rate_diagnostic(x),
    lt_rate_common_domain = validate_lt_rate_common_domain(x),
    lt_branch_state = validate_lt_branch_state(x),
    lt_association_domain = validate_lt_association_domain(x),
    lt_association = validate_lt_association(x),
    lt_branch_state_ensemble = validate_lt_branch_state_ensemble(x),
    lt_calibration_spec = validate_lt_calibration_spec(x),
    lt_calibration = validate_lt_calibration(x),
    lt_calibration_adjustment = validate_lt_calibration_adjustment(x),
    lt_log_relative = validate_lt_log_relative(x),
    lt_logGBI_encoded = validate_lt_logGBI_encoded(x),
    lt_boundary_association = validate_lt_boundary_association(x),
    lt_boundary_calibration = validate_lt_boundary_calibration(x),
    lt_boundary_adjustments = .lt_validate_boundary_adjustments(x),
    lt_legacy_encoded_association =
      .lt_validate_legacy_encoded_association(x),
    lt_legacy_encoded_calibration =
      .lt_validate_legacy_encoded_calibration(x),
    lt_legacy_encoded_adjustments =
      .lt_validate_legacy_encoded_adjustments(x),
    lt_qc_spec = validate_lt_qc_spec(x),
    lt_diagnostic = validate_lt_diagnostic(x),
    lt_measurement_eligibility = {
      if (!inherits(x, "lt_measurement_eligibility")) {
        .lt_abort("Invalid measurement-fit eligibility object.")
      }
      if (!is.matrix(x$status) || !is.matrix(x$reason) ||
          !identical(dim(x$status), dim(x$reason)) ||
          !identical(dimnames(x$status), dimnames(x$reason)) ||
          !identical(x$status_sha256, .lt_hash(x$status)) ||
          !identical(x$reason_sha256, .lt_hash(x$reason))) {
        .lt_abort("Serialized measurement-fit eligibility structure/hash mismatch.")
      }
      invisible(x)
    }
  )
  invisible(x)
}

#' Serialize or restore a La Terra object
#'
#' Uses base R RDS serialization so classes, ordered ledgers, masks, and nested
#' provenance containers round-trip without scientific transformation.
#'
#' @param x Supported La Terra S3 object.
#' @param path Destination or source `.rds` path.
#' @param validate Whether to validate before writing or after reading.
#' @param overwrite Whether an existing destination may be replaced. The
#'   default is fail-closed; replacement is never silent.
#'
#' @return `lt_write_object()` returns the normalized path invisibly;
#'   `lt_read_object()` returns the restored object.
#' @export
lt_write_object <- function(x, path, validate = TRUE, overwrite = FALSE) {
  .lt_assert_scalar_character(path, "path")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    .lt_abort("`overwrite` must be one non-missing logical value.")
  }
  if (validate) .lt_validate_supported_object(x)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(path) && !overwrite) {
    .lt_abort(sprintf("Destination already exists; set overwrite=TRUE explicitly: %s", path))
  }
  temporary <- tempfile(pattern = ".laterra-object-", tmpdir = parent,
                        fileext = ".rds")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(x, temporary, version = 3)
  if (file.exists(path) && !overwrite) {
    .lt_abort(sprintf("Destination appeared during write; refusing replacement: %s", path))
  }
  if (!file.rename(temporary, path)) {
    .lt_abort(sprintf("Atomic object finalization failed: %s", path))
  }
  invisible(normalizePath(path, mustWork = TRUE))
}

#' @rdname lt_write_object
#' @export
lt_read_object <- function(path, validate = TRUE) {
  .lt_assert_scalar_character(path, "path")
  if (!file.exists(path)) {
    .lt_abort(sprintf("Serialized object does not exist: %s", path))
  }
  x <- readRDS(path)
  if (validate) .lt_validate_supported_object(x)
  x
}
