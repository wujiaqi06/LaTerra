# Frozen-snapshot audit fingerprints and dependency-staleness checks.

#' Fingerprint an exact frozen object snapshot
#'
#' Audit fingerprints answer only whether two snapshots are byte-semantically
#' the same serialized La Terra object. They are not general user-facing
#' validity certificates: editing display labels or non-authoritative metadata
#' may change a snapshot fingerprint while leaving the object scientifically
#' usable.
#'
#' @param x A supported La Terra object.
#' @return An `lt_audit_fingerprint` containing class, algorithm, and SHA-256.
#' @export
lt_audit_fingerprint <- function(x) {
  .lt_validate_supported_object(x)
  structure(list(
    object_class = class(x),
    algorithm = "sha256_R_serialization_v3",
    snapshot_sha256 = .lt_hash(list(class = class(x), object = unclass(x)))
  ), class = "lt_audit_fingerprint")
}

.lt_validate_audit_fingerprint <- function(fingerprint) {
  if (!inherits(fingerprint, "lt_audit_fingerprint") ||
      !is.character(fingerprint$object_class) ||
      !identical(fingerprint$algorithm, "sha256_R_serialization_v3")) {
    .lt_abort("`fingerprint` must be an lt_audit_fingerprint.")
  }
  .lt_assert_sha256(fingerprint$snapshot_sha256, "snapshot_sha256")
  invisible(fingerprint)
}

#' Test or require exact frozen-snapshot identity
#'
#' A mismatch means `NOT_THE_SAME_FROZEN_SNAPSHOT`; it does not mean the object
#' is structurally or scientifically invalid. Use ordinary `validate_*()` for
#' structural/scientific validation.
#'
#' @param x A supported La Terra object.
#' @param fingerprint A prior `lt_audit_fingerprint`.
#' @return `lt_same_frozen_snapshot()` returns one logical value;
#'   `validate_lt_frozen_snapshot()` returns `x` invisibly or errors.
#' @export
lt_same_frozen_snapshot <- function(x, fingerprint) {
  .lt_validate_audit_fingerprint(fingerprint)
  current <- lt_audit_fingerprint(x)
  identical(current$object_class, fingerprint$object_class) &&
    identical(current$snapshot_sha256, fingerprint$snapshot_sha256)
}

#' @rdname lt_same_frozen_snapshot
#' @export
validate_lt_frozen_snapshot <- function(x, fingerprint) {
  if (!lt_same_frozen_snapshot(x, fingerprint)) {
    .lt_abort(paste(
      "NOT_THE_SAME_FROZEN_SNAPSHOT:",
      "the object may remain structurally valid; recertify only if audit identity is required."
    ))
  }
  invisible(x)
}

#' Check whether a rate still matches its scientific dependencies
#'
#' This operation-boundary check separates an editable input object from a
#' dependent result. Display-label and non-authoritative metadata edits do not
#' make a rate stale. Numeric values, scientific keys, coordinate state,
#' Layer-2 declarations, or baseline-fit identity changes do.
#'
#' @param rate An `lt_rate`.
#' @param matrix_object Candidate source `lt_matrix`.
#' @param baseline_fit Candidate `lt_c3_fit` or `lt_c2_fit`.
#' @param measurement_fit_eligibility Candidate Layer-2 object; defaults to the
#'   fit-domain declaration.
#' @return An `lt_dependency_status` with `CURRENT` or
#'   `STALE_RECOMPUTE_REQUIRED` state and machine-readable reasons.
#' @export
lt_rate_dependency_status <- function(
    rate, matrix_object, baseline_fit,
    measurement_fit_eligibility = baseline_fit$domain$eligibility) {
  validate_lt_rate(rate)
  validate_lt_matrix(matrix_object)
  validate_lt_measurement_eligibility(
    measurement_fit_eligibility, matrix_object
  )
  reasons <- character()
  source_identity <- .lt_hash(list(
    matrix_object$values, matrix_object$coord_state, matrix_object$value_reason,
    matrix_object$gene_ids, matrix_object$branch_ids
  ))
  if (!identical(source_identity,
                 rate$representation_provenance$source_matrix_identity_sha256)) {
    reasons <- c(reasons, "source_matrix_or_scientific_keys_changed")
  }
  if (!identical(
    .lt_hash(.lt_coordinate_provenance_semantic_payload(
      matrix_object$coordinate_provenance
    )), rate$output_hashes$coordinate_provenance_sha256
  )) reasons <- c(reasons, "coordinate_scientific_identity_changed")
  if (!identical(measurement_fit_eligibility$status_sha256,
                 rate$measurement_fit_eligibility$status_sha256) ||
      !identical(measurement_fit_eligibility$reason_sha256,
                 rate$measurement_fit_eligibility$reason_sha256) ||
      !identical(.lt_hash(measurement_fit_eligibility$provenance),
                 rate$measurement_fit_eligibility$provenance_sha256)) {
    reasons <- c(reasons, "measurement_layer2_changed")
  }
  if (inherits(baseline_fit, "lt_c3_fit")) {
    validate_lt_c3_fit(baseline_fit)
    mu_identity <- .lt_c3_mu_identity(baseline_fit)
    if (!identical(rate$baseline_estimand,
                   "C3_nonnegative_raw_scale_working_mean") ||
        !identical(rate$baseline_fit_id, baseline_fit$run_id) ||
        !identical(rate$baseline_mu_sha256,
                   mu_identity$certified_mu_sha256)) {
      reasons <- c(reasons, "baseline_fit_changed")
    }
  } else if (inherits(baseline_fit, "lt_c2_fit")) {
    validate_lt_c2_fit(baseline_fit)
    if (!identical(rate$baseline_estimand, "C2_positive_cell_log_OLS") ||
        !identical(rate$baseline_fit_id, baseline_fit$run_id) ||
        !identical(rate$baseline_mu_sha256,
                   baseline_fit$evaluation_mu_sha256)) {
      reasons <- c(reasons, "baseline_fit_changed")
    }
  } else {
    .lt_abort("`baseline_fit` must be lt_c3_fit or lt_c2_fit.")
  }
  reasons <- unique(reasons)
  structure(list(
    state = if (length(reasons)) "STALE_RECOMPUTE_REQUIRED" else "CURRENT",
    current = !length(reasons),
    reasons = reasons,
    action = if (length(reasons))
      "recompute dependent fit/rate or use a future explicit compatible rekey"
    else "safe_to_reuse",
    snapshot_identity_checked = FALSE
  ), class = "lt_dependency_status")
}

#' Require a current rate dependency
#'
#' @param ... Arguments passed to `lt_rate_dependency_status()`.
#' @return The dependency status invisibly, or an operation-boundary error.
#' @export
validate_lt_rate_dependency <- function(...) {
  status <- lt_rate_dependency_status(...)
  if (!isTRUE(status$current)) {
    .lt_abort(paste0(
      "STALE_RECOMPUTE_REQUIRED: ", paste(status$reasons, collapse = "; ")
    ))
  }
  invisible(status)
}

