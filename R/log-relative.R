# Boundary-aware log-relative representation and explicit historical migration.

.lt_ratio_state_levels <- function() {
  c("POSITIVE", "EXACT_ZERO_BOUNDARY", "UNAVAILABLE")
}

.lt_log_relative_reason_levels <- function() {
  c(
    "not_applicable",
    "baseline_zero",
    "baseline_unavailable",
    "numeric_indeterminate",
    "measurement_ineligible",
    "measurement_unresolved",
    "coordinate_unavailable",
    "observed_payload_unavailable",
    "cross_component_baseline_unavailable"
  )
}

.lt_log_relative_contract_version <- function() {
  "LaTerra_V1_log_relative_boundary_v1"
}

.lt_log_relative_canonicalization_version <- function() {
  "LaTerra_R_serialization_v3_raw_code_order_v1"
}

.lt_log_relative_identity_hash <- function(payload, canonicalization_version) {
  if (!identical(
    canonicalization_version,
    "LaTerra_R_serialization_v3_raw_code_order_v1"
  )) {
    .lt_abort_class(
      "Unknown log-relative canonicalization; identity computation fails closed.",
      "lt_error_log_relative_canonicalization_unknown"
    )
  }
  .lt_hash_serialization_v3(payload)
}

.lt_raw_matrix <- function(code, nrow, ncol, dimnames) {
  matrix(as.raw(code), nrow = nrow, ncol = ncol, dimnames = dimnames)
}

.lt_log_relative_source_identity <- function(gbi) {
  gbi$output_hashes$semantic_object_sha256 %||% .lt_hash(list(
    representation = gbi$representation,
    values = gbi$values,
    availability = gbi$representation_availability,
    gene_ids = gbi$ordered_gene_ledger,
    branch_ids = gbi$ordered_branch_ledger
  ))
}

.lt_log_relative_source_c3_identity <- function(gbi) {
  .lt_hash(list(
    baseline_fit_id = gbi$baseline_fit_id,
    baseline_mu_sha256 = gbi$baseline_mu_sha256,
    baseline_mu_all_admitted_sha256 =
      gbi$baseline_mu_all_admitted_sha256,
    baseline_estimand = gbi$baseline_estimand
  ))
}

.lt_log_relative_layer2_identity <- function(gbi) {
  .lt_hash(list(
    status_sha256 = gbi$measurement_fit_eligibility$status_sha256,
    reason_sha256 = gbi$measurement_fit_eligibility$reason_sha256,
    provenance_sha256 =
      gbi$measurement_fit_eligibility$provenance_sha256
  ))
}

.lt_log_relative_reason_labels_from_gbi <- function(gbi) {
  matrix(
    gbi$representation_reason_levels[gbi$representation_reason],
    nrow = nrow(gbi$values), ncol = ncol(gbi$values),
    dimnames = dimnames(gbi$values)
  )
}

.lt_gbi_zero_origin_from_embedded <- function(gbi) {
  if (is.null(gbi$zero_origin_state)) {
    .lt_abort_class(
      paste(
        "The source GBI does not carry certified numerator-zero origin.",
        "Use `lt_migrate_logGBI()` with the exact source matrix and C3 fit;",
        "post-division `GBI == 0` is not sufficient evidence."
      ),
      "lt_error_log_relative_zero_origin_required"
    )
  }
  gbi$zero_origin_state
}

.lt_certify_gbi_zero_origin <- function(gbi, source_matrix, c3_fit) {
  validate_lt_rate(gbi)
  if (!identical(gbi$representation, "GBI")) {
    .lt_abort("`source_gbi` must be a GBI lt_rate.")
  }
  mu_identity <- .lt_rate_identity_gate(
    source_matrix, c3_fit$domain$eligibility, c3_fit
  )
  if (!identical(gbi$ordered_gene_ledger, source_matrix$gene_ids) ||
      !identical(gbi$ordered_branch_ledger, source_matrix$branch_ids) ||
      !identical(gbi$baseline_fit_id, c3_fit$run_id) ||
      !identical(gbi$baseline_mu_sha256,
                 mu_identity$certified_mu_sha256) ||
      !identical(gbi$baseline_mu_all_admitted_sha256,
                 mu_identity$all_admitted_mu_sha256)) {
    .lt_abort_class(
      "Source GBI, source matrix, and C3 fit are not exact dependencies.",
      "lt_error_log_relative_source_mismatch"
    )
  }
  expected <- matrix(
    NA_real_, nrow(source_matrix$values), ncol(source_matrix$values),
    dimnames = dimnames(source_matrix$values)
  )
  edge <- c3_fit$edge_linear
  y <- c3_fit$observed_value
  mu <- c3_fit$mu
  available <- mu > 0
  expected[edge[available]] <- y[available] / mu[available]
  if (!identical(expected, gbi$values)) {
    .lt_abort_class(
      "Source GBI values are not the exact deterministic view of the supplied C3 authority.",
      "lt_error_log_relative_source_mismatch"
    )
  }
  origin <- .lt_raw_matrix(
    0L, nrow(expected), ncol(expected), dimnames(expected)
  )
  origin[edge[y > 0 & mu > 0]] <- as.raw(1L)
  origin[edge[y == 0 & mu > 0]] <- as.raw(2L)
  origin
}

.lt_log_relative_semantic_payload <- function(x) {
  list(
    representation_id = x$representation_id,
    positive_log_values = x$positive_log_values,
    ratio_state = x$ratio_state,
    unavailability_reason = x$unavailability_reason,
    ordered_gene_ledger = x$ordered_gene_ledger,
    ordered_branch_ledger = x$ordered_branch_ledger,
    domain_identity = x$domain_identity,
    source_GBI_identity = x$source_GBI_identity,
    source_baseline_identity = x$source_baseline_identity,
    source_coordinate_state_identity = x$source_coordinate_state_identity,
    measurement_eligibility_identity = x$measurement_eligibility_identity,
    ratio_state_vocabulary_sha256 = x$ratio_state_vocabulary_sha256,
    unavailability_reason_vocabulary_sha256 =
      x$unavailability_reason_vocabulary_sha256,
    representation_contract_version = x$representation_contract_version,
    canonicalization_version = x$canonicalization_version
  )
}

.lt_log_relative_audit_payload <- function(x) {
  list(
    semantic_identity = x$semantic_identity,
    source_GBI_identity = x$source_GBI_identity,
    source_GBI_values_sha256 = x$source_GBI_values_sha256,
    source_baseline_identity = x$source_baseline_identity,
    source_C3_fit_identity = x$source_C3_fit_identity,
    ordered_axes_identity = x$ordered_axes_identity,
    source_coordinate_state_identity = x$source_coordinate_state_identity,
    measurement_eligibility_identity = x$measurement_eligibility_identity,
    ratio_state_vocabulary_sha256 = x$ratio_state_vocabulary_sha256,
    unavailability_reason_vocabulary_sha256 =
      x$unavailability_reason_vocabulary_sha256,
    representation_contract_version = x$representation_contract_version,
    canonicalization_version = x$canonicalization_version,
    representation_provenance = x$representation_provenance,
    migration_provenance = x$migration_provenance
  )
}

.lt_new_log_relative <- function(gbi, zero_origin_state,
                                 migration_provenance = NULL,
                                 metadata = list()) {
  validate_lt_rate(gbi)
  if (!identical(gbi$representation, "GBI")) {
    .lt_abort("`gbi` must be a production GBI lt_rate.")
  }
  .lt_assert_named_list(metadata, "metadata")
  dn <- dimnames(gbi$values)
  if (!is.raw(zero_origin_state) || !is.matrix(zero_origin_state) ||
      !identical(dimnames(zero_origin_state), dn) ||
      any(as.integer(zero_origin_state) > 2L)) {
    .lt_abort("Certified zero-origin state is not a valid raw matrix on the GBI axes.")
  }
  if (any(gbi$values < 0, na.rm = TRUE) ||
      any(!is.finite(gbi$values[!is.na(gbi$values)]))) {
    .lt_abort_class(
      "GBI contains a negative or nonfinite available value.",
      "lt_error_log_relative_numeric_invalid"
    )
  }

  nr <- nrow(gbi$values); nc <- ncol(gbi$values)
  state <- .lt_raw_matrix(3L, nr, nc, dn)
  reason_labels <- .lt_log_relative_reason_labels_from_gbi(gbi)
  current_reason <- matrix(
    "numeric_indeterminate", nr, nc, dimnames = dn
  )
  reason_map <- reason_labels %in% .lt_log_relative_reason_levels()
  current_reason[reason_map] <- reason_labels[reason_map]

  positive <- gbi$representation_availability & gbi$values > 0
  certified_zero <- gbi$representation_availability & gbi$values == 0 &
    as.integer(zero_origin_state) == 2L
  uncertain_zero <- gbi$representation_availability & gbi$values == 0 &
    as.integer(zero_origin_state) != 2L
  state[positive] <- as.raw(1L)
  state[certified_zero] <- as.raw(2L)
  current_reason[positive | certified_zero] <- "not_applicable"
  current_reason[uncertain_zero] <- "numeric_indeterminate"

  positive_log <- matrix(NA_real_, nr, nc, dimnames = dn)
  positive_log[positive] <- log(gbi$values[positive])
  if (any(!is.finite(positive_log[positive]))) {
    .lt_abort_class(
      "A positive GBI value did not yield a finite mathematical log coordinate.",
      "lt_error_log_relative_numeric_invalid"
    )
  }
  reason_levels <- .lt_log_relative_reason_levels()
  reason_code <- match(current_reason, reason_levels)
  if (anyNA(reason_code)) {
    .lt_abort("An unavailable cell could not be mapped to the frozen reason vocabulary.")
  }
  reason_code <- matrix(
    as.raw(reason_code), nr, nc, dimnames = dn
  )
  ratio_available <- as.integer(state) %in% c(1L, 2L)
  dim(ratio_available) <- c(nr, nc); dimnames(ratio_available) <- dn
  zero_mask <- as.integer(state) == 2L
  dim(zero_mask) <- c(nr, nc); dimnames(zero_mask) <- dn
  positive_mask <- as.integer(state) == 1L
  dim(positive_mask) <- c(nr, nc); dimnames(positive_mask) <- dn

  m_ref <- if (any(positive_mask)) min(gbi$values[positive_mask]) else NA_real_
  L_min_plus <- if (any(positive_mask)) min(positive_log[positive_mask]) else
    NA_real_
  sentinel <- if (any(positive_mask)) log(m_ref) - 10 else NA_real_
  ratio_levels <- .lt_ratio_state_levels()
  source_gbi_identity <- .lt_log_relative_source_identity(gbi)
  source_baseline_identity <- .lt_hash(list(
    estimand = gbi$baseline_estimand,
    fit_id = gbi$baseline_fit_id,
    mu = gbi$baseline_mu_sha256,
    mu_all = gbi$baseline_mu_all_admitted_sha256
  ))
  coordinate_identity <- .lt_hash(list(
    coord_state = gbi$coord_state,
    coordinate_provenance = gbi$coordinate_provenance
  ))
  layer2_identity <- .lt_log_relative_layer2_identity(gbi)
  axes_identity <- .lt_hash(list(
    gene_ids = gbi$ordered_gene_ledger,
    branch_ids = gbi$ordered_branch_ledger
  ))
  representation_provenance <- list(
    operator = "LaTerra_lt_log_relative_from_GBI",
    exact_zero_origin = "certified_observed_numerator_exact_zero_with_positive_baseline",
    post_division_zero_inference = FALSE,
    floating_underflow_policy = "UNAVAILABLE_numeric_indeterminate",
    positive_log_rule = "natural_log_of_positive_GBI_only",
    default_finite_sentinel = FALSE,
    default_joint_p = FALSE,
    original_source_m_ref = m_ref,
    L_min_plus = L_min_plus,
    legacy_global_k10_sentinel = sentinel,
    positive_domain_identity = .lt_rate_keyed_mask_hash(
      positive_mask, gbi$ordered_gene_ledger, gbi$ordered_branch_ledger
    ),
    exact_zero_identity = .lt_rate_keyed_mask_hash(
      zero_mask, gbi$ordered_gene_ledger, gbi$ordered_branch_ledger
    ),
    legacy_mode_identity = .lt_hash(list(
      mode = "legacy_global_k10",
      contract = "LaTerra_V1_logGBI_zero_state_v1",
      m_ref = m_ref,
      L_min_plus = L_min_plus,
      sentinel = sentinel
    )),
    C3_or_rate_regeneration = FALSE
  )
  fields <- list(
    representation_id = "log_relative_boundary",
    positive_log_values = positive_log,
    ratio_state = state,
    ratio_state_levels = ratio_levels,
    unavailability_reason = reason_code,
    unavailability_reason_levels = reason_levels,
    ordered_gene_ledger = gbi$ordered_gene_ledger,
    ordered_branch_ledger = gbi$ordered_branch_ledger,
    gene_display_labels = gbi$gene_display_labels %||%
      gbi$ordered_gene_ledger,
    branch_display_labels = gbi$branch_display_labels %||%
      gbi$ordered_branch_ledger,
    domain_identity = .lt_rate_keyed_mask_hash(
      ratio_available, gbi$ordered_gene_ledger, gbi$ordered_branch_ledger
    ),
    source_GBI_identity = source_gbi_identity,
    source_GBI_values_sha256 = .lt_hash(gbi$values),
    source_baseline_identity = source_baseline_identity,
    source_C3_fit_identity = .lt_log_relative_source_c3_identity(gbi),
    source_coordinate_state_identity = coordinate_identity,
    measurement_fit_eligibility = gbi$measurement_fit_eligibility,
    measurement_eligibility_identity = layer2_identity,
    ordered_axes_identity = axes_identity,
    ratio_state_vocabulary_sha256 = .lt_hash(ratio_levels),
    unavailability_reason_vocabulary_sha256 = .lt_hash(reason_levels),
    representation_contract_version = .lt_log_relative_contract_version(),
    canonicalization_version = .lt_log_relative_canonicalization_version(),
    representation_provenance = representation_provenance,
    migration_provenance = migration_provenance,
    semantic_identity = NULL,
    audit_identity = NULL,
    metadata = metadata
  )
  out <- structure(fields, class = "lt_log_relative")
  out$semantic_identity <- .lt_log_relative_identity_hash(
    .lt_log_relative_semantic_payload(out), out$canonicalization_version
  )
  out$audit_identity <- .lt_log_relative_identity_hash(
    .lt_log_relative_audit_payload(out), out$canonicalization_version
  )
  validate_lt_log_relative(out)
  out
}

#' Construct a boundary-aware log-relative object
#'
#' `lt_log_relative()` represents positive log coordinates and true observed
#' exact-zero boundaries without assigning a finite sentinel. Exact-zero
#' classification requires certified numerator-zero provenance carried by the
#' source GBI.
#'
#' “Exact zero” means a numeric exact zero preserved in the admitted observed
#' payload under the frozen input/eligibility contract. It does not establish
#' an ontological biological rate of exactly zero, exclude a tiny positive
#' quantity before admitted serialization, or prove optimizer-estimated strict
#' zero. Binary analysis uses boundary-aware two-component estimands; V1 does
#' not implement a full generative hurdle model.
#'
#' @param gbi A production `lt_rate` with representation `GBI` and certified
#'   zero-origin state.
#' @param metadata Named, non-authoritative annotations.
#' @return An `lt_log_relative` object.
#' @export
lt_log_relative <- function(gbi, metadata = list()) {
  .lt_new_log_relative(
    gbi, .lt_gbi_zero_origin_from_embedded(gbi), metadata = metadata
  )
}

#' Validate a boundary-aware log-relative object
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_log_relative <- function(x) {
  if (!inherits(x, "lt_log_relative")) {
    .lt_abort("`x` must inherit from `lt_log_relative`.")
  }
  if (!identical(x$representation_id, "log_relative_boundary") ||
      !identical(x$ratio_state_levels, .lt_ratio_state_levels()) ||
      !identical(x$unavailability_reason_levels,
                 .lt_log_relative_reason_levels()) ||
      !identical(x$representation_contract_version,
                 .lt_log_relative_contract_version()) ||
      !identical(x$canonicalization_version,
                 .lt_log_relative_canonicalization_version())) {
    .lt_abort("Log-relative representation/vocabulary contract is unknown.")
  }
  .lt_assert_unique_ids(x$ordered_gene_ledger, "ordered_gene_ledger")
  .lt_assert_unique_ids(x$ordered_branch_ledger, "ordered_branch_ledger")
  dn <- list(x$ordered_gene_ledger, x$ordered_branch_ledger)
  matrices <- c("positive_log_values", "ratio_state",
                "unavailability_reason")
  if (any(!vapply(x[matrices], is.matrix, logical(1L))) ||
      any(!vapply(x[matrices], function(z) identical(dimnames(z), dn),
                  logical(1L))) ||
      !is.numeric(x$positive_log_values) || !is.raw(x$ratio_state) ||
      !is.raw(x$unavailability_reason)) {
    .lt_abort("Primary log-relative payload must use aligned numeric/raw matrices.")
  }
  state <- as.integer(x$ratio_state)
  reason <- as.integer(x$unavailability_reason)
  if (any(!state %in% 1:3) ||
      any(!reason %in% seq_along(x$unavailability_reason_levels))) {
    .lt_abort("Log-relative compact code is outside its frozen vocabulary.")
  }
  positive <- state == 1L
  boundary <- state == 2L
  unavailable <- state == 3L
  if (any(!is.finite(x$positive_log_values[positive])) ||
      any(!is.na(x$positive_log_values[!positive]))) {
    .lt_abort("Finite positive_log_values are permitted only on POSITIVE cells.")
  }
  reason_label <- x$unavailability_reason_levels[reason]
  if (any(reason_label[positive | boundary] != "not_applicable") ||
      any(reason_label[unavailable] == "not_applicable")) {
    .lt_abort("Mathematical ratio state and unavailability reason disagree.")
  }
  if (!identical(x$ratio_state_vocabulary_sha256,
                 .lt_hash(x$ratio_state_levels)) ||
      !identical(x$unavailability_reason_vocabulary_sha256,
                 .lt_hash(x$unavailability_reason_levels))) {
    .lt_abort("Log-relative vocabulary identity is stale.")
  }
  expected_axes <- .lt_hash(list(
    gene_ids = x$ordered_gene_ledger,
    branch_ids = x$ordered_branch_ledger
  ))
  ratio_available <- state %in% c(1L, 2L)
  dim(ratio_available) <- dim(x$ratio_state)
  dimnames(ratio_available) <- dn
  expected_domain <- .lt_rate_keyed_mask_hash(
    ratio_available, x$ordered_gene_ledger, x$ordered_branch_ledger
  )
  if (!identical(x$ordered_axes_identity, expected_axes) ||
      !identical(x$domain_identity, expected_domain)) {
    .lt_abort("Log-relative ordered-axis/domain identity is stale.")
  }
  .lt_assert_named_list(x$representation_provenance,
                        "representation_provenance")
  provenance <- x$representation_provenance
  required_legacy <- c(
    "original_source_m_ref", "L_min_plus",
    "legacy_global_k10_sentinel", "positive_domain_identity",
    "exact_zero_identity", "legacy_mode_identity"
  )
  if (length(setdiff(required_legacy, names(provenance)))) {
    .lt_abort("Log-relative legacy replay provenance is incomplete.")
  }
  positive_mask <- matrix(
    positive, nrow(x$ratio_state), ncol(x$ratio_state), dimnames = dn
  )
  zero_mask <- matrix(
    boundary, nrow(x$ratio_state), ncol(x$ratio_state), dimnames = dn
  )
  if (!identical(
    provenance$positive_domain_identity,
    .lt_rate_keyed_mask_hash(
      positive_mask, x$ordered_gene_ledger, x$ordered_branch_ledger
    )
  ) || !identical(
    provenance$exact_zero_identity,
    .lt_rate_keyed_mask_hash(
      zero_mask, x$ordered_gene_ledger, x$ordered_branch_ledger
    )
  )) {
    .lt_abort("Log-relative positive/zero provenance identities are stale.")
  }
  if (any(positive)) {
    if (!is.numeric(provenance$original_source_m_ref) ||
        length(provenance$original_source_m_ref) != 1L ||
        !is.finite(provenance$original_source_m_ref) ||
        provenance$original_source_m_ref <= 0 ||
        !identical(provenance$L_min_plus,
                   log(provenance$original_source_m_ref)) ||
        !identical(provenance$L_min_plus,
                   min(x$positive_log_values[positive])) ||
        !identical(provenance$legacy_global_k10_sentinel,
                   provenance$L_min_plus - 10)) {
      .lt_abort("Log-relative source m_ref/L_min_plus/sentinel provenance is inconsistent.")
    }
  } else if (!all(is.na(c(
    provenance$original_source_m_ref, provenance$L_min_plus,
    provenance$legacy_global_k10_sentinel
  )))) {
    .lt_abort("A log-relative object without positive support cannot define legacy m_ref.")
  }
  expected_legacy_mode_identity <- .lt_hash(list(
    mode = "legacy_global_k10",
    contract = "LaTerra_V1_logGBI_zero_state_v1",
    m_ref = provenance$original_source_m_ref,
    L_min_plus = provenance$L_min_plus,
    sentinel = provenance$legacy_global_k10_sentinel
  ))
  if (!identical(provenance$legacy_mode_identity,
                 expected_legacy_mode_identity)) {
    .lt_abort("Log-relative legacy mode/version identity is stale.")
  }
  eligibility <- x$measurement_fit_eligibility
  if (!is.list(eligibility) ||
      !identical(eligibility$status_sha256, .lt_hash(eligibility$status)) ||
      !identical(eligibility$reason_sha256, .lt_hash(eligibility$reason)) ||
      !identical(eligibility$provenance_sha256,
                 .lt_hash(eligibility$provenance)) ||
      !identical(dimnames(eligibility$status), dn) ||
      !identical(dimnames(eligibility$reason), dn) ||
      !identical(x$measurement_eligibility_identity, .lt_hash(list(
        status_sha256 = eligibility$status_sha256,
        reason_sha256 = eligibility$reason_sha256,
        provenance_sha256 = eligibility$provenance_sha256
      )))) {
    .lt_abort("Embedded Layer-2 dependency is malformed or stale.")
  }
  hashes <- c(
    "domain_identity", "source_GBI_identity", "source_GBI_values_sha256",
    "source_baseline_identity", "source_C3_fit_identity",
    "source_coordinate_state_identity", "measurement_eligibility_identity",
    "ordered_axes_identity", "ratio_state_vocabulary_sha256",
    "unavailability_reason_vocabulary_sha256", "semantic_identity",
    "audit_identity"
  )
  invisible(lapply(hashes, function(name) .lt_assert_sha256(x[[name]], name)))
  if (!is.null(x$migration_provenance)) {
    .lt_assert_named_list(x$migration_provenance, "migration_provenance")
  }
  .lt_assert_named_list(x$metadata, "metadata")
  if (!identical(x$semantic_identity,
                 .lt_log_relative_identity_hash(
                   .lt_log_relative_semantic_payload(x),
                   x$canonicalization_version
                 )) ||
      !identical(x$audit_identity,
                 .lt_log_relative_identity_hash(
                   .lt_log_relative_audit_payload(x),
                   x$canonicalization_version
                 ))) {
    .lt_abort("Log-relative semantic or audit identity is stale.")
  }
  invisible(x)
}

#' Materialize log-relative state and component accessors
#' @param x An `lt_log_relative` object.
#' @return An aligned matrix.
#' @export
ratio_state <- function(x) {
  validate_lt_log_relative(x)
  matrix(
    x$ratio_state_levels[as.integer(x$ratio_state)],
    nrow = nrow(x$ratio_state), ncol = ncol(x$ratio_state),
    dimnames = dimnames(x$ratio_state)
  )
}

#' @rdname ratio_state
#' @export
unavailability_reason <- function(x) {
  validate_lt_log_relative(x)
  matrix(
    x$unavailability_reason_levels[as.integer(x$unavailability_reason)],
    nrow = nrow(x$unavailability_reason),
    ncol = ncol(x$unavailability_reason),
    dimnames = dimnames(x$unavailability_reason)
  )
}

#' @rdname ratio_state
#' @export
exact_zero <- function(x) {
  validate_lt_log_relative(x)
  state <- as.integer(x$ratio_state)
  out <- rep(NA, length(state))
  out[state == 1L] <- FALSE
  out[state == 2L] <- TRUE
  matrix(out, nrow = nrow(x$ratio_state), ncol = ncol(x$ratio_state),
         dimnames = dimnames(x$ratio_state))
}

#' @rdname ratio_state
#' @export
ratio_availability <- function(x) {
  validate_lt_log_relative(x)
  out <- as.integer(x$ratio_state) %in% c(1L, 2L)
  matrix(out, nrow = nrow(x$ratio_state), ncol = ncol(x$ratio_state),
         dimnames = dimnames(x$ratio_state))
}

#' @rdname ratio_state
#' @export
positive_log_availability <- function(x) {
  validate_lt_log_relative(x)
  out <- as.integer(x$ratio_state) == 1L
  matrix(out, nrow = nrow(x$ratio_state), ncol = ncol(x$ratio_state),
         dimnames = dimnames(x$ratio_state))
}

#' @rdname ratio_state
#' @export
positive_log_values <- function(x) {
  validate_lt_log_relative(x)
  x$positive_log_values
}

#' Check freshness against a source GBI dependency
#' @param x An `lt_log_relative` object.
#' @param source_gbi Candidate source GBI.
#' @param error Whether to throw a classed stale-dependency error.
#' @return A one-row dependency status data frame.
#' @export
lt_log_relative_dependency_status <- function(x, source_gbi, error = FALSE) {
  validate_lt_log_relative(x)
  validate_lt_rate(source_gbi)
  current <- identical(source_gbi$representation, "GBI") &&
    identical(x$source_GBI_identity,
              .lt_log_relative_source_identity(source_gbi)) &&
    identical(x$source_GBI_values_sha256, .lt_hash(source_gbi$values)) &&
    identical(x$source_baseline_identity, .lt_hash(list(
      estimand = source_gbi$baseline_estimand,
      fit_id = source_gbi$baseline_fit_id,
      mu = source_gbi$baseline_mu_sha256,
      mu_all = source_gbi$baseline_mu_all_admitted_sha256
    ))) &&
    identical(x$source_C3_fit_identity,
              .lt_log_relative_source_c3_identity(source_gbi)) &&
    identical(x$measurement_eligibility_identity,
              .lt_log_relative_layer2_identity(source_gbi))
  status <- if (current) "CURRENT" else "STALE_DEPENDENCY"
  if (!current && isTRUE(error)) {
    .lt_abort_class(
      "The lt_log_relative object is stale relative to the supplied source GBI.",
      "lt_error_stale_log_relative_dependency"
    )
  }
  data.frame(
    status = status,
    expected_source_GBI_identity = x$source_GBI_identity,
    observed_source_GBI_identity = .lt_log_relative_source_identity(source_gbi),
    stringsAsFactors = FALSE
  )
}

#' Explicitly migrate historical logGBI objects
#'
#' Migration never infers a boundary from sentinel equality. Historical GBI
#' objects lacking embedded zero-origin provenance require the exact source
#' matrix and C3 fit.
#'
#' @param old Historical `lt_rate` with representation `logGBI` or
#'   `logGBI_strict`.
#' @param source_gbi Exact compatible GBI authority.
#' @param source_matrix Exact source `lt_matrix`, required when `source_gbi`
#'   lacks embedded zero-origin provenance.
#' @param c3_fit Exact source `lt_c3_fit`, required with `source_matrix`.
#' @param source_association_domain Optional exact historical keyed domain used
#'   only to retain replay identity; ordinary runtime migration does not require
#'   it.
#' @param source_reference Optional immutable source path/reference.
#' @param source_file_sha256 Optional source-file SHA-256.
#' @param metadata Named non-authoritative annotations.
#' @return A new `lt_log_relative`; `old` is never modified.
#' @export
lt_migrate_logGBI <- function(old, source_gbi, source_matrix = NULL,
                              c3_fit = NULL, source_association_domain = NULL,
                              source_reference = NA_character_,
                              source_file_sha256 = NA_character_,
                              metadata = list()) {
  if (missing(source_gbi) || is.null(source_gbi)) {
    .lt_abort_class(
      "Scientific logGBI migration requires an exact source GBI authority.",
      "lt_error_logGBI_migration_source_required"
    )
  }
  validate_lt_rate(old)
  validate_lt_rate(source_gbi)
  if (!old$representation %in% c("logGBI", "logGBI_strict") ||
      !identical(source_gbi$representation, "GBI")) {
    .lt_abort_class(
      "Historical object or source GBI contract is not migratable.",
      "lt_error_unclassified_historical_logGBI"
    )
  }
  if (!identical(old$ordered_gene_ledger, source_gbi$ordered_gene_ledger) ||
      !identical(old$ordered_branch_ledger, source_gbi$ordered_branch_ledger) ||
      !identical(old$baseline_mu_sha256, source_gbi$baseline_mu_sha256)) {
    .lt_abort_class(
      "Historical logGBI and source GBI are not exact compatible dependencies.",
      "lt_error_log_relative_source_mismatch"
    )
  }
  origin <- if (!is.null(source_gbi$zero_origin_state)) {
    source_gbi$zero_origin_state
  } else {
    if (is.null(source_matrix) || is.null(c3_fit)) {
      .lt_abort_class(
        paste(
          "This historical source GBI lacks certified zero-origin state;",
          "supply the exact source matrix and C3 fit."
        ),
        "lt_error_log_relative_zero_origin_required"
      )
    }
    .lt_certify_gbi_zero_origin(source_gbi, source_matrix, c3_fit)
  }
  positive <- source_gbi$representation_availability & source_gbi$values > 0
  if (any(positive)) {
    expected <- log(source_gbi$values[positive])
    if (!identical(old$values[positive], expected)) {
      .lt_abort_class(
        "Historical positive coordinates differ from log(source GBI).",
        "lt_error_log_relative_source_mismatch"
      )
    }
  }
  source_reason <- .lt_log_relative_reason_labels_from_gbi(old)
  replay_domain <- list(
    association_domain_identity = NA_character_,
    cell_domain_identity = NA_character_
  )
  if (!is.null(source_association_domain)) {
    validate_lt_association_domain(source_association_domain)
    if (!identical(source_association_domain$gene_ids,
                   source_gbi$ordered_gene_ledger) ||
        !identical(source_association_domain$branch_ids,
                   source_gbi$ordered_branch_ledger) ||
        !identical(source_association_domain$cell_domain_identity,
                   .lt_rate_keyed_mask_hash(
                     source_gbi$representation_availability,
                     source_gbi$ordered_gene_ledger,
                     source_gbi$ordered_branch_ledger
                   ))) {
      .lt_abort_class(
        "Historical association domain is not the exact compatible GBI cell domain.",
        "lt_error_log_relative_source_mismatch"
      )
    }
    replay_domain <- list(
      association_domain_identity =
        source_association_domain$association_domain_identity,
      cell_domain_identity = source_association_domain$cell_domain_identity
    )
  }
  detected <- .lt_rate_log_contract(old)
  migration <- list(
    operator = "lt_migrate_logGBI_explicit_v1",
    automatic = FALSE,
    original_overwritten = FALSE,
    detected_historical_contract = detected,
    source_reference = source_reference,
    source_file_sha256 = source_file_sha256,
    source_semantic_identity = old$output_hashes$semantic_object_sha256,
    source_reason_vocabulary = old$representation_reason_levels,
    source_reason_vocabulary_sha256 = .lt_hash(old$representation_reason_levels),
    source_representation_reason_sha256 = .lt_hash(source_reason),
    state_inferred_from_sentinel_equality = FALSE,
    source_GBI_required = TRUE,
    source_authority_certified = TRUE,
    source_association_domain_identity =
      replay_domain$association_domain_identity,
    source_cell_domain_identity = replay_domain$cell_domain_identity,
    warnings = character()
  )
  .lt_new_log_relative(
    source_gbi, origin, migration_provenance = migration, metadata = metadata
  )
}

.lt_legacy_encoded_payload <- function(x) {
  x[c(
    "view_name", "mode", "values", "ratio_state", "ratio_state_levels",
    "ordered_gene_ledger", "ordered_branch_ledger",
    "source_semantic_identity", "source_audit_identity",
    "source_representation_id", "source_representation_contract_version",
    "source_canonicalization_version", "source_domain_identity",
    "source_ordered_axes_identity", "source_positive_log_values_sha256",
    "source_ratio_state_sha256", "source_positive_domain_identity",
    "source_exact_zero_identity", "source_object_identity",
    "encoding_contract_version", "encoding_rule", "legacy_mode_identity",
    "m_ref", "L_min_plus", "sentinel", "m_ref_provenance",
    "m_ref_provenance_identity", "positive_log_values_identity",
    "positive_domain_identity", "exact_zero_identity", "view_domain_identity",
    "view_ordered_axes_identity",
    "source_association_domain_identity", "source_cell_domain_identity",
    "mathematical_log_at_encoded_zero",
    "default_continuous_inference_authorized", "provenance"
  )]
}

.lt_legacy_encoded_positive_identity <- function(
    values, ratio_state, gene_ids, branch_ids) {
  positive_values <- matrix(
    NA_real_, nrow(values), ncol(values),
    dimnames = list(gene_ids, branch_ids)
  )
  positive <- as.integer(ratio_state) == 1L
  positive_values[positive] <- values[positive]
  .lt_hash(list(
    positive_log_values = positive_values,
    gene_ids = gene_ids,
    branch_ids = branch_ids
  ))
}

.lt_legacy_source_object_identity <- function(x) {
  .lt_hash(list(
    representation_id = x$source_representation_id,
    representation_contract_version =
      x$source_representation_contract_version,
    canonicalization_version = x$source_canonicalization_version,
    semantic_identity = x$source_semantic_identity,
    audit_identity = x$source_audit_identity,
    domain_identity = x$source_domain_identity,
    ordered_axes_identity = x$source_ordered_axes_identity,
    positive_log_values_sha256 = x$source_positive_log_values_sha256,
    ratio_state_sha256 = x$source_ratio_state_sha256,
    positive_domain_identity = x$source_positive_domain_identity,
    exact_zero_identity = x$source_exact_zero_identity
  ))
}

.lt_legacy_m_ref_provenance <- function(x) {
  list(
    source_object_identity = x$source_object_identity,
    source_semantic_identity = x$source_semantic_identity,
    source_audit_identity = x$source_audit_identity,
    source_domain_identity = x$source_domain_identity,
    source_positive_domain_identity = x$source_positive_domain_identity,
    original_source_m_ref = x$m_ref,
    L_min_plus = x$L_min_plus,
    sentinel = x$sentinel,
    extraction_rule = "minimum_positive_GBI_from_full_source_authority",
    subset_recomputation = FALSE
  )
}

.lt_match_legacy_view_axis <- function(requested, source, name) {
  if (is.null(requested)) return(seq_along(source))
  .lt_assert_unique_ids(requested, name)
  if (!length(requested)) {
    .lt_abort(sprintf("`%s` must retain at least one scientific key.", name))
  }
  index <- match(requested, source)
  if (anyNA(index)) {
    .lt_abort(sprintf("`%s` contains unknown scientific keys.", name))
  }
  index
}

#' Construct the explicit historical finite-encoded logGBI view
#' @param object An `lt_log_relative`.
#' @param mode The sole supported historical mode, `legacy_global_k10`.
#' @param gene_ids Optional exact ordered scientific gene-key subset.
#' @param branch_ids Optional exact ordered scientific branch-key subset.
#' @return An `lt_logGBI_encoded` derived view. It is not a primary authority.
#' @export
logGBI_encoded <- function(object, mode = "legacy_global_k10",
                           gene_ids = NULL, branch_ids = NULL) {
  validate_lt_log_relative(object)
  if (!identical(mode, "legacy_global_k10")) {
    .lt_abort_class(
      "Only the frozen `legacy_global_k10` encoded mode is supported.",
      "lt_error_unsupported_logGBI_encoded_mode"
    )
  }
  m_ref <- object$representation_provenance$original_source_m_ref
  L_min_plus <- object$representation_provenance$L_min_plus
  sentinel <- object$representation_provenance$legacy_global_k10_sentinel
  if (!is.numeric(m_ref) || length(m_ref) != 1L || !is.finite(m_ref) ||
      m_ref <= 0 || !is.finite(L_min_plus) || !is.finite(sentinel)) {
    .lt_abort_class(
      "legacy_global_k10 is NONESTIMABLE because the object has no positive reference cell.",
      "lt_error_logGBI_encoded_nonestimable"
    )
  }
  gene_index <- .lt_match_legacy_view_axis(
    gene_ids, object$ordered_gene_ledger, "gene_ids"
  )
  branch_index <- .lt_match_legacy_view_axis(
    branch_ids, object$ordered_branch_ledger, "branch_ids"
  )
  view_gene_ids <- object$ordered_gene_ledger[gene_index]
  view_branch_ids <- object$ordered_branch_ledger[branch_index]
  ratio_state <- object$ratio_state[
    gene_index, branch_index, drop = FALSE
  ]
  values <- object$positive_log_values[
    gene_index, branch_index, drop = FALSE
  ]
  state <- as.integer(ratio_state)
  values[state == 2L] <- sentinel
  positive_mask <- state == 1L
  zero_mask <- state == 2L
  available_mask <- state %in% c(1L, 2L)
  dim(positive_mask) <- dim(values); dimnames(positive_mask) <- dimnames(values)
  dim(zero_mask) <- dim(values); dimnames(zero_mask) <- dimnames(values)
  dim(available_mask) <- dim(values); dimnames(available_mask) <- dimnames(values)
  source_fields <- list(
    source_semantic_identity = object$semantic_identity,
    source_audit_identity = object$audit_identity,
    source_representation_id = object$representation_id,
    source_representation_contract_version =
      object$representation_contract_version,
    source_canonicalization_version = object$canonicalization_version,
    source_domain_identity = object$domain_identity,
    source_ordered_axes_identity = object$ordered_axes_identity,
    source_positive_log_values_sha256 = .lt_hash(object$positive_log_values),
    source_ratio_state_sha256 = .lt_hash(object$ratio_state),
    source_positive_domain_identity =
      object$representation_provenance$positive_domain_identity,
    source_exact_zero_identity =
      object$representation_provenance$exact_zero_identity
  )
  source_fields$source_object_identity <-
    .lt_legacy_source_object_identity(source_fields)
  core <- list(
    view_name = "logGBI_encoded",
    mode = mode,
    values = values,
    ratio_state = ratio_state,
    ratio_state_levels = object$ratio_state_levels,
    ordered_gene_ledger = view_gene_ids,
    ordered_branch_ledger = view_branch_ids,
    source_semantic_identity = source_fields$source_semantic_identity,
    source_audit_identity = source_fields$source_audit_identity,
    source_representation_id = source_fields$source_representation_id,
    source_representation_contract_version =
      source_fields$source_representation_contract_version,
    source_canonicalization_version =
      source_fields$source_canonicalization_version,
    source_domain_identity = source_fields$source_domain_identity,
    source_ordered_axes_identity = source_fields$source_ordered_axes_identity,
    source_positive_log_values_sha256 =
      source_fields$source_positive_log_values_sha256,
    source_ratio_state_sha256 = source_fields$source_ratio_state_sha256,
    source_positive_domain_identity =
      source_fields$source_positive_domain_identity,
    source_exact_zero_identity = source_fields$source_exact_zero_identity,
    source_object_identity = source_fields$source_object_identity,
    encoding_contract_version = "LaTerra_V1_logGBI_zero_state_v1",
    encoding_rule = "global_positive_log_min_minus_10",
    legacy_mode_identity =
      object$representation_provenance$legacy_mode_identity,
    m_ref = m_ref,
    L_min_plus = L_min_plus,
    sentinel = sentinel,
    m_ref_provenance = NULL,
    m_ref_provenance_identity = NULL,
    positive_log_values_identity = .lt_legacy_encoded_positive_identity(
      values, ratio_state, view_gene_ids, view_branch_ids
    ),
    positive_domain_identity = .lt_rate_keyed_mask_hash(
      positive_mask, view_gene_ids, view_branch_ids
    ),
    exact_zero_identity = .lt_rate_keyed_mask_hash(
      zero_mask, view_gene_ids, view_branch_ids
    ),
    view_domain_identity = .lt_rate_keyed_mask_hash(
      available_mask, view_gene_ids, view_branch_ids
    ),
    view_ordered_axes_identity = .lt_hash(list(
      gene_ids = view_gene_ids, branch_ids = view_branch_ids
    )),
    source_association_domain_identity =
      object$migration_provenance$source_association_domain_identity %||%
      NA_character_,
    source_cell_domain_identity =
      object$migration_provenance$source_cell_domain_identity %||%
      object$domain_identity,
    mathematical_log_at_encoded_zero = FALSE,
    default_continuous_inference_authorized = FALSE,
    provenance = list(
      derived_view = TRUE,
      primary_authority = FALSE,
      explicit_request_required = TRUE,
      positive_values_reused_without_log_exp_roundtrip = TRUE,
      source_gene_count = length(object$ordered_gene_ledger),
      source_branch_count = length(object$ordered_branch_ledger),
      view_gene_count = length(view_gene_ids),
      view_branch_count = length(view_branch_ids),
      subset_applied = length(gene_index) !=
        length(object$ordered_gene_ledger) || length(branch_index) !=
        length(object$ordered_branch_ledger),
      m_ref_recomputed_after_subset = FALSE
    )
  )
  core$m_ref_provenance <- .lt_legacy_m_ref_provenance(core)
  core$m_ref_provenance_identity <- .lt_hash(core$m_ref_provenance)
  out <- structure(c(core, list(
    encoded_identity = .lt_hash(core), metadata = list()
  )), class = "lt_logGBI_encoded")
  validate_lt_logGBI_encoded(out)
  out
}

#' Structurally validate a historical finite-encoded view
#'
#' This checks local shape, state, numeric, and self-hash consistency only. It
#' does not authenticate source provenance, certify frozen replay identity, or
#' authorize inference. Use [validate_lt_logGBI_encoded_source()] with a
#' separately supplied boundary source for source-bound validation.
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_logGBI_encoded <- function(x) {
  if (!inherits(x, "lt_logGBI_encoded") ||
      !identical(x$view_name, "logGBI_encoded") ||
      !identical(x$mode, "legacy_global_k10") ||
      !identical(x$encoding_contract_version,
                 "LaTerra_V1_logGBI_zero_state_v1") ||
      !identical(x$encoding_rule, "global_positive_log_min_minus_10") ||
      !identical(x$ratio_state_levels, .lt_ratio_state_levels()) ||
      !identical(x$source_representation_id, "log_relative_boundary") ||
      !identical(x$source_representation_contract_version,
                 .lt_log_relative_contract_version()) ||
      !identical(x$source_canonicalization_version,
                 .lt_log_relative_canonicalization_version()) ||
      !identical(x$mathematical_log_at_encoded_zero, FALSE) ||
      !identical(x$default_continuous_inference_authorized, FALSE)) {
    .lt_abort("Legacy finite-encoded view contract is malformed.")
  }
  .lt_assert_unique_ids(x$ordered_gene_ledger, "ordered_gene_ledger")
  .lt_assert_unique_ids(x$ordered_branch_ledger, "ordered_branch_ledger")
  if (!is.matrix(x$values) || !is.numeric(x$values) ||
      !is.matrix(x$ratio_state) || !is.raw(x$ratio_state) ||
      !identical(dimnames(x$values),
                 list(x$ordered_gene_ledger, x$ordered_branch_ledger)) ||
      !identical(dimnames(x$ratio_state), dimnames(x$values)) ||
      !identical(dim(x$ratio_state), dim(x$values)) ||
      any(is.infinite(x$values), na.rm = TRUE) || any(is.nan(x$values))) {
    .lt_abort("Legacy finite-encoded values are malformed.")
  }
  state <- as.integer(x$ratio_state)
  if (any(!state %in% 1:3)) {
    .lt_abort("Legacy finite-encoded ratio state is outside its vocabulary.")
  }
  positive <- state == 1L
  boundary <- state == 2L
  unavailable <- state == 3L
  if (!is.numeric(x$m_ref) || length(x$m_ref) != 1L ||
      !is.finite(x$m_ref) || x$m_ref <= 0 ||
      !is.numeric(x$L_min_plus) || length(x$L_min_plus) != 1L ||
      !is.finite(x$L_min_plus) ||
      !is.numeric(x$sentinel) || length(x$sentinel) != 1L ||
      !is.finite(x$sentinel) ||
      !identical(x$L_min_plus, log(x$m_ref)) ||
      !identical(x$sentinel, x$L_min_plus - 10)) {
    .lt_abort("Legacy global-k10 m_ref/L_min_plus/sentinel invariants fail.")
  }
  if (any(!is.finite(x$values[positive])) ||
      (any(boundary) &&
       (anyNA(x$values[boundary]) ||
        any(x$values[boundary] != x$sentinel))) ||
      any(!is.na(x$values[unavailable]))) {
    .lt_abort("Legacy encoded values disagree with positive, boundary, or unavailable states.")
  }
  hashes <- c("source_semantic_identity", "source_audit_identity",
              "source_domain_identity", "source_ordered_axes_identity",
              "source_positive_log_values_sha256", "source_ratio_state_sha256",
              "source_positive_domain_identity", "source_exact_zero_identity",
              "source_object_identity", "legacy_mode_identity",
              "m_ref_provenance_identity", "positive_log_values_identity",
              "positive_domain_identity", "exact_zero_identity",
              "view_domain_identity", "view_ordered_axes_identity",
              "encoded_identity")
  invisible(lapply(hashes, function(name) .lt_assert_sha256(x[[name]], name)))
  optional_hashes <- c("source_association_domain_identity",
                       "source_cell_domain_identity")
  invisible(lapply(optional_hashes, function(name) {
    value <- x[[name]]
    if (!(is.character(value) && length(value) == 1L && is.na(value))) {
      .lt_assert_sha256(value, name)
    }
  }))
  positive_mask <- matrix(
    positive, nrow(x$values), ncol(x$values), dimnames = dimnames(x$values)
  )
  zero_mask <- matrix(
    boundary, nrow(x$values), ncol(x$values), dimnames = dimnames(x$values)
  )
  available_mask <- matrix(
    positive | boundary, nrow(x$values), ncol(x$values),
    dimnames = dimnames(x$values)
  )
  if (!identical(x$positive_log_values_identity,
                 .lt_legacy_encoded_positive_identity(
                   x$values, x$ratio_state, x$ordered_gene_ledger,
                   x$ordered_branch_ledger
                 )) ||
      !identical(x$positive_domain_identity,
                 .lt_rate_keyed_mask_hash(
                   positive_mask, x$ordered_gene_ledger,
                   x$ordered_branch_ledger
                 )) ||
      !identical(x$exact_zero_identity,
                 .lt_rate_keyed_mask_hash(
                   zero_mask, x$ordered_gene_ledger,
                   x$ordered_branch_ledger
                 )) ||
      !identical(x$view_domain_identity,
                 .lt_rate_keyed_mask_hash(
                   available_mask, x$ordered_gene_ledger,
                   x$ordered_branch_ledger
                 )) ||
      !identical(x$view_ordered_axes_identity, .lt_hash(list(
        gene_ids = x$ordered_gene_ledger,
        branch_ids = x$ordered_branch_ledger
      )))) {
    .lt_abort("Legacy encoded keyed value/state identities are stale.")
  }
  if (!identical(x$source_object_identity,
                 .lt_legacy_source_object_identity(x))) {
    .lt_abort("Legacy encoded source-object identity is stale.")
  }
  expected_mode_identity <- .lt_hash(list(
    mode = "legacy_global_k10",
    contract = x$encoding_contract_version,
    m_ref = x$m_ref,
    L_min_plus = x$L_min_plus,
    sentinel = x$sentinel
  ))
  if (!identical(x$legacy_mode_identity, expected_mode_identity)) {
    .lt_abort("Legacy encoded mode/version identity is inconsistent.")
  }
  .lt_assert_named_list(x$m_ref_provenance, "m_ref_provenance")
  if (!identical(x$m_ref_provenance,
                 .lt_legacy_m_ref_provenance(x)) ||
      !identical(x$m_ref_provenance_identity,
                 .lt_hash(x$m_ref_provenance))) {
    .lt_abort("Legacy encoded source-domain m_ref provenance is inconsistent.")
  }
  .lt_assert_named_list(x$provenance, "provenance")
  if (!identical(x$provenance$m_ref_recomputed_after_subset, FALSE)) {
    .lt_abort("Legacy encoded view may not recompute m_ref after subsetting.")
  }
  if (!identical(x$encoded_identity,
                 .lt_hash(.lt_legacy_encoded_payload(x)))) {
    .lt_abort("Legacy finite-encoded identity is stale.")
  }
  invisible(x)
}

#' Materialize legacy encoded numeric values for diagnostics or replay checks
#'
#' Materialization does not authorize statistical inference. Ordinary La Terra
#' inference routes reject the returned object's finite encoded representation.
#' @param x An `lt_logGBI_encoded` view.
#' @return The derived numeric matrix.
#' @export
lt_encoded_values <- function(x) {
  validate_lt_logGBI_encoded(x)
  x$values
}

#' @export
as.numeric.lt_log_relative <- function(x, ...) {
  .lt_abort_class(
    "lt_log_relative cannot be implicitly coerced to numeric; use an explicit component accessor.",
    "lt_error_log_relative_numeric_coercion"
  )
}

#' @export
as.double.lt_log_relative <- as.numeric.lt_log_relative

#' @export
as.matrix.lt_log_relative <- function(x, ...) {
  .lt_abort_class(
    "lt_log_relative cannot be implicitly coerced to a matrix; use `positive_log_values()` or a state accessor.",
    "lt_error_log_relative_numeric_coercion"
  )
}

#' @export
mean.lt_log_relative <- function(x, ...) {
  .lt_abort_class(
    "A scalar mean is undefined for the boundary object; select ZERO_MASS or POSITIVE_LOG explicitly.",
    "lt_error_log_relative_numeric_coercion"
  )
}

#' @importFrom stats t.test
#' @export
t.test.lt_log_relative <- function(x, ...) {
  .lt_abort_class(
    "t.test() is not an authorized boundary-object association route.",
    "lt_error_log_relative_numeric_coercion"
  )
}

#' @export
as.matrix.lt_logGBI_encoded <- function(x, ...) {
  lt_encoded_values(x)
}

#' @importFrom stats t.test
#' @export
t.test.lt_logGBI_encoded <- function(x, ...) {
  .lt_refuse_legacy_encoded_inference(x, "stats::t.test")
}

#' @export
print.lt_log_relative <- function(x, ...) {
  validate_lt_log_relative(x)
  state <- as.integer(x$ratio_state)
  cat("<lt_log_relative> log_relative_boundary\n")
  cat("  matrix:", length(x$ordered_gene_ledger), "genes x",
      length(x$ordered_branch_ledger), "branches\n")
  cat("  POSITIVE:", sum(state == 1L),
      " EXACT_ZERO_BOUNDARY:", sum(state == 2L),
      " UNAVAILABLE:", sum(state == 3L), "\n")
  cat("  default joint p: not defined\n")
  invisible(x)
}

#' @export
print.lt_logGBI_encoded <- function(x, ...) {
  validate_lt_logGBI_encoded(x)
  cat("<lt_logGBI_encoded>", x$mode, "\n")
  cat("  historical/diagnostic replay view; ordinary inference: REFUSED\n")
  invisible(x)
}
