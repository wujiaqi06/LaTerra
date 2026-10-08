# Rate specifications and C3-derived production representations.

.lt_production_rate_methods <- function() {
  c("ADD_LT", "GBI", "logGBI", "logGBI_strict")
}
.lt_c2_forbidden_rate_methods <- function() c(
  "C2", "ADD_C2", "GBI_C2", "logGBI_C2", "C2_positive_cell_log_OLS"
)

.lt_rate_spec_semantic_payload <- function(x) {
  x[c(
    "rate_id", "method", "parameters", "trait_independent", "input_scale",
    "output_scale", "production_baseline", "execution_status"
  )]
}

#' Construct a La Terra rate specification
#'
#' `ADD_LT`, `GBI`, and the historical positive-only `logGBI_strict` are
#' deterministic views of one certified C3 fitted mean and are always trait
#' independent. An ordinary `logGBI` request returns classed guidance to
#' construct `GBI` and then call [lt_log_relative()]. C2 remains available only
#' through the separate sensitivity API.
#'
#' @param rate_id Stable rate-specification identifier.
#' @param method Named rate representation. The V1 production methods are
#'   `ADD_LT`, `GBI`, and `logGBI_strict`. `logGBI` is retained only for
#'   validation/reading of historical objects and cannot be newly constructed.
#' @param parameters Named parameter list. Production representations accept no
#'   pseudocount, baseline switch, trimming, or trait parameter.
#' @param trait_independent Whether rate construction is trait-independent.
#' @param input_scale Description of the expected input scale.
#' @param output_scale Description of the declared output scale.
#' @param metadata Optional named non-authoritative annotation list. It is not
#'   part of `rate_spec_id`; changing scientific fields still invalidates the
#'   specification identity.
#' @return An object of class `lt_rate_spec`.
#' @export
lt_rate_spec <- function(rate_id,
                         method,
                         parameters = list(),
                         trait_independent = TRUE,
                         input_scale = "branch_associated_values",
                         output_scale = "declared_by_method",
                         metadata = list()) {
  .lt_assert_scalar_character(rate_id, "rate_id")
  .lt_assert_scalar_character(method, "method")
  if (identical(method, "logGBI")) {
    .lt_abort_class(
      paste(
        "Finite-sentinel `logGBI` construction moved out of the production",
        "rate API. Construct GBI, then call `lt_log_relative(gbi)`; use",
        "`logGBI_encoded(..., mode = \"legacy_global_k10\")` only for",
        "explicit historical replay."
      ),
      "lt_error_logGBI_constructor_moved",
      list(requested_method = method)
    )
  }
  if (method %in% .lt_c2_forbidden_rate_methods() ||
      identical(parameters$baseline, "C2")) {
    .lt_abort("C2 is a required sensitivity estimand, not a selectable production rate baseline.")
  }
  if (method %in% .lt_production_rate_methods()) {
    if (!isTRUE(trait_independent)) {
      .lt_abort("Production C3 rate specifications must be trait independent.")
    }
    forbidden <- intersect(
      names(parameters),
      c("pseudocount", "floor", "trim", "winsorize", "trait", "baseline")
    )
    if (length(forbidden)) {
      .lt_abort("Production rate parameters cannot alter baseline, zeros, QC, or traits.")
    }
    output_scale <- if (method == "ADD_LT") "same_physical_units_as_input" else
      "dimensionless"
  }
  core <- list(
    rate_id = rate_id,
    method = method,
    parameters = parameters,
    trait_independent = trait_independent,
    input_scale = input_scale,
    output_scale = output_scale,
    production_baseline = if (method %in% .lt_production_rate_methods())
      "C3_nonnegative_raw_scale_working_mean" else "declarative_only",
    execution_status = if (method %in% .lt_production_rate_methods())
      "T1B_PRODUCTION_RECOGNIZED" else "DECLARATIVE_ONLY_NOT_EXECUTABLE"
  )
  structure(
    c(core, list(metadata = metadata, rate_spec_id = .lt_hash(core))),
    class = "lt_rate_spec"
  )
}

#' Validate a La Terra rate specification
#'
#' Current production specifications are validated for ordinary use. A
#' historical deserialized `method = "logGBI"` specification is
#' read/migrate-only and fails closed here. Enclosing historical `lt_rate`
#' objects use a private structural validator solely so that explicit
#' [lt_migrate_logGBI()] remains possible.
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
.lt_validate_rate_spec <- function(x, allow_historical_logGBI = FALSE) {
  if (!inherits(x, "lt_rate_spec")) {
    .lt_abort("`x` must inherit from `lt_rate_spec`.")
  }
  .lt_assert_scalar_character(x$rate_id, "rate_id")
  .lt_assert_scalar_character(x$method, "method")
  if (identical(x$method, "logGBI") && !isTRUE(allow_historical_logGBI)) {
    .lt_abort_class(
      paste(
        "Historical `lt_rate_spec(method = \"logGBI\")` is read/migrate-only",
        "and cannot enter ordinary validation or execution. Construct GBI and",
        "call `lt_log_relative(gbi)`. `logGBI_encoded(...,",
        "mode = \"legacy_global_k10\")` is diagnostic/replay-only and cannot",
        "enter ordinary inference."
      ),
      "lt_error_logGBI_constructor_moved",
      list(requested_method = "logGBI", source = "historical_rate_spec")
    )
  }
  .lt_assert_named_list(x$parameters, "parameters")
  if (!is.logical(x$trait_independent) || length(x$trait_independent) != 1L ||
      is.na(x$trait_independent)) {
    .lt_abort("`trait_independent` must be TRUE or FALSE.")
  }
  .lt_assert_scalar_character(x$input_scale, "input_scale")
  .lt_assert_scalar_character(x$output_scale, "output_scale")
  .lt_assert_named_list(x$metadata, "metadata")
  if (x$method %in% .lt_c2_forbidden_rate_methods() ||
      identical(x$parameters$baseline, "C2")) {
    .lt_abort("C2 cannot be selected as a production baseline through `lt_rate_spec()`.")
  }
  if (x$method %in% .lt_production_rate_methods()) {
    if (!isTRUE(x$trait_independent) ||
        !identical(x$production_baseline, "C3_nonnegative_raw_scale_working_mean") ||
        !identical(x$execution_status, "T1B_PRODUCTION_RECOGNIZED")) {
      .lt_abort("Production rate specification violates the frozen C3/trait-independent contract.")
    }
  }
  if (!identical(x$rate_spec_id, .lt_hash(.lt_rate_spec_semantic_payload(x)))) {
    .lt_abort("Rate specification hash does not match its normalized fields.")
  }
  invisible(x)
}

#' @rdname validate_lt_rate_spec
#' @export
validate_lt_rate_spec <- function(x) {
  .lt_validate_rate_spec(x, allow_historical_logGBI = FALSE)
}

.lt_is_historical_logGBI_rate <- function(x) {
  inherits(x, "lt_rate") && identical(x$representation, "logGBI")
}

.lt_refuse_historical_logGBI_inference <- function(rate, route) {
  if (.lt_is_historical_logGBI_rate(rate)) {
    .lt_abort_class(
      paste(
        "Historical `lt_rate` representation `logGBI` is read/migrate-only",
        "and is forbidden in ordinary inference. Use `lt_migrate_logGBI(...)`",
        "to obtain `lt_log_relative`. Finite encoded views are retained only",
        "for diagnostic display and source-bound historical replay validation;",
        "use GBI for La Terra V1 ordinary inference."
      ),
      "lt_error_historical_logGBI_inference_forbidden",
      list(route = route, representation = "logGBI")
    )
  }
  invisible(rate)
}

.lt_ratio_baseline_contract <- function(y, mu) {
  if (!is.numeric(y) || !is.numeric(mu) || length(y) != length(mu) ||
      any(!is.finite(y)) || any(!is.finite(mu)) || any(y < 0) || any(mu < 0)) {
    .lt_abort("Ratio baseline contract requires aligned finite nonnegative Y and mu.")
  }
  contradiction <- y > 0 & mu == 0
  if (any(contradiction)) {
    .lt_abort_class(
      "fatal_positive_y_zero_mu: positive observed Y cannot have a zero ratio denominator.",
      "lt_error_positive_y_zero_mu",
      list(cell_indices = which(contradiction))
    )
  }
  list(
    ratio_available = mu > 0,
    reason = ifelse(mu > 0, "available", "baseline_zero")
  )
}

.lt_rate_reason_levels <- function() c(
  "available", "input_zero_log_undefined", "baseline_zero",
  "baseline_unavailable", "numeric_indeterminate", "measurement_ineligible",
  "measurement_unresolved", "coordinate_unavailable",
  "observed_payload_unavailable", "cross_component_baseline_unavailable"
)

.lt_rate_support_levels <- function() c(
  lt_c3_existence_states(), "BASELINE_UNAVAILABLE", "C2_SUPPORTED",
  "C2_UNSUPPORTED_VERTEX", "C2_CROSS_COMPONENT"
)

.lt_rate_keyed_mask_hash <- function(mask, gene_ids, branch_ids) {
  go <- order(gene_ids, method = "radix")
  bo <- order(branch_ids, method = "radix")
  canonical <- mask[go, bo, drop = FALSE]
  dimnames(canonical) <- list(gene_ids[go], branch_ids[bo])
  .lt_hash(list(gene_ids = gene_ids[go], branch_ids = branch_ids[bo],
                mask = canonical))
}

.lt_coordinate_provenance_semantic_payload <- function(provenance) {
  list(
    coordinate_system = provenance$coordinate_system,
    scientific_coordinate_id = provenance$scientific_coordinate_id %||% NA_character_,
    reference_tree_sha256 = provenance$reference_tree_sha256,
    ordered_taxon_ledger_sha256 = provenance$ordered_taxon_ledger_sha256
  )
}

.lt_representation_provenance_semantic_payload <- function(provenance) {
  keep <- intersect(c(
    "source_values_sha256", "source_coordinate_state_sha256",
    "source_value_reason_sha256", "source_matrix_identity_sha256",
    "source_coordinate_provenance_sha256",
    "operator", "baseline_fit_count", "deterministic_view_only",
    "trait_independent", "pseudocount", "numeric_trimming",
    "QC_flags_consulted", "units", "logarithm", "legacy_GBI_equivalence",
    "sensitivity_estimand", "production", "factor_gauge_warning",
    "representation_contract_version", "zero_policy", "log_base", "m_ref",
    "m_ref_source_domain_identity", "m_ref_source_domain_count",
    "L_min_plus", "zero_sentinel", "zero_encoded_count",
    "zero_mask_identity", "source_GBI_identity", "source_GBI_values_sha256",
    "source_C3_fit_identity", "prediction_use_status",
    "conditional_positive_geometry"
  ), names(provenance))
  provenance[keep]
}

.lt_rate_authoritative_hashes <- function(x) {
  status_hash <- .lt_hash(x$measurement_fit_eligibility$status)
  reason_hash <- .lt_hash(x$measurement_fit_eligibility$reason)
  provenance_hash <- .lt_hash(x$measurement_fit_eligibility$provenance)
  base <- list(
    values_sha256 = .lt_hash(x$values),
    coord_state_sha256 = .lt_hash(x$coord_state),
    representation_availability_sha256 = .lt_hash(x$representation_availability),
    representation_reason_sha256 = .lt_hash(x$representation_reason),
    baseline_support_state_sha256 = .lt_hash(x$baseline_support_state),
    measurement_eligibility_sha256 = status_hash,
    measurement_reason_sha256 = reason_hash,
    measurement_provenance_sha256 = provenance_hash,
    ordered_gene_ledger_sha256 = .lt_hash(x$ordered_gene_ledger),
    ordered_branch_ledger_sha256 = .lt_hash(x$ordered_branch_ledger),
    ordered_cell_identity_sha256 = .lt_hash(x$ordered_cell_identity),
    coordinate_provenance_sha256 = .lt_hash(
      .lt_coordinate_provenance_semantic_payload(x$coordinate_provenance)
    ),
    representation_provenance_sha256 = .lt_hash(
      .lt_representation_provenance_semantic_payload(x$representation_provenance)
    ),
    rate_spec_identity_sha256 = if (inherits(x$rate_spec, "lt_rate_spec")) {
      .lt_hash(.lt_rate_spec_semantic_payload(x$rate_spec))
    } else .lt_hash(x$rate_spec),
    representation_reason_vocabulary_sha256 = .lt_hash(x$representation_reason_levels),
    baseline_support_vocabulary_sha256 = .lt_hash(x$baseline_support_state_levels)
  )
  if (!is.null(x$zero_encoded)) {
    base$zero_encoded_sha256 <- .lt_hash(x$zero_encoded)
  }
  if (!is.null(x$zero_origin_state)) {
    base$zero_origin_state_sha256 <- .lt_hash(x$zero_origin_state)
    base$zero_origin_vocabulary_sha256 <- .lt_hash(x$zero_origin_state_levels)
  }
  semantic <- list(
    rate_spec_identity_sha256 = base$rate_spec_identity_sha256,
    representation = x$representation,
    values_sha256 = base$values_sha256,
    coordinate_state_sha256 = base$coord_state_sha256,
    measurement_status_sha256 = status_hash,
    measurement_reason_sha256 = reason_hash,
    measurement_provenance_sha256 = provenance_hash,
    baseline_estimand = x$baseline_estimand,
    baseline_fit_id = x$baseline_fit_id,
    baseline_mu_sha256 = x$baseline_mu_sha256,
    baseline_mu_all_admitted_sha256 = x$baseline_mu_all_admitted_sha256,
    baseline_support_state_sha256 = base$baseline_support_state_sha256,
    representation_availability_sha256 = base$representation_availability_sha256,
    representation_reason_sha256 = base$representation_reason_sha256,
    ordered_gene_ledger_sha256 = base$ordered_gene_ledger_sha256,
    ordered_branch_ledger_sha256 = base$ordered_branch_ledger_sha256,
    ordered_cell_identity_sha256 = base$ordered_cell_identity_sha256,
    coordinate_provenance_sha256 = base$coordinate_provenance_sha256,
    representation_provenance_sha256 = base$representation_provenance_sha256,
    representation_reason_vocabulary_sha256 =
      base$representation_reason_vocabulary_sha256,
    baseline_support_vocabulary_sha256 =
      base$baseline_support_vocabulary_sha256,
    production = isTRUE(x$metadata$production),
    sensitivity = isTRUE(x$metadata$sensitivity),
    sensitivity_estimand = x$metadata$sensitivity_estimand %||% NA_character_
  )
  if (!is.null(x$zero_encoded)) {
    semantic$zero_encoded_sha256 <- base$zero_encoded_sha256
  }
  if (!is.null(x$zero_origin_state)) {
    semantic$zero_origin_state_sha256 <- base$zero_origin_state_sha256
    semantic$zero_origin_vocabulary_sha256 <-
      base$zero_origin_vocabulary_sha256
  }
  base$semantic_object_sha256 <- .lt_hash(semantic)
  base$output_sha256 <- base$semantic_object_sha256
  base
}

.lt_c3_mu_identity <- function(fit) {
  validate_lt_c3_fit(fit)
  if (identical(fit$state, "NUMERICALLY_INDETERMINATE")) {
    .lt_abort("Required C3 component is NUMERICALLY_INDETERMINATE; rate construction fails closed.")
  }
  reconstructed <- rep(NA_real_, length(fit$edge_linear))
  boundary <- fit$baseline_support_state %in% c(
    "ZERO_TOTAL_COMPONENT", "ZERO_MARGIN_BOUNDARY", "FACIAL_BOUNDARY"
  )
  reconstructed[boundary] <- 0
  edge_pos <- integer(length(fit$domain$matrix$values))
  edge_pos[fit$edge_linear] <- seq_along(fit$edge_linear)
  if (length(fit$component_fits) != length(fit$certificates)) {
    .lt_abort("C3 component fits and certificates are not one-to-one.")
  }
  component_hashes <- character(length(fit$component_fits))
  component_ids <- character(length(fit$component_fits))
  for (i in seq_along(fit$component_fits)) {
    z <- fit$component_fits[[i]]
    cert <- fit$certificates[[i]]
    observed_hash <- .lt_hash(z$mu)
    if (!identical(observed_hash, z$mu_sha256) ||
        !identical(observed_hash, cert$mu_sha256) ||
        !identical(cert$mu_sha256, cert$replay_mu_sha256)) {
      .lt_abort("Certified C3 baseline mu hash mismatch.")
    }
    reconstructed[edge_pos[z$global_edge_linear]] <- z$mu
    component_hashes[[i]] <- observed_hash
    component_ids[[i]] <- cert$component_id
  }
  if (anyNA(reconstructed) || !identical(reconstructed, fit$mu)) {
    .lt_abort("Released C3 baseline mu differs from certified component/boundary reconstruction.")
  }
  certified <- if (length(component_hashes) == 1L) component_hashes[[1L]] else
    .lt_hash(data.frame(component_id = component_ids, mu_sha256 = component_hashes,
                        stringsAsFactors = FALSE))
  list(certified_mu_sha256 = certified,
       all_admitted_mu_sha256 = .lt_hash(fit$mu))
}

.lt_rate_identity_gate <- function(x, eligibility, c3_fit) {
  validate_lt_matrix(x)
  validate_lt_c3_fit(c3_fit)
  domain <- c3_fit$domain
  if (!identical(x$gene_ids, domain$matrix$gene_ids)) {
    .lt_abort("Ordered gene ledger differs from the certified C3 fit; explicit keyed mapping is required.")
  }
  if (!identical(x$branch_ids, domain$matrix$branch_ids)) {
    .lt_abort("Ordered branch ledger differs from the certified C3 fit; explicit keyed mapping is required.")
  }
  validate_lt_measurement_eligibility(eligibility, x)
  if (!identical(.lt_hash(x$values), domain$hashes$semantic_values_sha256)) {
    .lt_abort("Input matrix values differ from the exact matrix used for C3 fitting.")
  }
  if (!identical(.lt_hash(x$coord_state),
                 domain$hashes$semantic_coordinate_state_sha256) ||
      !identical(
        .lt_hash(.lt_coordinate_provenance_semantic_payload(x$coordinate_provenance)),
        .lt_hash(.lt_coordinate_provenance_semantic_payload(
          domain$matrix$coordinate_provenance
        ))
      )) {
    .lt_abort("Coordinate state/provenance differs from the certified C3 fit input.")
  }
  if (!identical(eligibility$status_sha256,
                 domain$hashes$measurement_eligibility_sha256) ||
      !identical(eligibility$reason_sha256,
                 domain$hashes$measurement_reason_sha256) ||
      !identical(eligibility$provenance, domain$eligibility$provenance)) {
    .lt_abort("Layer-2 eligibility identity differs from the certified C3 fit.")
  }
  .lt_c3_mu_identity(c3_fit)
}

.lt_rate_compact_codes <- function(labels, levels, dimnames) {
  code <- match(labels, levels)
  if (anyNA(code)) .lt_abort("Internal rate-state code is outside its frozen vocabulary.")
  matrix(as.integer(code), nrow = length(dimnames[[1L]]),
         ncol = length(dimnames[[2L]]), dimnames = dimnames)
}

.lt_new_rate <- function(rate_spec, representation, values, x, eligibility,
                         baseline_estimand, baseline_fit_id,
                         baseline_mu_sha256, baseline_mu_all_admitted_sha256,
                         baseline_support_labels, reason_labels,
                         representation_provenance, metadata = list(),
                         zero_encoded = NULL, zero_origin_state = NULL) {
  dn <- dimnames(x$values)
  representation_provenance <- utils::modifyList(list(
    source_values_sha256 = .lt_hash(x$values),
    source_coordinate_state_sha256 = .lt_hash(x$coord_state),
    source_value_reason_sha256 = .lt_hash(x$value_reason),
    source_matrix_identity_sha256 = .lt_hash(list(
      x$values, x$coord_state, x$value_reason, x$gene_ids, x$branch_ids
    )),
    source_coordinate_provenance_sha256 = .lt_hash(
      .lt_coordinate_provenance_semantic_payload(x$coordinate_provenance)
    )
  ), representation_provenance)
  availability <- !is.na(values)
  reason_code <- .lt_rate_compact_codes(
    reason_labels, .lt_rate_reason_levels(), dn
  )
  support_code <- .lt_rate_compact_codes(
    baseline_support_labels, .lt_rate_support_levels(), dn
  )
  cell_identity <- list(
    gene_ids = x$gene_ids, branch_ids = x$branch_ids,
    full_grid_sha256 = .lt_hash(list(x$gene_ids, x$branch_ids)),
    available_cell_sha256 = .lt_rate_keyed_mask_hash(
      availability, x$gene_ids, x$branch_ids
    )
  )
  fields <- list(
    rate_spec = rate_spec,
    representation = representation,
    values = values,
    coord_state = x$coord_state,
    value_reason = reason_code,
    value_reason_levels = .lt_rate_reason_levels(),
    measurement_fit_eligibility = list(
      status_sha256 = eligibility$status_sha256,
      reason_sha256 = eligibility$reason_sha256,
      status = eligibility$status,
      reason = eligibility$reason,
      provenance = eligibility$provenance,
      provenance_sha256 = .lt_hash(eligibility$provenance)
    ),
    baseline_estimand = baseline_estimand,
    baseline_fit_id = baseline_fit_id,
    baseline_mu_sha256 = baseline_mu_sha256,
    baseline_mu_all_admitted_sha256 = baseline_mu_all_admitted_sha256,
    baseline_support_state = support_code,
    baseline_support_state_levels = .lt_rate_support_levels(),
    representation_availability = availability,
    representation_reason = reason_code,
    representation_reason_levels = .lt_rate_reason_levels(),
    ordered_gene_ledger = x$gene_ids,
    ordered_branch_ledger = x$branch_ids,
    gene_display_labels = .lt_gene_display_labels(x),
    branch_display_labels = .lt_branch_display_labels(x),
    ordered_cell_identity = cell_identity,
    coordinate_provenance = x$coordinate_provenance,
    representation_provenance = representation_provenance,
    output_hashes = NULL,
    metadata = metadata
  )
  if (!is.null(zero_encoded)) {
    if (!is.logical(zero_encoded) || !is.matrix(zero_encoded) ||
        !identical(dimnames(zero_encoded), dn)) {
      .lt_abort("`zero_encoded` must be a logical matrix on the exact rate axes.")
    }
    fields$zero_encoded <- zero_encoded
  }
  if (!is.null(zero_origin_state)) {
    if (!is.raw(zero_origin_state) || !is.matrix(zero_origin_state) ||
        !identical(dimnames(zero_origin_state), dn)) {
      .lt_abort("`zero_origin_state` must be a raw matrix on the exact rate axes.")
    }
    fields$zero_origin_state <- zero_origin_state
    fields$zero_origin_state_levels <- c(
      "NOT_CERTIFIED_OR_UNAVAILABLE",
      "CERTIFIED_POSITIVE_NUMERATOR",
      "CERTIFIED_EXACT_NUMERATOR_ZERO"
    )
  }
  out <- structure(fields, class = "lt_rate")
  out$output_hashes <- .lt_rate_authoritative_hashes(out)
  validate_lt_rate(out)
  out
}

#' Construct a production C3 rate representation
#'
#' C3 is not refit: this operator verifies the exact input, Layer-2, axis, and
#' fitted-mean identities and then materializes one deterministic representation.
#' Exact zeros are retained in GBI together with certified numerator-origin
#' state. Current boundary-aware log-relative objects are then constructed by
#' [lt_log_relative()]. `logGBI_strict` retains historical positive-only
#' conditional geometry. No pseudocount, trimming, QC threshold, trait, or
#' biological result is inspected.
#'
#' @param x Exact `lt_matrix` used to construct the certified C3 fit domain.
#' @param c3_fit Certified `lt_c3_fit`.
#' @param rate_spec Production `lt_rate_spec` for `ADD_LT`, `GBI`, or
#'   `logGBI_strict`. Historical `logGBI` requests fail with guidance.
#' @param measurement_fit_eligibility Exact Layer-2 declaration used by C3.
#' @param metadata Optional named metadata. Only production/sensitivity
#'   contract flags are authoritative; ordinary annotations and display labels
#'   remain outside `semantic_object_sha256`.
#' @return An `lt_rate` object.
#' @export
lt_rate <- function(x, c3_fit, rate_spec,
                    measurement_fit_eligibility = c3_fit$domain$eligibility,
                    metadata = list()) {
  validate_lt_rate_spec(rate_spec)
  if (identical(rate_spec$method, "logGBI")) {
    .lt_abort_class(
      paste(
        "`lt_rate(..., method = \"logGBI\")` is no longer a production",
        "constructor. Construct GBI and call `lt_log_relative(gbi)`."
      ),
      "lt_error_logGBI_constructor_moved",
      list(requested_method = rate_spec$method)
    )
  }
  if (!rate_spec$method %in% .lt_production_rate_methods()) {
    .lt_abort(paste(
      "Only ADD_LT, GBI, logGBI, and logGBI_strict are executable",
      "production rate representations."
    ))
  }
  .lt_assert_named_list(metadata, "metadata")
  mu_identity <- .lt_rate_identity_gate(
    x, measurement_fit_eligibility, c3_fit
  )
  dn <- dimnames(x$values)
  ncell <- length(x$values)
  reason <- rep("coordinate_unavailable", ncell)
  observed <- x$coord_state == "observed"
  reason[observed & measurement_fit_eligibility$status == "ineligible"] <-
    "measurement_ineligible"
  reason[observed & measurement_fit_eligibility$status == "unresolved"] <-
    "measurement_unresolved"
  reason[observed & measurement_fit_eligibility$status == "eligible" &
           is.na(x$values)] <- "observed_payload_unavailable"
  support <- rep("BASELINE_UNAVAILABLE", ncell)
  edge <- c3_fit$edge_linear
  y <- c3_fit$observed_value
  mu <- c3_fit$mu
  support[edge] <- c3_fit$baseline_support_state
  ratio_contract <- .lt_ratio_baseline_contract(y, mu)
  positive_edge <- y > 0 & mu > 0
  encoded_zero_edge <- y == 0 & mu > 0
  available_edge <- switch(
    rate_spec$method,
    ADD_LT = rep(TRUE, length(edge)),
    GBI = ratio_contract$ratio_available,
    logGBI = ratio_contract$ratio_available,
    logGBI_strict = positive_edge
  )
  edge_reason <- rep("available", length(edge))
  edge_reason[!available_edge & mu == 0] <- "baseline_zero"
  edge_reason[!available_edge & mu > 0 & y == 0] <-
    "input_zero_log_undefined"
  reason[edge] <- edge_reason
  out_values <- rep(NA_real_, ncell)
  zero_encoded <- NULL
  zero_origin_state <- NULL
  log_provenance <- list()
  if (rate_spec$method == "ADD_LT") {
    out_values[edge] <- y - mu
  } else if (rate_spec$method == "GBI") {
    out_values[edge[available_edge]] <- y[available_edge] / mu[available_edge]
    zero_origin_state <- matrix(
      as.raw(0L), nrow(x$values), ncol(x$values), dimnames = dn
    )
    zero_origin_state[edge[positive_edge]] <- as.raw(1L)
    zero_origin_state[edge[encoded_zero_edge]] <- as.raw(2L)
  } else if (rate_spec$method == "logGBI_strict") {
    out_values[edge[available_edge]] <- log(y[available_edge] / mu[available_edge])
    zero_encoded <- matrix(NA, nrow(x$values), ncol(x$values), dimnames = dn)
    zero_encoded[edge[available_edge]] <- FALSE
    strict_zero_mask <- matrix(
      FALSE, nrow(x$values), ncol(x$values), dimnames = dn
    )
    log_provenance <- list(
      representation_contract_version = "LaTerra_V1_logGBI_strict_v1",
      zero_policy = "strict_positive_only_unavailable_at_zero",
      log_base = "natural",
      zero_encoded_count = 0L,
      zero_mask_identity = .lt_rate_keyed_mask_hash(
        strict_zero_mask, x$gene_ids, x$branch_ids
      ),
      conditional_positive_geometry = TRUE,
      prediction_use_status = "NOT_APPLICABLE_STRICT_POSITIVE_ONLY"
    )
  } else {
    if (!any(positive_edge)) {
      .lt_abort("Default logGBI requires at least one positive GBI reference cell.")
    }
    positive_gbi <- y[positive_edge] / mu[positive_edge]
    m_ref <- min(positive_gbi)
    L_min_plus <- log(m_ref)
    zero_sentinel <- log(m_ref) - 10
    out_values[edge[positive_edge]] <- log(y[positive_edge] / mu[positive_edge])
    out_values[edge[encoded_zero_edge]] <- zero_sentinel
    zero_encoded <- matrix(NA, nrow(x$values), ncol(x$values), dimnames = dn)
    zero_encoded[edge[available_edge]] <- FALSE
    zero_encoded[edge[encoded_zero_edge]] <- TRUE
    positive_mask <- matrix(FALSE, nrow(x$values), ncol(x$values), dimnames = dn)
    positive_mask[edge[positive_edge]] <- TRUE
    zero_mask <- matrix(FALSE, nrow(x$values), ncol(x$values), dimnames = dn)
    zero_mask[edge[encoded_zero_edge]] <- TRUE
    gbi_values <- rep(NA_real_, ncell)
    gbi_values[edge[mu > 0]] <- y[mu > 0] / mu[mu > 0]
    dim(gbi_values) <- dim(x$values); dimnames(gbi_values) <- dn
    gbi_mask <- !is.na(gbi_values)
    source_gbi_values_sha256 <- .lt_hash(gbi_values)
    source_gbi_domain_identity <- .lt_rate_keyed_mask_hash(
      gbi_mask, x$gene_ids, x$branch_ids
    )
    source_gbi_identity <- .lt_hash(list(
      representation = "GBI",
      values_sha256 = source_gbi_values_sha256,
      available_cell_sha256 = source_gbi_domain_identity,
      baseline_mu_sha256 = mu_identity$certified_mu_sha256
    ))
    source_c3_fit_identity <- .lt_hash(list(
      run_id = c3_fit$run_id,
      baseline_mu_sha256 = mu_identity$certified_mu_sha256,
      baseline_mu_all_admitted_sha256 = mu_identity$all_admitted_mu_sha256,
      source_values_sha256 = c3_fit$domain$hashes$semantic_values_sha256,
      source_coordinate_state_sha256 =
        c3_fit$domain$hashes$semantic_coordinate_state_sha256,
      measurement_eligibility_sha256 =
        c3_fit$domain$hashes$measurement_eligibility_sha256,
      ordered_cell_ledger_sha256 =
        c3_fit$domain$hashes$ordered_cell_ledger_sha256
    ))
    log_provenance <- list(
      representation_contract_version = "LaTerra_V1_logGBI_zero_state_v1",
      zero_policy = "global_positive_log_min_minus_10",
      log_base = "natural",
      m_ref = m_ref,
      m_ref_source_domain_identity = .lt_rate_keyed_mask_hash(
        positive_mask, x$gene_ids, x$branch_ids
      ),
      m_ref_source_domain_count = as.integer(sum(positive_mask)),
      L_min_plus = L_min_plus,
      zero_sentinel = zero_sentinel,
      zero_encoded_count = as.integer(sum(zero_mask)),
      zero_mask_identity = .lt_rate_keyed_mask_hash(
        zero_mask, x$gene_ids, x$branch_ids
      ),
      source_GBI_identity = source_gbi_identity,
      source_GBI_values_sha256 = source_gbi_values_sha256,
      source_C3_fit_identity = source_c3_fit_identity,
      prediction_use_status = paste0(
        "ASSOCIATION_OK_TRAINING_FREEZE_REQUIRED_FOR_",
        "HELDOUT_PREDICTION"
      ),
      conditional_positive_geometry = FALSE
    )
  }
  if (any(is.infinite(out_values), na.rm = TRUE) || any(is.nan(out_values))) {
    .lt_abort("Rate representation cannot store Inf, -Inf, or NaN.")
  }
  dim(out_values) <- dim(x$values); dimnames(out_values) <- dn
  dim(reason) <- dim(x$values); dimnames(reason) <- dn
  dim(support) <- dim(x$values); dimnames(support) <- dn
  units <- if (rate_spec$method == "ADD_LT") x$payload$units else "dimensionless"
  provenance <- utils::modifyList(list(
    operator = paste0("LaTerra_T1B_C3_", rate_spec$method),
    baseline_fit_count = 1L,
    deterministic_view_only = TRUE,
    trait_independent = TRUE,
    pseudocount = "none",
    numeric_trimming = FALSE,
    QC_flags_consulted = FALSE,
    units = units,
    logarithm = if (rate_spec$method %in% c("logGBI", "logGBI_strict"))
      "natural" else "not_applicable",
    legacy_GBI_equivalence = FALSE
  ), log_provenance)
  .lt_new_rate(
    rate_spec = rate_spec,
    representation = rate_spec$method,
    values = out_values,
    x = x,
    eligibility = measurement_fit_eligibility,
    baseline_estimand = "C3_nonnegative_raw_scale_working_mean",
    baseline_fit_id = c3_fit$run_id,
    baseline_mu_sha256 = mu_identity$certified_mu_sha256,
    baseline_mu_all_admitted_sha256 = mu_identity$all_admitted_mu_sha256,
    baseline_support_labels = support,
    reason_labels = reason,
    representation_provenance = provenance,
    metadata = utils::modifyList(list(
      units = units, production = TRUE, sensitivity = FALSE,
      scientific_pass = FALSE,
      representation_semantics = if (rate_spec$method == "logGBI")
        "positive_log_relative_magnitude_plus_explicit_zero_state" else
        if (rate_spec$method == "logGBI_strict")
          "conditional_positive_log_relative_geometry" else
          "direct_C3_deterministic_view"
    ), metadata),
    zero_encoded = zero_encoded
    , zero_origin_state = zero_origin_state
  )
}

.lt_rate_log_contract <- function(x) {
  if (identical(x$representation, "logGBI") &&
      is.null(x$zero_encoded) &&
      identical(x$representation_provenance$operator,
                "LaTerra_T1B_C3_logGBI") &&
      is.null(x$representation_provenance$representation_contract_version)) {
    return("LEGACY_DEVELOPMENT_LOGGBI_STRICT_NOT_REINTERPRETED")
  }
  if (identical(x$representation, "logGBI")) {
    return("DEFAULT_ZERO_STATE")
  }
  if (identical(x$representation, "logGBI_strict")) {
    return("STRICT_POSITIVE_ONLY")
  }
  "NOT_LOGGBI"
}

#' Validate a production or C2-sensitivity rate object
#'
#' Validation recomputes the semantic identity from all authoritative fields,
#' including embedded Layer-2 matrices and provenance. Display labels and
#' ordinary non-authoritative metadata are deliberately excluded. Use
#' [lt_audit_fingerprint()] when exact frozen-snapshot identity is required.
#'
#' @param x An `lt_rate`.
#' @return `x`, invisibly.
#' @export
validate_lt_rate <- function(x) {
  if (!inherits(x, "lt_rate")) .lt_abort("`x` must inherit from lt_rate.")
  .lt_assert_scalar_character(x$representation, "representation")
  .lt_assert_scalar_character(x$baseline_estimand, "baseline_estimand")
  .lt_assert_scalar_character(x$baseline_fit_id, "baseline_fit_id")
  .lt_assert_sha256(x$baseline_mu_sha256, "baseline_mu_sha256")
  .lt_assert_sha256(x$baseline_mu_all_admitted_sha256,
                    "baseline_mu_all_admitted_sha256")
  dn <- list(x$ordered_gene_ledger, x$ordered_branch_ledger)
  required_matrix <- c(
    "values", "coord_state", "value_reason", "baseline_support_state",
    "representation_availability", "representation_reason"
  )
  if (any(!vapply(x[required_matrix], is.matrix, logical(1L))) ||
      any(!vapply(x[required_matrix], function(z) identical(dimnames(z), dn), logical(1L)))) {
    .lt_abort("Rate value/state layers must be matrices on the exact ordered axes.")
  }
  .lt_assert_unique_ids(x$ordered_gene_ledger, "ordered_gene_ledger")
  .lt_assert_unique_ids(x$ordered_branch_ledger, "ordered_branch_ledger")
  gene_display <- x$gene_display_labels %||% x$ordered_gene_ledger
  branch_display <- x$branch_display_labels %||% x$ordered_branch_ledger
  if (!is.character(gene_display) || length(gene_display) != length(dn[[1L]]) ||
      anyNA(gene_display) || any(!nzchar(gene_display)) ||
      !is.character(branch_display) || length(branch_display) != length(dn[[2L]]) ||
      anyNA(branch_display) || any(!nzchar(branch_display))) {
    .lt_abort("Rate display labels must align to scientific keys; duplicates are permitted.")
  }
  if (!is.numeric(x$values) || any(is.infinite(x$values), na.rm = TRUE) ||
      any(is.nan(x$values))) {
    .lt_abort("Rate values must be numeric without Inf, -Inf, or NaN.")
  }
  if (!is.logical(x$representation_availability) ||
      anyNA(x$representation_availability) ||
      !identical(x$representation_availability, !is.na(x$values))) {
    .lt_abort("Rate availability must exactly equal the finite materialized-value domain.")
  }
  log_contract <- .lt_rate_log_contract(x)
  if (!is.null(x$zero_encoded)) {
    if (!is.logical(x$zero_encoded) || !is.matrix(x$zero_encoded) ||
        !identical(dimnames(x$zero_encoded), dn) ||
        !identical(dim(x$zero_encoded), dim(x$values)) ||
        any(is.na(x$zero_encoded[x$representation_availability])) ||
        any(!is.na(x$zero_encoded[!x$representation_availability]))) {
      .lt_abort(paste(
        "zero_encoded must be TRUE/FALSE on available cells and NA on",
        "unavailable cells, on the exact rate axes."
      ))
    }
  }
  if (!is.null(x$zero_origin_state)) {
    expected_levels <- c(
      "NOT_CERTIFIED_OR_UNAVAILABLE",
      "CERTIFIED_POSITIVE_NUMERATOR",
      "CERTIFIED_EXACT_NUMERATOR_ZERO"
    )
    if (!identical(x$representation, "GBI") ||
        !is.raw(x$zero_origin_state) || !is.matrix(x$zero_origin_state) ||
        !identical(dimnames(x$zero_origin_state), dn) ||
        any(as.integer(x$zero_origin_state) > 2L) ||
        !identical(x$zero_origin_state_levels, expected_levels)) {
      .lt_abort("Certified GBI zero-origin state is malformed or attached to a non-GBI rate.")
    }
    origin <- as.integer(x$zero_origin_state)
    available <- as.vector(x$representation_availability)
    values <- as.vector(x$values)
    if (any(origin == 2L & (!available | values != 0)) ||
        any(origin == 1L & !available)) {
      .lt_abort("Certified GBI zero-origin state conflicts with GBI availability/value storage.")
    }
  }
  if (!is.integer(x$representation_reason) ||
      anyNA(x$representation_reason) ||
      any(!x$representation_reason %in% seq_along(x$representation_reason_levels)) ||
      !identical(x$value_reason, x$representation_reason)) {
    .lt_abort("Compact representation reason codes violate their frozen vocabulary.")
  }
  reason <- x$representation_reason_levels[x$representation_reason]
  if (any(x$representation_availability & reason != "available") ||
      any(!x$representation_availability & reason == "available")) {
    .lt_abort("Rate availability and representation reason disagree.")
  }
  if (!is.integer(x$baseline_support_state) ||
      anyNA(x$baseline_support_state) ||
      any(!x$baseline_support_state %in% seq_along(x$baseline_support_state_levels))) {
    .lt_abort("Compact baseline support codes violate their frozen vocabulary.")
  }

  eligibility <- x$measurement_fit_eligibility
  if (!is.list(eligibility) || !is.matrix(eligibility$status) ||
      !is.character(eligibility$status) || !is.matrix(eligibility$reason) ||
      !is.character(eligibility$reason) ||
      !identical(dimnames(eligibility$status), dn) ||
      !identical(dimnames(eligibility$reason), dn) ||
      !identical(dim(eligibility$status), dim(x$values)) ||
      !identical(dim(eligibility$reason), dim(x$values))) {
    .lt_abort("Embedded Layer-2 status/reason must be character matrices on the exact rate axes.")
  }
  .lt_assert_named_list(
    eligibility$provenance, "measurement_fit_eligibility$provenance"
  )
  actual_status_sha256 <- .lt_hash(eligibility$status)
  actual_reason_sha256 <- .lt_hash(eligibility$reason)
  actual_provenance_sha256 <- .lt_hash(eligibility$provenance)
  if (!identical(eligibility$status_sha256, actual_status_sha256) ||
      !identical(eligibility$reason_sha256, actual_reason_sha256) ||
      !identical(eligibility$provenance_sha256, actual_provenance_sha256)) {
    .lt_abort("Embedded Layer-2 status/reason/provenance identity is stale or inconsistent.")
  }
  if (anyNA(eligibility$status) || anyNA(eligibility$reason) ||
      any(!nzchar(eligibility$status)) || any(!nzchar(eligibility$reason))) {
    .lt_abort("Embedded Layer-2 status/reason labels must be explicit.")
  }
  if (length(setdiff(unique(as.vector(eligibility$status)),
                     .lt_eligibility_states())) ||
      length(setdiff(unique(as.vector(eligibility$reason)),
                     .lt_eligibility_reasons()))) {
    .lt_abort("Embedded Layer-2 status/reason uses an unknown frozen label.")
  }
  observed_coordinate <- x$coord_state == "observed"
  if (any(!observed_coordinate & eligibility$status != "not_applicable") ||
      any(observed_coordinate & eligibility$status == "not_applicable") ||
      any(eligibility$status == "not_applicable" &
            eligibility$reason != "coordinate_not_observed") ||
      any(eligibility$status == "eligible" &
            eligibility$reason != "eligible_declared") ||
      any(eligibility$status == "unresolved" &
            eligibility$reason != "authority_unresolved") ||
      any(eligibility$status == "ineligible" & eligibility$reason %in% c(
        "eligible_declared", "authority_unresolved", "coordinate_not_observed"
      ))) {
    .lt_abort("Embedded Layer-2 coordinate-scoped state machine is inconsistent.")
  }

  if (x$representation %in% .lt_production_rate_methods()) {
    .lt_validate_rate_spec(
      x$rate_spec,
      allow_historical_logGBI = identical(x$representation, "logGBI")
    )
    if (!identical(x$rate_spec$method, x$representation) ||
        !identical(x$baseline_estimand,
                   "C3_nonnegative_raw_scale_working_mean") ||
        !isTRUE(x$metadata$production) || isTRUE(x$metadata$sensitivity)) {
      .lt_abort("Production lt_rate cannot masquerade as another representation or estimand.")
    }
  } else if (x$representation %in% c("ADD_C2", "GBI_C2", "logGBI_C2")) {
    if (!inherits(x$rate_spec, "lt_c2_sensitivity_spec") ||
        !identical(x$rate_spec$method, x$representation) ||
        !identical(x$baseline_estimand, "C2_positive_cell_log_OLS") ||
        !identical(x$rate_spec$sensitivity_estimand,
                   "C2_positive_cell_log_OLS") ||
        !identical(x$metadata$sensitivity_estimand,
                   "C2_positive_cell_log_OLS") ||
        isTRUE(x$metadata$production) || !isTRUE(x$metadata$sensitivity)) {
      .lt_abort("C2 sensitivity lt_rate cannot masquerade as a production C3 object.")
    }
  } else {
    .lt_abort("Unknown lt_rate representation contract.")
  }

  provenance <- x$representation_provenance
  if (identical(log_contract, "DEFAULT_ZERO_STATE")) {
    required <- c(
      "representation_contract_version", "zero_policy", "log_base", "m_ref",
      "m_ref_source_domain_identity", "m_ref_source_domain_count",
      "L_min_plus", "zero_sentinel", "zero_encoded_count",
      "zero_mask_identity", "source_GBI_identity", "source_GBI_values_sha256",
      "source_C3_fit_identity", "prediction_use_status"
    )
    if (is.null(x$zero_encoded) || any(!required %in% names(provenance)) ||
        !identical(provenance$representation_contract_version,
                   "LaTerra_V1_logGBI_zero_state_v1") ||
        !identical(provenance$zero_policy,
                   "global_positive_log_min_minus_10") ||
        !identical(provenance$log_base, "natural") ||
        !is.numeric(provenance$m_ref) || length(provenance$m_ref) != 1L ||
        !is.finite(provenance$m_ref) || provenance$m_ref <= 0 ||
        !identical(provenance$L_min_plus, log(provenance$m_ref)) ||
        !identical(provenance$zero_sentinel, log(provenance$m_ref) - 10) ||
        !identical(provenance$zero_encoded_count,
                   as.integer(sum(x$zero_encoded, na.rm = TRUE))) ||
        !identical(provenance$m_ref_source_domain_count,
                   as.integer(sum(x$representation_availability &
                                    !x$zero_encoded, na.rm = TRUE))) ||
        !identical(provenance$prediction_use_status,
          "ASSOCIATION_OK_TRAINING_FREEZE_REQUIRED_FOR_HELDOUT_PREDICTION")) {
      .lt_abort("Default logGBI zero-state provenance contract is incomplete or inconsistent.")
    }
    zero_mask <- x$zero_encoded %in% TRUE
    dim(zero_mask) <- dim(x$values); dimnames(zero_mask) <- dn
    positive_mask <- x$representation_availability & !zero_mask
    if (!identical(provenance$zero_mask_identity,
                   .lt_rate_keyed_mask_hash(
                     zero_mask, x$ordered_gene_ledger,
                     x$ordered_branch_ledger
                   )) ||
        !identical(provenance$m_ref_source_domain_identity,
                   .lt_rate_keyed_mask_hash(
                     positive_mask, x$ordered_gene_ledger,
                     x$ordered_branch_ledger
                   )) ||
        !is.character(provenance$source_GBI_identity) ||
        length(provenance$source_GBI_identity) != 1L ||
        !grepl("^[0-9a-f]{64}$", provenance$source_GBI_identity) ||
        !is.character(provenance$source_GBI_values_sha256) ||
        length(provenance$source_GBI_values_sha256) != 1L ||
        !grepl("^[0-9a-f]{64}$", provenance$source_GBI_values_sha256) ||
        !is.character(provenance$source_C3_fit_identity) ||
        length(provenance$source_C3_fit_identity) != 1L ||
        !grepl("^[0-9a-f]{64}$", provenance$source_C3_fit_identity)) {
      .lt_abort("Default logGBI zero/positive domain identity contract failed.")
    }
    if (any(zero_mask) &&
        (!identical(unique(x$values[zero_mask]), provenance$zero_sentinel) ||
         !(provenance$zero_sentinel < min(x$values[positive_mask])))) {
      .lt_abort("Default logGBI must use one global sentinel below every positive cell.")
    }
    expected_gbi_identity <- .lt_hash(list(
      representation = "GBI",
      values_sha256 = provenance$source_GBI_values_sha256,
      available_cell_sha256 = .lt_rate_keyed_mask_hash(
        x$representation_availability, x$ordered_gene_ledger,
        x$ordered_branch_ledger
      ),
      baseline_mu_sha256 = x$baseline_mu_sha256
    ))
    if (!any(positive_mask) ||
        !identical(provenance$L_min_plus, min(x$values[positive_mask])) ||
        !identical(provenance$source_GBI_identity, expected_gbi_identity)) {
      .lt_abort("Default logGBI reference minimum/source GBI identity is inconsistent.")
    }
  } else if (identical(log_contract, "STRICT_POSITIVE_ONLY")) {
    if (is.null(x$zero_encoded) || any(x$zero_encoded, na.rm = TRUE) ||
        !identical(provenance$representation_contract_version,
                   "LaTerra_V1_logGBI_strict_v1") ||
        !identical(provenance$zero_policy,
                   "strict_positive_only_unavailable_at_zero") ||
        !identical(provenance$log_base, "natural") ||
        !identical(provenance$zero_encoded_count, 0L) ||
        !identical(provenance$zero_mask_identity,
                   .lt_rate_keyed_mask_hash(
                     matrix(FALSE, nrow(x$values), ncol(x$values),
                            dimnames = dn),
                     x$ordered_gene_ledger, x$ordered_branch_ledger
                   )) ||
        !identical(provenance$conditional_positive_geometry, TRUE)) {
      .lt_abort("logGBI_strict positive-only companion contract is inconsistent.")
    }
  } else if (identical(
    log_contract, "LEGACY_DEVELOPMENT_LOGGBI_STRICT_NOT_REINTERPRETED"
  )) {
    if (any(x$representation_availability & reason != "available") ||
        any(!x$representation_availability & reason == "available")) {
      .lt_abort("Legacy strict logGBI availability/reason contract failed.")
    }
  } else if (!is.null(x$zero_encoded)) {
    .lt_abort("zero_encoded is defined only for logGBI representations.")
  }

  observed_hashes <- .lt_rate_authoritative_hashes(x)
  if (!identical(x$output_hashes, observed_hashes)) {
    .lt_abort("Rate semantic identity is stale or inconsistent with authoritative fields.")
  }
  invisible(x)
}

#' Materialize representation reasons or an identity-keyed rate ledger
#'
#' `lt_rate_reason()` expands the compact reason codes as an aligned character
#' matrix. `lt_rate_ledger()` exposes one keyed row per cell; callers should use
#' it selectively for large matrices because the compact `lt_rate` object is the
#' primary memory-safe representation.
#'
#' @param x An `lt_rate`.
#' @param available_only Whether to materialize only available cells.
#' @return A character matrix or data frame.
#' @export
lt_rate_reason <- function(x) {
  validate_lt_rate(x)
  matrix(
    x$representation_reason_levels[x$representation_reason],
    nrow = nrow(x$values), ncol = ncol(x$values),
    dimnames = dimnames(x$values)
  )
}

#' @rdname lt_rate_reason
#' @export
lt_rate_ledger <- function(x, available_only = FALSE) {
  validate_lt_rate(x)
  if (!is.logical(available_only) || length(available_only) != 1L ||
      is.na(available_only)) .lt_abort("`available_only` must be TRUE or FALSE.")
  linear <- if (available_only) which(x$representation_availability) else
    seq_along(x$values)
  nr <- nrow(x$values)
  gi <- ((linear - 1L) %% nr) + 1L
  bi <- ((linear - 1L) %/% nr) + 1L
  data.frame(
    gene_id = x$ordered_gene_ledger[gi],
    branch_id = x$ordered_branch_ledger[bi],
    value = x$values[linear],
    coord_state = x$coord_state[linear],
    measurement_status = x$measurement_fit_eligibility$status[linear],
    measurement_reason = x$measurement_fit_eligibility$reason[linear],
    baseline_support_state = x$baseline_support_state_levels[
      x$baseline_support_state[linear]
    ],
    representation_available = x$representation_availability[linear],
    representation_reason = x$representation_reason_levels[
      x$representation_reason[linear]
    ],
    zero_encoded = if (is.null(x$zero_encoded))
      rep(NA, length(linear)) else x$zero_encoded[linear],
    stringsAsFactors = FALSE
  )
}

#' @export
print.lt_rate <- function(x, ...) {
  validate_lt_rate(x)
  cat("<lt_rate> ", x$representation, "\n", sep = "")
  cat("  baseline: ", x$baseline_estimand, "\n", sep = "")
  cat("  matrix: ", nrow(x$values), " genes x ", ncol(x$values), " branches\n", sep = "")
  cat("  available cells: ", sum(x$representation_availability), "\n", sep = "")
  if (.lt_rate_log_contract(x) != "NOT_LOGGBI") {
    cat("  log contract: ", .lt_rate_log_contract(x), "\n", sep = "")
    if (!is.null(x$zero_encoded)) {
      cat("  zero-encoded cells: ", sum(x$zero_encoded, na.rm = TRUE),
          "\n", sep = "")
    }
  }
  cat("  production: ", isTRUE(x$metadata$production), "\n", sep = "")
  invisible(x)
}
