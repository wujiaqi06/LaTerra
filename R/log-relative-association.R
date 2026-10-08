# Boundary-aware binary association, external-null calibration, and adjustment.

.lt_boundary_component_names <- function() c("ZERO_MASS", "POSITIVE_LOG")
.lt_boundary_engine_token <- new.env(parent = emptyenv())

.lt_boundary_state_index <- function(object, branch_state, branch_ids,
                                     mismatch) {
  validate_lt_log_relative(object)
  validate_lt_branch_state(branch_state)
  if (!identical(branch_state$type, "binary")) {
    .lt_abort_class(
      "Boundary-aware two-component association is defined only for binary branch state in V1.",
      "lt_error_log_relative_continuous_trait_out_of_scope"
    )
  }
  mismatch <- match.arg(mismatch, c("warn", "error"))
  rate_ids <- object$ordered_branch_ledger
  rate_only <- setdiff(rate_ids, branch_state$branch_ids)
  state_only <- setdiff(branch_state$branch_ids, rate_ids)
  if (length(rate_only) || length(state_only)) {
    msg <- paste0(
      "Exact branch-key mismatch: ", length(rate_only), " object-only and ",
      length(state_only), " state-only branches."
    )
    if (mismatch == "error") .lt_abort(msg)
    warning(msg, call. = FALSE)
  }
  index <- match(rate_ids, branch_state$branch_ids)
  matched <- !is.na(index)
  available <- rep(FALSE, length(rate_ids))
  available[matched] <-
    branch_state$availability[index[matched]] == "available"
  user_allowed <- rep(TRUE, length(rate_ids))
  if (!is.null(branch_ids)) {
    .lt_assert_unique_ids(branch_ids, "branch_ids")
    unknown <- setdiff(branch_ids, rate_ids)
    if (length(unknown)) {
      .lt_abort("Requested branch_ids contain unknown scientific keys.")
    }
    user_allowed <- rate_ids %in% branch_ids
  }
  state_values <- rep(NA, length(rate_ids))
  state_values[matched] <- branch_state$values[index[matched]]
  list(
    index = as.integer(index),
    state_values = state_values,
    branch_allowed = available & user_allowed,
    user_allowed = user_allowed,
    rate_only = rate_only,
    state_only = state_only,
    mismatch = mismatch
  )
}

.lt_boundary_component <- function(component, statistic, status, reason,
                                   counts, object, branch_state,
                                   association_id) {
  core <- list(
    component = component,
    estimand = if (component == "ZERO_MASS")
      "delta_z_focal_minus_reference" else
      "delta_L_positive_focal_minus_reference",
    gene_ids = object$ordered_gene_ledger,
    statistic = statistic,
    status = status,
    reason = reason,
    counts = counts,
    source_semantic_identity = object$semantic_identity,
    source_audit_identity = object$audit_identity,
    branch_state_identity = branch_state$branch_state_identity,
    association_id = association_id
  )
  structure(c(core, list(component_identity = .lt_hash(core))),
            class = "lt_boundary_component_association")
}

.lt_boundary_score <- function(object, state_values, branch_allowed,
                               branch_state, .validated_token = NULL) {
  if (!identical(.validated_token, .lt_boundary_engine_token)) {
    .lt_abort_class(
      "Direct boundary-engine bypass is forbidden; use a validated public association/calibration route.",
      "lt_error_log_relative_internal_bypass"
    )
  }
  ref_key <- .lt_level_key(branch_state$coding$reference_level)
  focal_key <- .lt_level_key(branch_state$coding$focal_level)
  keys <- rep(NA_character_, length(state_values))
  present <- !is.na(state_values)
  keys[present] <- vapply(
    state_values[present], .lt_level_key, character(1L)
  )
  ref_cols <- which(keys == ref_key & branch_allowed)
  focal_cols <- which(keys == focal_key & branch_allowed)
  state <- as.integer(object$ratio_state)
  dim(state) <- dim(object$ratio_state)
  dimnames(state) <- dimnames(object$ratio_state)
  ratio_ok <- state %in% c(1L, 2L)
  positive <- state == 1L
  zero <- state == 2L
  dim(ratio_ok) <- dim(state); dimnames(ratio_ok) <- dimnames(state)
  dim(positive) <- dim(state); dimnames(positive) <- dimnames(state)
  dim(zero) <- dim(state); dimnames(zero) <- dimnames(state)
  logs <- object$positive_log_values
  logs[!positive] <- 0
  summarize <- function(columns) {
    ratio_part <- ratio_ok[, columns, drop = FALSE]
    positive_part <- positive[, columns, drop = FALSE]
    zero_part <- zero[, columns, drop = FALSE]
    list(
      ratio_n = as.integer(rowSums(ratio_part)),
      zero_n = as.integer(rowSums(zero_part)),
      positive_n = as.integer(rowSums(positive_part)),
      positive_sum = rowSums(logs[, columns, drop = FALSE])
    )
  }
  ref <- summarize(ref_cols)
  focal <- summarize(focal_cols)
  dz <- focal$zero_n / focal$ratio_n - ref$zero_n / ref$ratio_n
  dl <- focal$positive_sum / focal$positive_n -
    ref$positive_sum / ref$positive_n
  zero_valid <- ref$ratio_n > 0L & focal$ratio_n > 0L & is.finite(dz)
  positive_valid <- ref$positive_n > 0L & focal$positive_n > 0L &
    is.finite(dl)
  dz[!zero_valid] <- NA_real_
  dl[!positive_valid] <- NA_real_
  dz <- unname(dz)
  dl <- unname(dl)
  list(
    ZERO_MASS = list(
      statistic = dz,
      status = unname(ifelse(zero_valid, "ESTIMABLE", "NONESTIMABLE")),
      reason = unname(ifelse(
        zero_valid, "estimable_zero_fraction_contrast",
        "focal_or_reference_group_absent"
      )),
      counts = data.frame(
        gene_id = object$ordered_gene_ledger,
        n_reference = ref$ratio_n, n_focal = focal$ratio_n,
        reference_exact_zero_count = ref$zero_n,
        focal_exact_zero_count = focal$zero_n,
        stringsAsFactors = FALSE
      )
    ),
    POSITIVE_LOG = list(
      statistic = dl,
      status = unname(ifelse(positive_valid, "ESTIMABLE", "NONESTIMABLE")),
      reason = unname(ifelse(
        positive_valid, "estimable_conditional_positive_contrast",
        "focal_or_reference_positive_support_absent"
      )),
      counts = data.frame(
        gene_id = object$ordered_gene_ledger,
        n_reference = ref$positive_n, n_focal = focal$positive_n,
        reference_exact_zero_count = ref$zero_n,
        focal_exact_zero_count = focal$zero_n,
        stringsAsFactors = FALSE
      )
    )
  )
}

.lt_boundary_association_payload <- function(x) {
  x[c(
    "association_id", "source_semantic_identity", "source_audit_identity",
    "branch_state_identity", "trait_id", "gene_ids", "branch_ids",
    "state_index", "branch_allowed", "components",
    "joint_inference_status", "provenance"
  )]
}

#' Associate a boundary-aware log-relative object with a binary branch state
#'
#' Returns separate ZERO_MASS and POSITIVE_LOG estimands. No joint statistic or
#' joint p-value is defined.
#'
#' @param object An `lt_log_relative`.
#' @param branch_state A binary `lt_branch_state`.
#' @param branch_ids Optional exact scientific-key restriction.
#' @param mismatch Exact-key mismatch policy.
#' @param association_id Stable result identifier.
#' @param null_ensembles Optional named list of frozen Null A/B ensembles.
#' @param adjust Optional separate adjustment methods, normally `BH` and `BY`.
#' @return An `lt_boundary_association`, or an `lt_log_relative_result` when
#'   null ensembles are supplied.
#' @export
lt_associate_log_relative <- function(
    object, branch_state, branch_ids = NULL, mismatch = c("warn", "error"),
    association_id = paste0("log_relative_boundary_",
                            branch_state$state_id, "_association"),
    null_ensembles = NULL, adjust = c("BH", "BY")) {
  validate_lt_log_relative(object)
  joined <- .lt_boundary_state_index(
    object, branch_state, branch_ids, mismatch
  )
  scored <- .lt_boundary_score(
    object, joined$state_values, joined$branch_allowed, branch_state,
    .validated_token = .lt_boundary_engine_token
  )
  components <- lapply(.lt_boundary_component_names(), function(name) {
    z <- scored[[name]]
    .lt_boundary_component(
      name, z$statistic, z$status, z$reason, z$counts,
      object, branch_state, association_id
    )
  })
  names(components) <- .lt_boundary_component_names()
  core <- list(
    association_id = association_id,
    source_semantic_identity = object$semantic_identity,
    source_audit_identity = object$audit_identity,
    branch_state_identity = branch_state$branch_state_identity,
    trait_id = branch_state$trait_id,
    gene_ids = object$ordered_gene_ledger,
    branch_ids = object$ordered_branch_ledger,
    state_index = joined$index,
    branch_allowed = joined$branch_allowed,
    components = components,
    joint_inference_status = "NOT_DEFINED_BY_CONTRACT",
    provenance = list(
      operator = "LaTerra_boundary_aware_two_component_binary_association",
      exact_key_join = TRUE,
      sentinel_used = FALSE,
      imputation = FALSE,
      component_gene_sets_forced_equal = FALSE,
      continuous_trait_analogue = "OUT_OF_SCOPE"
    )
  )
  association <- structure(c(core, list(
    association_identity = .lt_hash(core), metadata = list()
  )), class = "lt_boundary_association")
  validate_lt_boundary_association(association)
  if (is.null(null_ensembles)) return(association)
  calibration <- lt_calibrate_log_relative(
    object, branch_state, association, null_ensembles
  )
  adjustments <- lt_adjust_log_relative(calibration, methods = adjust)
  result_core <- list(
    association = association,
    calibration = calibration,
    adjustments = adjustments,
    joint_inference_status = "NOT_DEFINED_BY_CONTRACT"
  )
  structure(c(result_core, list(result_identity = .lt_hash(result_core))),
            class = "lt_log_relative_result")
}

#' Validate a boundary association
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_boundary_association <- function(x) {
  if (!inherits(x, "lt_boundary_association") ||
      !identical(names(x$components), .lt_boundary_component_names()) ||
      !identical(x$joint_inference_status, "NOT_DEFINED_BY_CONTRACT") ||
      "joint_p" %in% names(x)) {
    .lt_abort("Boundary association shape violates the no-joint-p contract.")
  }
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  .lt_assert_unique_ids(x$branch_ids, "branch_ids")
  n <- length(x$gene_ids)
  for (name in .lt_boundary_component_names()) {
    z <- x$components[[name]]
    if (!inherits(z, "lt_boundary_component_association") ||
        !identical(z$component, name) || !identical(z$gene_ids, x$gene_ids) ||
        length(z$statistic) != n || length(z$status) != n ||
        length(z$reason) != n ||
        any(!z$status %in% c("ESTIMABLE", "NONESTIMABLE")) ||
        any(!is.na(z$statistic[z$status != "ESTIMABLE"])) ||
        !identical(z$component_identity,
                   .lt_hash(z[setdiff(names(z), "component_identity")])) ) {
      .lt_abort("Boundary component association is malformed or stale.")
    }
  }
  .lt_assert_sha256(x$association_identity, "association_identity")
  if (!identical(x$association_identity,
                 .lt_hash(.lt_boundary_association_payload(x)))) {
    .lt_abort("Boundary association identity is stale.")
  }
  invisible(x)
}

#' Extract one boundary-aware result component
#' @param x A boundary association or calibration bundle.
#' @param name `ZERO_MASS` or `POSITIVE_LOG`.
#' @param ... Additional method arguments.
#' @return The requested component result.
#' @export
component <- function(x, name, ...) UseMethod("component")

#' @export
component.lt_boundary_association <- function(x, name, ...) {
  validate_lt_boundary_association(x)
  if (!name %in% names(x$components)) {
    .lt_abort("Unknown boundary component; use ZERO_MASS or POSITIVE_LOG.")
  }
  x$components[[name]]
}

.lt_boundary_calibration_payload <- function(x) {
  x[c(
    "calibration_id", "component", "null_hypothesis", "gene_ids",
    "observed_statistic", "observed_status", "null_ensemble_identity",
    "null_replicate_ids", "null_replicate_count",
    "valid_null_replicate_count", "null_nonestimable_count",
    "exceedance_count", "p_empirical", "p_resolution", "monte_carlo_se",
    "calibration_status", "reason", "null_replicate_nonestimable_gene_count",
    "alternative", "domain_identity", "cell_domain_identity", "provenance"
  )]
}

.lt_calibrate_boundary_component <- function(
    object, branch_state, association, ensemble, component_name,
    null_name, alternative, calibration_id) {
  validate_lt_branch_state_ensemble(ensemble)
  if (!identical(ensemble$type, "binary") ||
      !identical(ensemble$trait_id, branch_state$trait_id)) {
    .lt_abort("Boundary null ensemble type/trait differs from observed state.")
  }
  obs_coding <- vapply(
    branch_state$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  null_coding <- vapply(
    ensemble$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  if (!identical(obs_coding, null_coding)) {
    .lt_abort("Observed and null binary coding differ.")
  }
  idx <- match(association$branch_ids, ensemble$branch_ids)
  if (anyNA(idx) || length(setdiff(ensemble$branch_ids,
                                   association$branch_ids))) {
    .lt_abort("Null ensemble branch scientific keys differ from the association.")
  }
  component <- association$components[[component_name]]
  observed <- component$statistic
  n_gene <- length(association$gene_ids)
  r_count <- length(ensemble$replicate_ids)
  active <- component$status == "ESTIMABLE" & is.finite(observed)
  valid_count <- integer(n_gene)
  invalid_count <- integer(n_gene)
  exceedance <- integer(n_gene)
  replicate_bad <- integer(r_count)
  null_values <- ensemble$values[idx, , drop = FALSE]
  null_available <- ensemble$availability[idx, , drop = FALSE] == "available"
  for (j in seq_len(r_count)) {
    values <- null_values[, j]
    values[!null_available[, j]] <- NA
    scored <- .lt_boundary_score(
      object, values, association$branch_allowed, branch_state,
      .validated_token = .lt_boundary_engine_token
    )[[component_name]]
    ok <- active & scored$status == "ESTIMABLE" & is.finite(scored$statistic)
    bad <- active & !ok
    valid_count[ok] <- valid_count[ok] + 1L
    invalid_count[bad] <- invalid_count[bad] + 1L
    replicate_bad[[j]] <- sum(bad)
    hits <- switch(
      alternative,
      greater = scored$statistic >= observed,
      less = scored$statistic <= observed,
      two_sided = abs(scored$statistic) >= abs(observed)
    )
    exceedance[ok & hits] <- exceedance[ok & hits] + 1L
  }
  status <- rep("OBSERVED_DIAGNOSTIC_NONESTIMABLE", n_gene)
  reason <- component$reason
  complete <- active & invalid_count == 0L & valid_count == r_count
  incomplete <- active & !complete
  status[complete] <- "CALIBRATED"
  reason[complete] <- "empirical_null_calibration_complete"
  status[incomplete] <- "NULL_CALIBRATION_INCOMPLETE"
  reason[incomplete] <- paste0(
    "NULL_REPLICATE_NONESTIMABLE_count_", invalid_count[incomplete]
  )
  p <- resolution <- mcse <- rep(NA_real_, n_gene)
  exceedance_out <- rep(NA_integer_, n_gene)
  p[complete] <- (1 + exceedance[complete]) / (r_count + 1)
  resolution[complete] <- 1 / (r_count + 1)
  mcse[complete] <- sqrt(
    p[complete] * (1 - p[complete]) / (r_count + 1)
  )
  exceedance_out[complete] <- as.integer(exceedance[complete])
  replay_domain_identity <-
    object$migration_provenance$source_association_domain_identity %||%
    object$domain_identity
  if (length(replay_domain_identity) != 1L ||
      is.na(replay_domain_identity)) replay_domain_identity <- object$domain_identity
  core <- list(
    calibration_id = calibration_id,
    component = component_name,
    null_hypothesis = null_name,
    gene_ids = association$gene_ids,
    observed_statistic = observed,
    observed_status = component$status,
    null_ensemble_identity = ensemble$null_ensemble_identity,
    null_replicate_ids = ensemble$replicate_ids,
    null_replicate_count = rep.int(as.integer(r_count), n_gene),
    valid_null_replicate_count = as.integer(valid_count),
    null_nonestimable_count = as.integer(invalid_count),
    exceedance_count = exceedance_out,
    p_empirical = p,
    p_resolution = resolution,
    monte_carlo_se = mcse,
    calibration_status = status,
    reason = reason,
    null_replicate_nonestimable_gene_count = as.double(replicate_bad),
    alternative = alternative,
    domain_identity = replay_domain_identity,
    cell_domain_identity = object$domain_identity,
    provenance = list(
      operator = "LaTerra_boundary_component_external_null_calibration",
      one_null_world_applied_to_all_genes = TRUE,
      inclusive_ties = TRUE,
      plus_one_empirical_p = TRUE,
      fail_closed_null_nonestimability = TRUE,
      null_generation = "not_performed",
      sentinel_used = FALSE,
      pooled_components = FALSE,
      pooled_null_hypotheses = FALSE
    )
  )
  out <- structure(c(core, list(
    calibration_identity = .lt_hash(core), metadata = list()
  )), class = "lt_boundary_component_calibration")
  validate_lt_boundary_component_calibration(out)
  out
}

.lt_finalize_boundary_component <- function(
    component, component_name, null_name, ensemble, object,
    valid_count, invalid_count, exceedance, replicate_bad,
    alternative, calibration_id, exact_boundary_recalculation_count = 0L,
    accelerator = "blocked_matrix_products", block_size = 128L,
    tie_tolerance_multiplier = 1e-10) {
  observed <- component$statistic
  n_gene <- length(component$gene_ids)
  r_count <- length(ensemble$replicate_ids)
  active <- component$status == "ESTIMABLE" & is.finite(observed)
  complete <- active & invalid_count == 0L & valid_count == r_count
  incomplete <- active & !complete
  status <- rep("OBSERVED_DIAGNOSTIC_NONESTIMABLE", n_gene)
  reason <- component$reason
  invalid_count[!active] <- r_count
  status[complete] <- "CALIBRATED"
  reason[complete] <- "empirical_null_calibration_complete"
  status[incomplete] <- "NULL_CALIBRATION_INCOMPLETE"
  reason[incomplete] <- paste0(
    "NULL_REPLICATE_NONESTIMABLE_count_", invalid_count[incomplete]
  )
  p <- resolution <- mcse <- rep(NA_real_, n_gene)
  exceedance_out <- rep(NA_integer_, n_gene)
  p[complete] <- (1 + exceedance[complete]) / (r_count + 1)
  resolution[complete] <- 1 / (r_count + 1)
  mcse[complete] <- sqrt(
    p[complete] * (1 - p[complete]) / (r_count + 1)
  )
  exceedance_out[complete] <- as.integer(exceedance[complete])
  replay_domain_identity <-
    object$migration_provenance$source_association_domain_identity %||%
    object$domain_identity
  if (length(replay_domain_identity) != 1L ||
      is.na(replay_domain_identity)) replay_domain_identity <- object$domain_identity
  core <- list(
    calibration_id = calibration_id,
    component = component_name,
    null_hypothesis = null_name,
    gene_ids = component$gene_ids,
    observed_statistic = observed,
    observed_status = component$status,
    null_ensemble_identity = ensemble$null_ensemble_identity,
    null_replicate_ids = ensemble$replicate_ids,
    null_replicate_count = rep.int(as.integer(r_count), n_gene),
    valid_null_replicate_count = as.integer(valid_count),
    null_nonestimable_count = as.integer(invalid_count),
    exceedance_count = exceedance_out,
    p_empirical = p,
    p_resolution = resolution,
    monte_carlo_se = mcse,
    calibration_status = status,
    reason = reason,
    null_replicate_nonestimable_gene_count = as.double(replicate_bad),
    alternative = alternative,
    domain_identity = replay_domain_identity,
    cell_domain_identity = object$domain_identity,
    provenance = list(
      operator = "LaTerra_boundary_component_external_null_calibration",
      one_null_world_applied_to_all_genes = TRUE,
      shared_world_score_pass_for_both_components = TRUE,
      inclusive_ties = TRUE,
      plus_one_empirical_p = TRUE,
      fail_closed_null_nonestimability = TRUE,
      null_generation = "not_performed",
      sentinel_used = FALSE,
      pooled_components = FALSE,
      pooled_null_hypotheses = FALSE,
      accelerator = accelerator,
      accelerator_block_size = as.integer(block_size),
      exact_boundary_recalculation_count =
        as.integer(exact_boundary_recalculation_count),
      boundary_candidate_tolerance_multiplier =
        tie_tolerance_multiplier
    )
  )
  out <- structure(c(core, list(
    calibration_identity = .lt_hash(core), metadata = list()
  )), class = "lt_boundary_component_calibration")
  validate_lt_boundary_component_calibration(out)
  out
}

.lt_calibrate_boundary_null <- function(
    object, branch_state, association, ensemble, null_name,
    alternative, calibration_id) {
  validate_lt_branch_state_ensemble(ensemble)
  if (!identical(ensemble$type, "binary") ||
      !identical(ensemble$trait_id, branch_state$trait_id)) {
    .lt_abort("Boundary null ensemble type/trait differs from observed state.")
  }
  obs_coding <- vapply(
    branch_state$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  null_coding <- vapply(
    ensemble$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  if (!identical(obs_coding, null_coding)) {
    .lt_abort("Observed and null binary coding differ.")
  }
  idx <- match(association$branch_ids, ensemble$branch_ids)
  if (anyNA(idx) || length(setdiff(ensemble$branch_ids,
                                   association$branch_ids))) {
    .lt_abort("Null ensemble branch scientific keys differ from the association.")
  }
  n_gene <- length(association$gene_ids)
  r_count <- length(ensemble$replicate_ids)
  trackers <- lapply(.lt_boundary_component_names(), function(name) {
    component <- association$components[[name]]
    list(
      active = component$status == "ESTIMABLE" &
        is.finite(component$statistic),
      valid = integer(n_gene), invalid = integer(n_gene),
      exceedance = integer(n_gene), replicate_bad = integer(r_count),
      exact_boundary_recalculation_count = 0L
    )
  })
  names(trackers) <- .lt_boundary_component_names()
  null_values <- ensemble$values[idx, , drop = FALSE]
  null_available <- ensemble$availability[idx, , drop = FALSE] == "available"
  state <- as.integer(object$ratio_state)
  dim(state) <- dim(object$ratio_state)
  dimnames(state) <- dimnames(object$ratio_state)
  E_logical <- state %in% c(1L, 2L)
  Z_logical <- state == 2L
  A_logical <- state == 1L
  dim(E_logical) <- dim(state); dimnames(E_logical) <- dimnames(state)
  dim(Z_logical) <- dim(state); dimnames(Z_logical) <- dimnames(state)
  dim(A_logical) <- dim(state); dimnames(A_logical) <- dimnames(state)
  E_logical[, !association$branch_allowed] <- FALSE
  Z_logical[, !association$branch_allowed] <- FALSE
  A_logical[, !association$branch_allowed] <- FALSE
  P <- object$positive_log_values
  P[!A_logical] <- 0
  E <- E_logical; Z <- Z_logical; A <- A_logical
  storage.mode(E) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(A) <- "double"
  storage.mode(P) <- "double"
  total_E <- rowSums(E)
  total_Z <- rowSums(Z)
  total_A <- rowSums(A)
  total_P <- rowSums(P)
  reference_key <- .lt_level_key(ensemble$coding$reference_level)
  focal_key <- .lt_level_key(ensemble$coding$focal_level)
  block_size <- 128L
  tie_tolerance_multiplier <- 1e-10
  constant_availability <- all(null_available)
  exact_positive_score <- function(i, j, focal_mask, reference_mask) {
    focal_cols <- which(focal_mask[, j])
    reference_cols <- which(reference_mask[, j])
    af <- A_logical[i, focal_cols]
    ar <- A_logical[i, reference_cols]
    if (!any(af) || !any(ar)) return(NA_real_)
    yf <- P[i, focal_cols, drop = FALSE]
    yr <- P[i, reference_cols, drop = FALSE]
    yf[!af] <- 0
    yr[!ar] <- 0
    rowSums(yf) / sum(af) - rowSums(yr) / sum(ar)
  }
  for (lo in seq.int(1L, r_count, by = block_size)) {
    hi <- min(r_count, lo + block_size - 1L)
    values <- null_values[, lo:hi, drop = FALSE]
    available <- null_available[, lo:hi, drop = FALSE]
    key_matrix <- matrix(
      vapply(as.vector(values), .lt_level_key, character(1L)),
      nrow = nrow(values), ncol = ncol(values)
    )
    focal_mask <- key_matrix == focal_key & available &
      association$branch_allowed
    F <- focal_mask
    storage.mode(F) <- "double"
    nF <- E %*% F
    zF <- Z %*% F
    aF <- A %*% F
    pF <- P %*% F
    if (constant_availability) {
      reference_mask <- !focal_mask & association$branch_allowed
      nR <- total_E - nF
      zR <- total_Z - zF
      aR <- total_A - aF
      pR <- total_P - pF
    } else {
      reference_mask <- key_matrix == reference_key & available &
        association$branch_allowed
      R <- reference_mask
      storage.mode(R) <- "double"
      nR <- E %*% R
      zR <- Z %*% R
      aR <- A %*% R
      pR <- P %*% R
    }
    score_list <- list(
      ZERO_MASS = zF / nF - zR / nR,
      POSITIVE_LOG = pF / aF - pR / aR
    )
    valid_list <- list(
      ZERO_MASS = nF > 0 & nR > 0 & is.finite(score_list$ZERO_MASS),
      POSITIVE_LOG = aF > 0 & aR > 0 &
        is.finite(score_list$POSITIVE_LOG)
    )
    for (name in .lt_boundary_component_names()) {
      tracker <- trackers[[name]]
      component <- association$components[[name]]
      score <- score_list[[name]]
      valid <- valid_list[[name]]
      if (identical(name, "POSITIVE_LOG")) {
        near <- valid & tracker$active &
          abs(abs(score) - abs(component$statistic)) <=
          tie_tolerance_multiplier * pmax(
            1, abs(score), abs(component$statistic)
          )
        if (any(near)) {
          ij <- which(near, arr.ind = TRUE)
          for (qq in seq_len(nrow(ij))) {
            score[ij[qq, 1L], ij[qq, 2L]] <- exact_positive_score(
              ij[qq, 1L], ij[qq, 2L], focal_mask, reference_mask
            )
          }
          valid <- aF > 0 & aR > 0 & is.finite(score)
          tracker$exact_boundary_recalculation_count <-
            tracker$exact_boundary_recalculation_count + nrow(ij)
        }
      }
      ok <- valid & tracker$active
      bad <- !valid & tracker$active
      hits <- switch(
        alternative,
        greater = ok & score >= component$statistic,
        less = ok & score <= component$statistic,
        two_sided = ok & abs(score) >= abs(component$statistic)
      )
      hits[is.na(hits)] <- FALSE
      tracker$valid <- tracker$valid + rowSums(ok)
      tracker$invalid <- tracker$invalid + rowSums(bad)
      tracker$exceedance <- tracker$exceedance + rowSums(hits)
      tracker$replicate_bad[lo:hi] <- colSums(bad)
      trackers[[name]] <- tracker
    }
  }
  out <- lapply(.lt_boundary_component_names(), function(name) {
    tracker <- trackers[[name]]
    .lt_finalize_boundary_component(
      association$components[[name]], name, null_name, ensemble, object,
      tracker$valid, tracker$invalid, tracker$exceedance,
      tracker$replicate_bad, alternative,
      paste(calibration_id, name, null_name, sep = "__"),
      exact_boundary_recalculation_count =
        tracker$exact_boundary_recalculation_count,
      accelerator = "blocked_matrix_products", block_size = block_size,
      tie_tolerance_multiplier = tie_tolerance_multiplier
    )
  })
  names(out) <- paste(.lt_boundary_component_names(), null_name, sep = "::")
  out
}

#' Calibrate both boundary components against frozen external null ensembles
#' @param object An `lt_log_relative`.
#' @param branch_state Observed binary `lt_branch_state`.
#' @param association Matching `lt_boundary_association`.
#' @param null_ensembles Named list of external `lt_branch_state_ensemble`s.
#' @param alternative Empirical alternative; production authority uses
#'   `two_sided`.
#' @param fail_closed Must be `TRUE`.
#' @param calibration_id Stable bundle identifier.
#' @return An `lt_boundary_calibration` with one member per component and null.
#' @export
lt_calibrate_log_relative <- function(
    object, branch_state, association, null_ensembles,
    alternative = c("two_sided", "greater", "less"), fail_closed = TRUE,
    calibration_id = paste0(association$association_id, "_calibration")) {
  validate_lt_log_relative(object)
  validate_lt_boundary_association(association)
  validate_lt_branch_state(branch_state)
  alternative <- match.arg(alternative)
  if (!isTRUE(fail_closed)) {
    .lt_abort("Boundary calibration supports only fail_closed null handling.")
  }
  if (!identical(object$semantic_identity,
                 association$source_semantic_identity) ||
      !identical(object$audit_identity, association$source_audit_identity) ||
      !identical(branch_state$branch_state_identity,
                 association$branch_state_identity)) {
    .lt_abort_class(
      "Boundary calibration dependencies are stale or incompatible.",
      "lt_error_stale_log_relative_dependency"
    )
  }
  if (!is.list(null_ensembles) || !length(null_ensembles) ||
      is.null(names(null_ensembles)) || any(!nzchar(names(null_ensembles))) ||
      anyDuplicated(names(null_ensembles))) {
    .lt_abort("`null_ensembles` must be a uniquely named non-empty list.")
  }
  calibrations <- do.call(c, unname(lapply(names(null_ensembles), function(null_name) {
    .lt_calibrate_boundary_null(
      object, branch_state, association, null_ensembles[[null_name]],
      null_name, alternative, calibration_id
    )
  })))
  ledger_certificate <- lapply(names(null_ensembles), function(name) {
    ids <- null_ensembles[[name]]$replicate_ids
    list(null_hypothesis = name, replicate_ids = ids,
         replicate_ledger_identity = .lt_hash(ids),
         shared_by_components = TRUE)
  })
  names(ledger_certificate) <- names(null_ensembles)
  core <- list(
    calibration_id = calibration_id,
    source_association_identity = association$association_identity,
    source_audit_identity = object$audit_identity,
    calibrations = calibrations,
    shared_ledger_certificate = ledger_certificate,
    joint_inference_status = "NOT_DEFINED_BY_CONTRACT"
  )
  out <- structure(c(core, list(
    calibration_identity = .lt_hash(core), metadata = list()
  )), class = "lt_boundary_calibration")
  validate_lt_boundary_calibration(out)
  out
}

#' Validate one boundary-component calibration
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_boundary_component_calibration <- function(x) {
  if (!inherits(x, "lt_boundary_component_calibration") ||
      !x$component %in% .lt_boundary_component_names() ||
      !x$alternative %in% c("two_sided", "greater", "less")) {
    .lt_abort("Boundary component calibration contract is malformed.")
  }
  n <- length(x$gene_ids)
  fields <- c(
    "observed_statistic", "observed_status", "null_replicate_count",
    "valid_null_replicate_count", "null_nonestimable_count",
    "exceedance_count", "p_empirical", "p_resolution", "monte_carlo_se",
    "calibration_status", "reason"
  )
  if (any(vapply(x[fields], length, integer(1L)) != n)) {
    .lt_abort("Boundary calibration gene ledgers are misaligned.")
  }
  calibrated <- x$calibration_status == "CALIBRATED"
  if (any(x$valid_null_replicate_count[calibrated] !=
          x$null_replicate_count[calibrated]) ||
      any(x$null_nonestimable_count[calibrated] != 0L) ||
      any(!is.finite(x$p_empirical[calibrated])) ||
      any(!is.na(x$p_empirical[!calibrated])) ||
      any(!is.na(x$exceedance_count[!calibrated]))) {
    .lt_abort("Boundary calibration fail-closed semantics are inconsistent.")
  }
  expected <- (1 + x$exceedance_count[calibrated]) /
    (x$null_replicate_count[calibrated] + 1)
  if (!identical(x$p_empirical[calibrated], expected)) {
    .lt_abort("Boundary empirical p does not use inclusive plus-one semantics.")
  }
  .lt_assert_sha256(x$calibration_identity, "calibration_identity")
  if (!identical(x$calibration_identity,
                 .lt_hash(.lt_boundary_calibration_payload(x)))) {
    .lt_abort("Boundary component calibration identity is stale.")
  }
  invisible(x)
}

#' Validate a boundary calibration bundle
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_boundary_calibration <- function(x) {
  if (!inherits(x, "lt_boundary_calibration") ||
      !is.list(x$calibrations) || !length(x$calibrations) ||
      !identical(x$joint_inference_status, "NOT_DEFINED_BY_CONTRACT") ||
      "joint_p" %in% names(x)) {
    .lt_abort("Boundary calibration violates the separate-family contract.")
  }
  invisible(lapply(x$calibrations,
                   validate_lt_boundary_component_calibration))
  for (name in names(x$shared_ledger_certificate)) {
    cert <- x$shared_ledger_certificate[[name]]
    members <- x$calibrations[vapply(
      x$calibrations, function(z) identical(z$null_hypothesis, name),
      logical(1L)
    )]
    if (length(members) != 2L ||
        !all(vapply(members, function(z)
          identical(z$null_replicate_ids, cert$replicate_ids),
          logical(1L)))) {
      .lt_abort("Components did not consume the identical null replicate ledger.")
    }
  }
  .lt_assert_sha256(x$calibration_identity, "calibration_identity")
  core <- x[setdiff(names(x), c("calibration_identity", "metadata"))]
  if (!identical(x$calibration_identity, .lt_hash(core))) {
    .lt_abort("Boundary calibration bundle identity is stale.")
  }
  invisible(x)
}

#' @export
component.lt_boundary_calibration <- function(x, name, null_hypothesis = NULL,
                                              ...) {
  validate_lt_boundary_calibration(x)
  if (is.null(null_hypothesis)) {
    out <- x$calibrations[vapply(
      x$calibrations, function(z) identical(z$component, name), logical(1L)
    )]
    if (!length(out)) .lt_abort("Unknown boundary component.")
    return(out)
  }
  key <- paste(name, null_hypothesis, sep = "::")
  if (is.null(x$calibrations[[key]])) .lt_abort("Unknown component/null pair.")
  x$calibrations[[key]]
}

#' Apply separate BH/BY families to boundary calibrations
#' @param calibration An `lt_boundary_calibration`.
#' @param methods Adjustment methods; only explicit supported methods are used.
#' @return An `lt_boundary_adjustments` bundle.
#' @export
lt_adjust_log_relative <- function(calibration, methods = c("BH", "BY")) {
  validate_lt_boundary_calibration(calibration)
  methods <- unique(methods)
  if (!length(methods) || any(!methods %in% c("BH", "BY"))) {
    .lt_abort("Boundary production adjustments support explicit BH and/or BY.")
  }
  adjustments <- list()
  for (key in names(calibration$calibrations)) {
    z <- calibration$calibrations[[key]]
    selected <- z$calibration_status == "CALIBRATED"
    for (method in methods) {
      adjusted <- rep(NA_real_, length(z$gene_ids))
      adjusted[selected] <- stats::p.adjust(
        z$p_empirical[selected], method = method
      )
      family <- list(
        component = z$component,
        null_hypothesis = z$null_hypothesis,
        alternative = z$alternative,
        method = method,
        gene_ids = z$gene_ids[selected]
      )
      core <- list(
        component = z$component,
        null_hypothesis = z$null_hypothesis,
        method = method,
        calibration_identity = z$calibration_identity,
        family_identity = .lt_hash(family),
        gene_ids = z$gene_ids,
        n_tests = as.integer(sum(selected)),
        adjusted_p = adjusted,
        provenance = list(
          separate_component_family = TRUE,
          separate_null_hypothesis_family = TRUE,
          joint_p = "not_defined"
        )
      )
      name <- paste(key, method, sep = "::")
      adjustments[[name]] <- structure(c(core, list(
        adjustment_identity = .lt_hash(core)
      )), class = "lt_boundary_adjustment")
    }
  }
  core <- list(
    calibration_identity = calibration$calibration_identity,
    adjustments = adjustments,
    joint_inference_status = "NOT_DEFINED_BY_CONTRACT"
  )
  structure(c(core, list(adjustment_bundle_identity = .lt_hash(core))),
            class = "lt_boundary_adjustments")
}

.lt_legacy_encoded_association_payload <- function(x) {
  x[setdiff(names(x), "association_identity")]
}

.lt_validate_legacy_encoded_association <- function(x, encoded = NULL) {
  if (!inherits(x, "lt_legacy_encoded_association") ||
      !identical(x$mode, "legacy_global_k10") ||
      !identical(x$default_continuous_inference_authorized, FALSE) ||
      !isTRUE(x$acknowledgement_recorded) ||
      !is.character(x$gene_ids) || anyNA(x$gene_ids) ||
      length(x$beta) != length(x$gene_ids) ||
      length(x$status) != length(x$gene_ids) ||
      length(x$reason) != length(x$gene_ids)) {
    .lt_abort("Legacy encoded association is malformed.")
  }
  .lt_assert_sha256(x$source_encoded_identity,
                    "source_encoded_identity")
  # Historical 9016 evidence can be read structurally. Its former
  # acknowledgement fields are inert metadata and confer no authority.
  if (!is.null(encoded)) {
    validate_lt_logGBI_encoded(encoded)
    if (!identical(x$source_encoded_identity, encoded$encoded_identity) ||
        !identical(x$gene_ids, encoded$ordered_gene_ledger)) {
      .lt_abort("Legacy encoded association is stale or incompatible.")
    }
  }
  .lt_assert_sha256(x$association_identity, "association_identity")
  if (!identical(
    x$association_identity,
    .lt_hash(.lt_legacy_encoded_association_payload(x))
  )) {
    .lt_abort("Legacy encoded association identity is stale.")
  }
  invisible(x)
}

.lt_associate_legacy_encoded <- function(encoded, branch_state,
                                         association_id,
                                         authorization = NULL) {
  .lt_refuse_legacy_encoded_inference(
    encoded, ".lt_associate_legacy_encoded"
  )
  validate_lt_logGBI_encoded(encoded)
  validate_lt_branch_state(branch_state)
  if (!identical(branch_state$type, "binary")) {
    .lt_abort("Legacy encoded replay currently supports binary branch state only.")
  }
  idx <- match(encoded$ordered_branch_ledger, branch_state$branch_ids)
  available <- !is.na(idx)
  available[available] <-
    branch_state$availability[idx[available]] == "available"
  state <- rep(NA, length(idx))
  state[!is.na(idx)] <- branch_state$values[idx[!is.na(idx)]]
  keys <- rep(NA_character_, length(state))
  keys[!is.na(state)] <- vapply(
    state[!is.na(state)], .lt_level_key, character(1L)
  )
  ref <- which(keys == .lt_level_key(branch_state$coding$reference_level) &
                 available)
  focal <- which(keys == .lt_level_key(branch_state$coding$focal_level) &
                   available)
  summarize <- function(cols) {
    y <- encoded$values[, cols, drop = FALSE]
    ok <- !is.na(y)
    y[!ok] <- 0
    list(n = as.integer(rowSums(ok)), sum = rowSums(y))
  }
  r <- summarize(ref); f <- summarize(focal)
  beta <- f$sum / f$n - r$sum / r$n
  groups_present <- r$n > 0L & f$n > 0L
  valid <- groups_present & is.finite(beta)
  beta[!valid] <- NA_real_
  beta <- unname(beta)
  status <- rep("NUMERICALLY_INDETERMINATE", length(beta))
  reason <- rep("nonfinite_binary_point_estimate", length(beta))
  no_eligible <- r$n == 0L & f$n == 0L
  reference_absent <- !no_eligible & r$n == 0L
  focal_absent <- !no_eligible & !reference_absent & f$n == 0L
  status[no_eligible] <- "NO_ELIGIBLE_BRANCHES"
  reason[no_eligible] <- "no_cell_passes_layer5_eligibility"
  status[reference_absent] <- "BINARY_REFERENCE_ABSENT"
  reason[reference_absent] <- "no_eligible_reference_branch"
  status[focal_absent] <- "BINARY_FOCAL_ABSENT"
  reason[focal_absent] <- "no_eligible_focal_branch"
  status[valid] <- "ESTIMABLE"
  reason[valid] <- "point_estimate_available"
  core <- list(
    association_id = association_id,
    mode = encoded$mode,
    gene_ids = encoded$ordered_gene_ledger,
    beta = beta,
    n_reference = r$n,
    n_focal = f$n,
    status = status,
    reason = reason,
    source_encoded_identity = encoded$encoded_identity,
    default_continuous_inference_authorized = FALSE,
    acknowledgement_recorded = FALSE,
    provenance = list(
      operator = "LaTerra_explicit_legacy_global_k10_binary_replay",
      sentinel_primary = FALSE,
      boundary_state_inferred_from_sentinel = FALSE,
      purpose = "HISTORICAL_REPLAY_VALIDATION",
      source_m_ref_provenance_identity =
        encoded$m_ref_provenance_identity,
      source_object_identity = encoded$source_object_identity,
      source_domain_identity = encoded$source_domain_identity,
      m_ref = encoded$m_ref,
      L_min_plus = encoded$L_min_plus,
      sentinel = encoded$sentinel,
      m_ref_recomputed_after_subset = FALSE
    )
  )
  out <- structure(c(core, list(association_identity = .lt_hash(core))),
                   class = "lt_legacy_encoded_association")
  .lt_validate_legacy_encoded_association(out, encoded)
  out
}

.lt_calibrate_legacy_encoded_null <- function(
    encoded, branch_state, association, ensemble, null_name, alternative,
    calibration_id, authorization = NULL) {
  .lt_refuse_legacy_encoded_inference(
    encoded, ".lt_calibrate_legacy_encoded_null"
  )
  validate_lt_branch_state_ensemble(ensemble)
  if (!identical(ensemble$type, "binary") ||
      !identical(ensemble$trait_id, branch_state$trait_id)) {
    .lt_abort("Legacy encoded null ensemble type/trait is incompatible.")
  }
  observed_coding <- vapply(
    branch_state$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  null_coding <- vapply(
    ensemble$coding[c("reference_level", "focal_level")],
    .lt_level_key, character(1L)
  )
  if (!identical(observed_coding, null_coding)) {
    .lt_abort("Legacy encoded observed/null coding differs.")
  }
  idx <- match(encoded$ordered_branch_ledger, ensemble$branch_ids)
  if (anyNA(idx) || length(setdiff(ensemble$branch_ids,
                                   encoded$ordered_branch_ledger))) {
    .lt_abort("Legacy encoded null branch keys differ from the encoded view.")
  }
  state_idx <- match(encoded$ordered_branch_ledger, branch_state$branch_ids)
  branch_allowed <- !is.na(state_idx)
  branch_allowed[branch_allowed] <-
    branch_state$availability[state_idx[branch_allowed]] == "available"
  E_logical <- !is.na(encoded$values)
  E_logical[, !branch_allowed] <- FALSE
  Y <- encoded$values
  Y[!E_logical] <- 0
  E <- E_logical
  storage.mode(E) <- "double"
  storage.mode(Y) <- "double"
  total_E <- rowSums(E)
  total_Y <- rowSums(Y)
  null_values <- ensemble$values[idx, , drop = FALSE]
  null_available <- ensemble$availability[idx, , drop = FALSE] == "available"
  r_count <- length(ensemble$replicate_ids)
  n_gene <- length(encoded$ordered_gene_ledger)
  active <- association$status == "ESTIMABLE" & is.finite(association$beta)
  valid_count <- integer(n_gene)
  invalid_count <- integer(n_gene)
  exceedance <- integer(n_gene)
  replicate_bad <- integer(r_count)
  block_size <- 128L
  tie_tolerance_multiplier <- 1e-10
  exact_count <- 0L
  constant_availability <- all(null_available)
  ref_key <- .lt_level_key(ensemble$coding$reference_level)
  focal_key <- .lt_level_key(ensemble$coding$focal_level)
  exact_score <- function(i, j, focal_mask, reference_mask) {
    focal_cols <- which(focal_mask[, j])
    reference_cols <- which(reference_mask[, j])
    ef <- E_logical[i, focal_cols]
    er <- E_logical[i, reference_cols]
    if (!any(ef) || !any(er)) return(NA_real_)
    yf <- Y[i, focal_cols, drop = FALSE]
    yr <- Y[i, reference_cols, drop = FALSE]
    yf[!ef] <- 0
    yr[!er] <- 0
    rowSums(yf) / sum(ef) - rowSums(yr) / sum(er)
  }
  for (lo in seq.int(1L, r_count, by = block_size)) {
    hi <- min(r_count, lo + block_size - 1L)
    values <- null_values[, lo:hi, drop = FALSE]
    available <- null_available[, lo:hi, drop = FALSE]
    keys <- matrix(
      vapply(as.vector(values), .lt_level_key, character(1L)),
      nrow = nrow(values), ncol = ncol(values)
    )
    focal_mask <- keys == focal_key & available & branch_allowed
    F <- focal_mask
    storage.mode(F) <- "double"
    nF <- E %*% F
    yF <- Y %*% F
    if (constant_availability) {
      reference_mask <- !focal_mask & branch_allowed
      nR <- total_E - nF
      yR <- total_Y - yF
    } else {
      reference_mask <- keys == ref_key & available & branch_allowed
      R <- reference_mask
      storage.mode(R) <- "double"
      nR <- E %*% R
      yR <- Y %*% R
    }
    score <- yF / nF - yR / nR
    valid <- nF > 0 & nR > 0 & is.finite(score)
    near <- valid & active &
      abs(abs(score) - abs(association$beta)) <=
      tie_tolerance_multiplier * pmax(
        1, abs(score), abs(association$beta)
      )
    if (any(near)) {
      ij <- which(near, arr.ind = TRUE)
      for (qq in seq_len(nrow(ij))) {
        score[ij[qq, 1L], ij[qq, 2L]] <- exact_score(
          ij[qq, 1L], ij[qq, 2L], focal_mask, reference_mask
        )
      }
      valid <- nF > 0 & nR > 0 & is.finite(score)
      exact_count <- exact_count + nrow(ij)
    }
    ok <- valid & active
    bad <- !valid & active
    hits <- switch(
      alternative,
      greater = ok & score >= association$beta,
      less = ok & score <= association$beta,
      two_sided = ok & abs(score) >= abs(association$beta)
    )
    hits[is.na(hits)] <- FALSE
    valid_count <- valid_count + rowSums(ok)
    invalid_count <- invalid_count + rowSums(bad)
    exceedance <- exceedance + rowSums(hits)
    replicate_bad[lo:hi] <- colSums(bad)
  }
  complete <- active & invalid_count == 0L & valid_count == r_count
  incomplete <- active & !complete
  status <- rep("OBSERVED_DIAGNOSTIC_NONESTIMABLE", n_gene)
  reason <- association$reason
  status[incomplete] <- "NULL_CALIBRATION_INCOMPLETE"
  reason[incomplete] <- paste0(
    "NULL_REPLICATE_NONESTIMABLE_count_", invalid_count[incomplete]
  )
  status[complete] <- "CALIBRATED"
  reason[complete] <- "empirical_null_calibration_complete"
  invalid_count[!active] <- r_count
  exceedance_out <- rep(NA_integer_, n_gene)
  exceedance_out[complete] <- as.integer(exceedance[complete])
  p <- resolution <- mcse <- rep(NA_real_, n_gene)
  p[complete] <- (1 + exceedance_out[complete]) / (r_count + 1)
  resolution[complete] <- 1 / (r_count + 1)
  mcse[complete] <- sqrt(p[complete] * (1 - p[complete]) / (r_count + 1))
  domain_identity <- encoded$source_association_domain_identity
  if (is.na(domain_identity)) {
    domain_identity <- .lt_hash(list(
      encoded_identity = encoded$encoded_identity,
      branch_state_identity = branch_state$branch_state_identity,
      eligible = E_logical
    ))
  }
  cell_domain_identity <- encoded$source_cell_domain_identity
  if (is.na(cell_domain_identity)) {
    cell_domain_identity <- .lt_rate_keyed_mask_hash(
      E_logical, encoded$ordered_gene_ledger, encoded$ordered_branch_ledger
    )
  }
  core <- list(
    calibration_id = paste(calibration_id, null_name, sep = "__"),
    mode = "legacy_global_k10",
    gene_ids = encoded$ordered_gene_ledger,
    observed_statistic = association$beta,
    observed_status = unname(ifelse(
      association$status == "ESTIMABLE", "ESTIMABLE", "NONESTIMABLE"
    )),
    null_hypothesis = null_name,
    null_ensemble_identity = ensemble$null_ensemble_identity,
    null_replicate_ids = ensemble$replicate_ids,
    null_replicate_count = rep.int(as.integer(r_count), n_gene),
    valid_null_replicate_count = as.integer(valid_count),
    null_nonestimable_count = as.integer(invalid_count),
    exceedance_count = exceedance_out,
    p_empirical = p,
    p_resolution = resolution,
    monte_carlo_se = mcse,
    calibration_status = status,
    reason = reason,
    null_replicate_nonestimable_gene_count = as.double(replicate_bad),
    alternative = alternative,
    domain_identity = domain_identity,
    cell_domain_identity = cell_domain_identity,
    provenance = list(
      explicit_encoded_mode = TRUE,
      purpose = "HISTORICAL_REPLAY_VALIDATION",
      default_continuous_inference_authorized = FALSE,
      inclusive_ties = TRUE,
      plus_one_empirical_p = TRUE,
      fail_closed_null_nonestimability = TRUE,
      null_generation = "not_performed",
      accelerator = "blocked_matrix_products",
      accelerator_block_size = block_size,
      exact_boundary_recalculation_count = as.integer(exact_count),
      boundary_candidate_tolerance_multiplier = tie_tolerance_multiplier
    )
  )
  structure(c(core, list(calibration_identity = .lt_hash(core))),
            class = "lt_legacy_encoded_component_calibration")
}

#' Calibrate the explicit historical finite-encoded view
#'
#' This retired entry point always fails closed. Source-bound replay identity
#' validation is available through [validate_legacy_logGBI_replay()], which
#' computes no association statistic.
#'
#' @param encoded An `lt_logGBI_encoded` view.
#' @param branch_state Observed binary `lt_branch_state`.
#' @param association Matching `lt_legacy_encoded_association`.
#' @param null_ensembles Named frozen external null ensembles.
#' @param encoded_mode Deprecated compatibility argument; it cannot authorize
#'   inference.
#' @param acknowledge_encoded_zero Deprecated compatibility argument. Its
#'   value has no effect on the mandatory refusal.
#' @param alternative Empirical alternative.
#' @param calibration_id Stable bundle identifier.
#' @return An `lt_legacy_encoded_calibration`.
#' @export
lt_calibrate_logGBI_encoded <- function(
    encoded, branch_state, association, null_ensembles,
    encoded_mode = NULL, acknowledge_encoded_zero = FALSE,
    alternative = c("two_sided", "greater", "less"),
    calibration_id = paste0(association$association_id, "_calibration")) {
  .lt_refuse_legacy_encoded_inference(
    encoded, "lt_calibrate_logGBI_encoded"
  )
  validate_lt_branch_state(branch_state)
  alternative <- match.arg(alternative)
  .lt_validate_legacy_encoded_association(association, encoded)
  if (!is.list(null_ensembles) || !length(null_ensembles) ||
      is.null(names(null_ensembles)) || any(!nzchar(names(null_ensembles))) ||
      anyDuplicated(names(null_ensembles))) {
    .lt_abort("`null_ensembles` must be a uniquely named non-empty list.")
  }
  calibrations <- lapply(names(null_ensembles), function(name) {
    .lt_calibrate_legacy_encoded_null(
      encoded, branch_state, association, null_ensembles[[name]], name,
      alternative, calibration_id
    )
  })
  names(calibrations) <- names(null_ensembles)
  core <- list(
    calibration_id = calibration_id,
    source_encoded_identity = encoded$encoded_identity,
    source_association_identity = association$association_identity,
    calibrations = calibrations,
    explicit_mode = "legacy_global_k10",
    acknowledgement_recorded = FALSE,
    default_continuous_inference_authorized = FALSE
  )
  out <- structure(c(core, list(calibration_identity = .lt_hash(core))),
                   class = "lt_legacy_encoded_calibration")
  .lt_validate_legacy_encoded_calibration(out)
  out
}

#' Adjust explicit historical encoded calibration families
#' @param calibration An `lt_legacy_encoded_calibration`.
#' @param methods Separate methods, normally `BH` and `BY`.
#' @return An `lt_legacy_encoded_adjustments` bundle.
#' @export
lt_adjust_logGBI_encoded <- function(calibration, methods = c("BH", "BY")) {
  .lt_abort_class(
    .lt_legacy_encoded_inference_message(),
    "lt_error_legacy_encoded_inference_not_supported",
    list(route = "lt_adjust_logGBI_encoded",
         representation = "legacy_global_k10")
  )
  .lt_validate_legacy_encoded_calibration(calibration)
  methods <- unique(methods)
  if (!length(methods) || any(!methods %in% c("BH", "BY"))) {
    .lt_abort("Legacy encoded replay supports explicit BH and/or BY.")
  }
  out <- list()
  for (null_name in names(calibration$calibrations)) {
    z <- calibration$calibrations[[null_name]]
    selected <- z$calibration_status == "CALIBRATED"
    for (method in methods) {
      adjusted <- rep(NA_real_, length(z$gene_ids))
      adjusted[selected] <- stats::p.adjust(z$p_empirical[selected], method)
      core <- list(
        mode = "legacy_global_k10",
        null_hypothesis = null_name,
        method = method,
        calibration_identity = z$calibration_identity,
        gene_ids = z$gene_ids,
        n_tests = as.integer(sum(selected)),
        adjusted_p = adjusted,
        provenance = list(
          explicit_encoded_mode = TRUE,
          separate_null_hypothesis_family = TRUE
        )
      )
      key <- paste(null_name, method, sep = "::")
      out[[key]] <- structure(c(core, list(
        adjustment_identity = .lt_hash(core)
      )), class = "lt_legacy_encoded_adjustment")
    }
  }
  core <- list(
    calibration_identity = calibration$calibration_identity,
    adjustments = out,
    explicit_mode = "legacy_global_k10"
  )
  result <- structure(c(core, list(
    adjustment_bundle_identity = .lt_hash(core)
  )), class = "lt_legacy_encoded_adjustments")
  .lt_validate_legacy_encoded_adjustments(result)
  result
}

.lt_validate_boundary_adjustments <- function(x) {
  if (!inherits(x, "lt_boundary_adjustments") ||
      !is.list(x$adjustments) || !length(x$adjustments) ||
      !identical(x$joint_inference_status, "NOT_DEFINED_BY_CONTRACT")) {
    .lt_abort("Boundary adjustment bundle is malformed.")
  }
  core <- x[setdiff(names(x), "adjustment_bundle_identity")]
  if (!identical(x$adjustment_bundle_identity, .lt_hash(core))) {
    .lt_abort("Boundary adjustment bundle identity is stale.")
  }
  invisible(x)
}

.lt_validate_legacy_encoded_calibration <- function(x) {
  if (!inherits(x, "lt_legacy_encoded_calibration") ||
      !identical(x$explicit_mode, "legacy_global_k10") ||
      !isTRUE(x$acknowledgement_recorded) ||
      !identical(x$default_continuous_inference_authorized, FALSE) ||
      !is.list(x$calibrations) || !length(x$calibrations)) {
    .lt_abort("Legacy encoded calibration bundle is malformed.")
  }
  # Historical 9016 acknowledgement fields are inert compatibility metadata.
  for (z in x$calibrations) {
    if (!inherits(z, "lt_legacy_encoded_component_calibration") ||
        !identical(z$mode, "legacy_global_k10") ||
        length(z$gene_ids) != length(z$p_empirical) ||
        any(!is.na(z$p_empirical[z$calibration_status != "CALIBRATED"]))) {
      .lt_abort("Legacy encoded component calibration is malformed.")
    }
    core <- z[setdiff(names(z), "calibration_identity")]
    if (!identical(z$calibration_identity, .lt_hash(core))) {
      .lt_abort("Legacy encoded component calibration identity is stale.")
    }
  }
  core <- x[setdiff(names(x), "calibration_identity")]
  if (!identical(x$calibration_identity, .lt_hash(core))) {
    .lt_abort("Legacy encoded calibration bundle identity is stale.")
  }
  invisible(x)
}

.lt_validate_legacy_encoded_adjustments <- function(x) {
  if (!inherits(x, "lt_legacy_encoded_adjustments") ||
      !identical(x$explicit_mode, "legacy_global_k10") ||
      !is.list(x$adjustments) || !length(x$adjustments)) {
    .lt_abort("Legacy encoded adjustment bundle is malformed.")
  }
  core <- x[setdiff(names(x), "adjustment_bundle_identity")]
  if (!identical(x$adjustment_bundle_identity, .lt_hash(core))) {
    .lt_abort("Legacy encoded adjustment bundle identity is stale.")
  }
  invisible(x)
}

#' @export
print.lt_boundary_association <- function(x, ...) {
  validate_lt_boundary_association(x)
  cat("<lt_boundary_association>", x$association_id, "\n")
  for (name in .lt_boundary_component_names()) {
    z <- x$components[[name]]
    cat(" ", name, ":", sum(z$status == "ESTIMABLE"), "/",
        length(z$status), "estimable genes\n")
  }
  cat("  joint p: not defined\n")
  invisible(x)
}

#' @export
print.lt_boundary_calibration <- function(x, ...) {
  validate_lt_boundary_calibration(x)
  cat("<lt_boundary_calibration>", x$calibration_id, "\n")
  for (z in x$calibrations) {
    cat(" ", z$component, "x", z$null_hypothesis, ":",
        sum(z$calibration_status == "CALIBRATED"), "calibrated genes\n")
  }
  cat("  joint p: not defined\n")
  invisible(x)
}

#' @export
print.lt_log_relative_result <- function(x, ...) {
  cat("<lt_log_relative_result>\n")
  print(x$association)
  cat("  calibrated component/null families:",
      length(x$calibration$calibrations), "\n")
  cat("  joint p: not defined\n")
  invisible(x)
}
