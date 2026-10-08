# Empirical calibration against externally supplied null worlds.

.lt_calibration_statistics <- function() {
  c("beta", "pearson_r", "spearman_rho", "point_biserial_r")
}

.lt_adjustment_methods <- function() {
  c("none", "BH", "BY", "Holm", "Bonferroni")
}

.lt_calibration_spec_payload <- function(x) {
  x[c("statistic", "alternative", "multiple_testing", "domain_mode",
      "nonestimable_null_policy")]
}

#' Declare an empirical association-calibration specification
#'
#' @param statistic Explicit T2A statistic.
#' @param alternative `greater`, `less`, or `two_sided`.
#' @param multiple_testing Requested later adjustment method. Calibration never
#'   applies it silently; call [lt_adjust_calibration()] explicitly.
#' @param domain_mode Native or common association-domain mode.
#' @param nonestimable_null_policy T2B currently supports only `fail_closed`.
#' @param metadata Named non-authoritative metadata list.
#' @return An `lt_calibration_spec`.
#' @export
lt_calibration_spec <- function(
    statistic = "beta", alternative = c("two_sided", "greater", "less"),
    multiple_testing = c("none", "BH", "BY", "Holm", "Bonferroni"),
    domain_mode = c("NATIVE_ASSOCIATION_DOMAIN",
                    "COMMON_ASSOCIATION_DOMAIN"),
    nonestimable_null_policy = "fail_closed", metadata = list()) {
  alternative <- match.arg(alternative)
  multiple_testing <- match.arg(multiple_testing)
  domain_mode <- match.arg(domain_mode)
  core <- list(
    statistic = statistic,
    alternative = alternative,
    multiple_testing = multiple_testing,
    domain_mode = domain_mode,
    nonestimable_null_policy = nonestimable_null_policy
  )
  out <- structure(c(core, list(
    calibration_spec_identity = .lt_hash(core), metadata = metadata
  )), class = "lt_calibration_spec")
  validate_lt_calibration_spec(out)
  out
}

#' Validate a calibration specification
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_calibration_spec <- function(x) {
  if (!inherits(x, "lt_calibration_spec")) {
    .lt_abort("`x` must inherit from `lt_calibration_spec`.")
  }
  if (!is.character(x$statistic) || length(x$statistic) != 1L ||
      !x$statistic %in% .lt_calibration_statistics()) {
    .lt_abort("Unknown T2B calibration statistic.")
  }
  if (!x$alternative %in% c("two_sided", "greater", "less")) {
    .lt_abort("Unknown T2B alternative.")
  }
  if (!x$multiple_testing %in% .lt_adjustment_methods()) {
    .lt_abort("Unknown multiple-testing method.")
  }
  if (!x$domain_mode %in% c("NATIVE_ASSOCIATION_DOMAIN",
                            "COMMON_ASSOCIATION_DOMAIN")) {
    .lt_abort("Unknown calibration domain mode.")
  }
  if (!identical(x$nonestimable_null_policy, "fail_closed")) {
    .lt_abort("T2B supports only fail_closed non-estimable-null handling.")
  }
  .lt_assert_named_list(x$metadata, "metadata")
  .lt_assert_sha256(x$calibration_spec_identity,
                    "calibration_spec_identity")
  if (!identical(x$calibration_spec_identity,
                 .lt_hash(.lt_calibration_spec_payload(x)))) {
    .lt_abort("Calibration-spec semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

.lt_observed_statistic <- function(association, statistic) {
  switch(
    statistic,
    beta = association$estimates$beta,
    pearson_r = association$descriptive_scores$pearson_r,
    spearman_rho = association$descriptive_scores$spearman_rho,
    point_biserial_r = association$descriptive_scores$point_biserial_r
  )
}

.lt_validate_statistic_type <- function(type, statistic) {
  allowed <- if (type == "binary") c("beta", "point_biserial_r") else
    c("beta", "pearson_r", "spearman_rho")
  if (!statistic %in% allowed) {
    .lt_abort(paste0("Statistic `", statistic,
                    "` is not defined for observed branch-state type `",
                    type, "`."))
  }
}

.lt_null_score <- function(rate, domain, state_values, type, coding,
                           statistic) {
  .lt_refuse_legacy_encoded_inference(rate, ".lt_null_score")
  .lt_refuse_historical_logGBI_inference(rate, ".lt_null_score")
  n_gene <- length(domain$gene_ids)
  score <- rep(NA_real_, n_gene)
  status <- rep("NULL_REPLICATE_NONESTIMABLE", n_gene)
  if (type == "binary") {
    state_keys <- vapply(state_values, .lt_level_key, character(1L))
    ref_key <- .lt_level_key(coding$reference_level)
    focal_key <- .lt_level_key(coding$focal_level)
    ref_cols <- which(state_keys == ref_key)
    focal_cols <- which(state_keys == focal_key)
    summarize <- function(columns) {
      ok <- domain$eligible[, columns, drop = FALSE]
      y <- rate$values[, columns, drop = FALSE]
      n <- rowSums(ok)
      y[!ok] <- 0
      list(n = as.integer(n), sum = rowSums(y), ss = rowSums(y * y))
    }
    ref <- summarize(ref_cols); focal <- summarize(focal_cols)
    valid_groups <- ref$n > 0L & focal$n > 0L
    beta <- focal$sum / focal$n - ref$sum / ref$n
    if (statistic == "beta") {
      score[valid_groups] <- beta[valid_groups]
    } else {
      n <- ref$n + focal$n
      total_sum <- ref$sum + focal$sum
      total_ss <- ref$ss + focal$ss - total_sum^2 / pmax(n, 1L)
      total_ss[total_ss < 0 & total_ss > -1e-12 *
        pmax(ref$ss + focal$ss, 1)] <- 0
      r <- beta * sqrt(ref$n * focal$n) / sqrt(n * total_ss)
      score[valid_groups] <- r[valid_groups]
    }
    finite <- valid_groups & is.finite(score)
    status[finite] <- "ESTIMABLE"
    status[valid_groups & !finite] <- "NUMERICALLY_INDETERMINATE"
    return(list(score = score, status = status))
  }
  for (i in seq_len(n_gene)) {
    keep <- domain$eligible[i, ]
    x <- as.double(state_values[keep])
    y <- as.double(rate$values[i, keep])
    if (length(x) < 2L || length(unique(x)) < 2L) next
    value <- switch(
      statistic,
      beta = {
        xm <- mean(x); ym <- mean(y)
        sum((x - xm) * (y - ym)) / sum((x - xm)^2)
      },
      pearson_r = .lt_safe_cor(x, y, "pearson"),
      spearman_rho = .lt_safe_cor(x, y, "spearman")
    )
    if (is.finite(value)) {
      score[[i]] <- value
      status[[i]] <- "ESTIMABLE"
    } else {
      status[[i]] <- "NUMERICALLY_INDETERMINATE"
    }
  }
  list(score = score, status = status)
}

.lt_calibration_payload <- function(x) {
  payload <- x[c(
    "calibration_id", "association_identity", "rate_identity",
    "branch_state_identity", "null_ensemble_identity",
    "statistic", "alternative", "trait_id",
    "rate_representation", "domain_mode", "gene_ids",
    "observed_statistic", "null_replicate_count",
    "valid_null_replicate_count", "exceedance_count", "p_empirical",
    "p_resolution", "monte_carlo_se", "calibration_status", "reason",
    "null_nonestimable_count", "null_replicate_status",
    "joint_null_max_abs", "domain_identity", "provenance"
  )]
  payload$calibration_spec_identity <-
    x$calibration_spec$calibration_spec_identity
  payload
}

#' Calibrate a T2A association against external null worlds
#'
#' The association domain is frozen before null scores are inspected. The
#' engine contains no random-number generation and applies each branch-state
#' replicate to every gene. Any non-estimable null replicate fails the affected
#' gene calibration closed. Historical finite-sentinel `lt_rate` objects whose
#' representation is `logGBI` are forbidden in this ordinary route. Finite
#' encoded views are diagnostic/replay-only and are also refused. Use GBI for
#' La Terra V1 inference.
#'
#' @param rate Frozen `lt_rate` used by the observed association.
#' @param branch_state Observed `lt_branch_state`.
#' @param association Frozen observed `lt_association`.
#' @param domain Frozen `lt_association_domain`.
#' @param null_ensemble External `lt_branch_state_ensemble`.
#' @param calibration_spec An `lt_calibration_spec`.
#' @param calibration_id Stable result identifier.
#' @param metadata Named non-authoritative metadata list.
#' @return An `lt_calibration` retaining every gene.
#' @export
lt_calibrate_association <- function(
    rate, branch_state, association, domain, null_ensemble,
    calibration_spec = lt_calibration_spec(
      domain_mode = domain$domain_mode
    ),
    calibration_id = paste0(association$association_id, "_",
                            calibration_spec$statistic, "_",
                            calibration_spec$alternative, "_calibration"),
    metadata = list()) {
  .lt_refuse_legacy_encoded_inference(rate, "lt_calibrate_association")
  .lt_refuse_historical_logGBI_inference(rate, "lt_calibrate_association")
  validate_lt_rate(rate)
  validate_lt_branch_state(branch_state)
  validate_lt_association(association)
  validate_lt_association_domain(domain)
  validate_lt_branch_state_ensemble(null_ensemble)
  validate_lt_calibration_spec(calibration_spec)
  .lt_assert_scalar_character(calibration_id, "calibration_id")
  .lt_assert_named_list(metadata, "metadata")
  .lt_validate_statistic_type(branch_state$type,
                              calibration_spec$statistic)
  if (!identical(branch_state$type, null_ensemble$type)) {
    .lt_abort("Observed and null branch-state types must match exactly.")
  }
  if (!identical(branch_state$trait_id, null_ensemble$trait_id)) {
    .lt_abort("Observed and null trait identities must match exactly.")
  }
  if (branch_state$type == "binary") {
    observed_coding <- vapply(
      branch_state$coding[c("reference_level", "focal_level")],
      .lt_level_key, character(1L)
    )
    null_coding <- vapply(
      null_ensemble$coding[c("reference_level", "focal_level")],
      .lt_level_key, character(1L)
    )
    if (!identical(observed_coding, null_coding)) {
      .lt_abort("Observed and null binary coding must match exactly.")
    }
  }
  if (length(setdiff(branch_state$branch_ids, null_ensemble$branch_ids)) ||
      length(setdiff(null_ensemble$branch_ids, branch_state$branch_ids))) {
    .lt_abort("Observed and null branch IDs must be the same exact scientific-key set.")
  }
  if (!identical(rate$ordered_gene_ledger, domain$gene_ids) ||
      !identical(rate$ordered_branch_ledger, domain$branch_ids) ||
      !identical(association$gene_ids, domain$gene_ids) ||
      !identical(association$association_domain_identity,
                 domain$association_domain_identity) ||
      !identical(association$branch_state_identity,
                 branch_state$branch_state_identity) ||
      !identical(association$rate_identity,
                 .lt_rate_association_identity(rate)$semantic_object_sha256)) {
    .lt_abort("Calibration inputs are not the exact scientific dependencies of the frozen T2A association/domain.")
  }
  if (!identical(calibration_spec$domain_mode, domain$domain_mode)) {
    .lt_abort("Calibration-spec domain mode differs from the frozen association domain.")
  }
  ensemble_index <- match(domain$branch_ids, null_ensemble$branch_ids)
  if (anyNA(ensemble_index)) {
    .lt_abort("Null ensemble cannot be joined to the domain by exact branch key.")
  }
  state_index <- match(domain$branch_ids, branch_state$branch_ids)
  if (anyNA(state_index)) {
    .lt_abort("Observed branch state cannot be joined to the domain by exact branch key.")
  }
  observed <- as.double(.lt_observed_statistic(
    association, calibration_spec$statistic
  ))
  n_gene <- length(domain$gene_ids)
  r_count <- length(null_ensemble$replicate_ids)
  status <- rep("CALIBRATED", n_gene)
  reason <- rep("empirical_null_calibration_complete", n_gene)
  observed_bad <- association$status != "ESTIMABLE"
  status[observed_bad] <- "OBSERVED_ASSOCIATION_NONESTIMABLE"
  reason[observed_bad] <- association$reason[observed_bad]
  numeric_bad <- !observed_bad & !is.finite(observed)
  status[numeric_bad] <- "NUMERICALLY_INDETERMINATE"
  reason[numeric_bad] <- "selected_observed_statistic_is_nonfinite"

  null_available <- null_ensemble$availability[ensemble_index, , drop = FALSE] ==
    "available"
  support_matches <- vapply(seq_len(n_gene), function(i) {
    required <- domain$eligible[i, ]
    !any(required) || all(null_available[required, , drop = FALSE])
  }, logical(1L))
  provenance <- list(
    operator = "LaTerra_T2B_external_null_empirical_calibration",
    calibration_computation_validated = TRUE,
    null_generator_scientific_validity = "EXTERNAL_NOT_ADJUDICATED",
    one_null_world_applied_to_all_genes = TRUE,
    domain_frozen_before_null_scores = TRUE,
    exact_key_join = TRUE,
    hidden_rng = FALSE,
    empirical_p_ties = "inclusive",
    empirical_p_correction = "plus_one",
    neutral_point = 0,
    streaming = "replicate_wise_no_gene_by_replicate_matrix",
    streaming_passes = 2L,
    joint_null_summary = "maximum_absolute_statistic_across_complete_genes",
    automatic_multiple_testing = FALSE,
    null_generation = "not_performed",
    rate_preprocessing = "none",
    zero_encoded_cells_consumed_as_is = FALSE,
    zero_sentinel_recomputed = FALSE
  )
  domain_bad <- status == "CALIBRATED" & !support_matches
  status[domain_bad] <- "NULL_DOMAIN_MISMATCH"
  reason[domain_bad] <-
    "null_availability_does_not_cover_gene_frozen_association_domain"
  null_values <- null_ensemble$values[ensemble_index, , drop = FALSE]
  # On unavailable branches the score never enters the frozen domain; replace
  # the storage NA only to keep type-specific score helpers total.
  if (branch_state$type == "continuous") {
    null_values[!null_available] <- 0
  } else {
    null_values[!null_available] <- null_ensemble$coding$reference_level
  }
  valid_count <- integer(n_gene)
  invalid_count <- integer(n_gene)
  exceedance <- integer(n_gene)
  replicate_bad <- integer(r_count)
  active <- status == "CALIBRATED"
  for (j in seq_len(r_count)) {
    scored <- .lt_null_score(
      rate, domain, null_values[, j], branch_state$type,
      branch_state$coding, calibration_spec$statistic
    )
    ok <- active & scored$status == "ESTIMABLE" & is.finite(scored$score)
    bad <- active & !ok
    valid_count[ok] <- valid_count[ok] + 1L
    invalid_count[bad] <- invalid_count[bad] + 1L
    replicate_bad[[j]] <- sum(bad)
    hits <- switch(
      calibration_spec$alternative,
      greater = scored$score >= observed,
      less = scored$score <= observed,
      two_sided = abs(scored$score) >= abs(observed)
    )
    exceedance[ok & hits] <- exceedance[ok & hits] + 1L
  }
  incomplete <- active & invalid_count > 0L
  status[incomplete] <- "NULL_CALIBRATION_INCOMPLETE"
  reason[incomplete] <- paste0(
    "NULL_REPLICATE_NONESTIMABLE_count_", invalid_count[incomplete]
  )
  complete <- active & invalid_count == 0L & valid_count == r_count
  p <- rep(NA_real_, n_gene)
  resolution <- rep(NA_real_, n_gene)
  mcse <- rep(NA_real_, n_gene)
  p[complete] <- (1 + exceedance[complete]) / (r_count + 1)
  resolution[complete] <- 1 / (r_count + 1)
  mcse[complete] <- sqrt(p[complete] * (1 - p[complete]) / (r_count + 1))
  exceedance_out <- rep(NA_integer_, n_gene)
  exceedance_out[complete] <- exceedance[complete]
  joint_max <- rep(NA_real_, r_count)
  if (any(complete)) {
    for (j in seq_len(r_count)) {
      scored <- .lt_null_score(
        rate, domain, null_values[, j], branch_state$type,
        branch_state$coding, calibration_spec$statistic
      )
      joint_max[[j]] <- max(abs(scored$score[complete]))
    }
  }
  null_replicate_status <- data.frame(
    replicate_id = null_ensemble$replicate_ids,
    status = ifelse(replicate_bad == 0L, "ESTIMABLE_FOR_ALL_ACTIVE_GENES",
                    "NULL_REPLICATE_NONESTIMABLE"),
    nonestimable_gene_count = as.integer(replicate_bad),
    stringsAsFactors = FALSE
  )
  core <- list(
    calibration_id = calibration_id,
    association_identity = association$association_identity,
    rate_identity = .lt_rate_association_identity(rate)$semantic_object_sha256,
    branch_state_identity = branch_state$branch_state_identity,
    null_ensemble_identity = null_ensemble$null_ensemble_identity,
    calibration_spec = calibration_spec,
    statistic = calibration_spec$statistic,
    alternative = calibration_spec$alternative,
    trait_id = association$trait_id,
    rate_representation = rate$representation,
    domain_mode = domain$domain_mode,
    gene_ids = domain$gene_ids,
    observed_statistic = observed,
    null_replicate_count = rep.int(as.integer(r_count), n_gene),
    valid_null_replicate_count = as.integer(valid_count),
    exceedance_count = exceedance_out,
    p_empirical = p,
    p_resolution = resolution,
    monte_carlo_se = mcse,
    calibration_status = status,
    reason = reason,
    null_nonestimable_count = as.integer(invalid_count),
    null_replicate_status = null_replicate_status,
    joint_null_max_abs = joint_max,
    domain_identity = domain$association_domain_identity,
    provenance = provenance
  )
  proto <- structure(c(core, list(metadata = metadata)),
                     class = "lt_calibration")
  out <- structure(c(core, list(
    calibration_identity = .lt_hash(.lt_calibration_payload(proto)),
    metadata = metadata
  )), class = "lt_calibration")
  validate_lt_calibration(out)
  out
}

.lt_calibration_status_levels <- function() {
  c(
    "CALIBRATED", "OBSERVED_ASSOCIATION_NONESTIMABLE",
    "NULL_DOMAIN_MISMATCH", "NULL_REPLICATE_NONESTIMABLE",
    "NULL_CALIBRATION_INCOMPLETE", "NUMERICALLY_INDETERMINATE"
  )
}

#' Validate an empirical calibration result
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_calibration <- function(x) {
  if (!inherits(x, "lt_calibration")) {
    .lt_abort("`x` must inherit from `lt_calibration`.")
  }
  .lt_assert_scalar_character(x$calibration_id, "calibration_id")
  .lt_assert_sha256(x$association_identity, "association_identity")
  .lt_assert_sha256(x$rate_identity, "rate_identity")
  .lt_assert_sha256(x$branch_state_identity, "branch_state_identity")
  .lt_assert_sha256(x$null_ensemble_identity, "null_ensemble_identity")
  validate_lt_calibration_spec(x$calibration_spec)
  if (!identical(x$statistic, x$calibration_spec$statistic) ||
      !identical(x$alternative, x$calibration_spec$alternative) ||
      !identical(x$domain_mode, x$calibration_spec$domain_mode)) {
    .lt_abort("Calibration result and calibration specification disagree.")
  }
  .lt_assert_scalar_character(x$trait_id, "trait_id")
  .lt_assert_scalar_character(x$rate_representation, "rate_representation")
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  n <- length(x$gene_ids)
  vector_n <- c(
    "observed_statistic", "null_replicate_count",
    "valid_null_replicate_count", "exceedance_count", "p_empirical",
    "p_resolution", "monte_carlo_se", "calibration_status", "reason",
    "null_nonestimable_count"
  )
  if (any(vapply(x[vector_n], length, integer(1L)) != n)) {
    .lt_abort("Calibration gene-level fields must retain the exact gene ledger.")
  }
  if (!is.numeric(x$observed_statistic) ||
      !is.character(x$calibration_status) || !is.character(x$reason) ||
      anyNA(x$calibration_status) ||
      any(!x$calibration_status %in% .lt_calibration_status_levels()) ||
      anyNA(x$reason) || any(!nzchar(x$reason))) {
    .lt_abort("Calibration observed-statistic/status/reason ledger is malformed.")
  }
  calibrated <- x$calibration_status == "CALIBRATED"
  if (any(!is.finite(x$p_empirical[calibrated])) ||
      any(x$p_empirical[calibrated] <= 0 | x$p_empirical[calibrated] > 1) ||
      any(!is.finite(x$p_resolution[calibrated])) ||
      any(!is.finite(x$monte_carlo_se[calibrated])) ||
      any(is.na(x$exceedance_count[calibrated])) ||
      any(x$valid_null_replicate_count[calibrated] !=
            x$null_replicate_count[calibrated]) ||
      any(x$null_nonestimable_count[calibrated] != 0L)) {
    .lt_abort("Calibrated genes have incomplete or invalid empirical-p evidence.")
  }
  counted <- x$calibration_status %in%
    c("CALIBRATED", "NULL_CALIBRATION_INCOMPLETE")
  if (!is.integer(x$null_replicate_count) ||
      !is.integer(x$valid_null_replicate_count) ||
      !is.integer(x$exceedance_count) ||
      !is.integer(x$null_nonestimable_count) ||
      any(x$null_replicate_count < 1L) ||
      any(x$valid_null_replicate_count < 0L) ||
      any(x$null_nonestimable_count < 0L) ||
      any(x$valid_null_replicate_count[counted] +
            x$null_nonestimable_count[counted] !=
            x$null_replicate_count[counted]) ||
      any(x$valid_null_replicate_count[!counted] != 0L) ||
      any(x$null_nonestimable_count[!counted] != 0L)) {
    .lt_abort("Calibration replicate counts are malformed or inconsistent.")
  }
  incomplete <- x$calibration_status == "NULL_CALIBRATION_INCOMPLETE"
  if (any(x$null_nonestimable_count[incomplete] < 1L)) {
    .lt_abort("Incomplete calibrations require at least one non-estimable null replicate.")
  }
  expected_p <- (1 + x$exceedance_count[calibrated]) /
    (x$null_replicate_count[calibrated] + 1)
  expected_resolution <- 1 / (x$null_replicate_count[calibrated] + 1)
  expected_mcse <- sqrt(expected_p * (1 - expected_p) /
                          (x$null_replicate_count[calibrated] + 1))
  if (!identical(x$p_empirical[calibrated], expected_p) ||
      !identical(x$p_resolution[calibrated], expected_resolution) ||
      !identical(x$monte_carlo_se[calibrated], expected_mcse)) {
    .lt_abort("Empirical p-value or Monte-Carlo precision formula is stale.")
  }
  if (any(!is.na(x$p_empirical[!calibrated])) ||
      any(!is.na(x$p_resolution[!calibrated])) ||
      any(!is.na(x$monte_carlo_se[!calibrated])) ||
      any(!is.na(x$exceedance_count[!calibrated]))) {
    .lt_abort("Non-calibrated genes must not release empirical-p evidence.")
  }
  if (!is.data.frame(x$null_replicate_status) || !identical(
    names(x$null_replicate_status),
    c("replicate_id", "status", "nonestimable_gene_count")
  ) || !is.character(x$null_replicate_status$replicate_id) ||
      anyNA(x$null_replicate_status$replicate_id) ||
      any(!nzchar(x$null_replicate_status$replicate_id)) ||
      anyDuplicated(x$null_replicate_status$replicate_id) ||
      !is.character(x$null_replicate_status$status) ||
      anyNA(x$null_replicate_status$status) ||
      any(!x$null_replicate_status$status %in%
            c("ESTIMABLE_FOR_ALL_ACTIVE_GENES",
              "NULL_REPLICATE_NONESTIMABLE")) ||
      !is.integer(x$null_replicate_status$nonestimable_gene_count) ||
      any(x$null_replicate_status$nonestimable_gene_count < 0L)) {
    .lt_abort("Null-replicate status ledger is malformed.")
  }
  if (length(unique(x$null_replicate_count)) != 1L ||
      unique(x$null_replicate_count) != nrow(x$null_replicate_status)) {
    .lt_abort("Calibration replicate ledger and gene-level replicate counts differ.")
  }
  if (!is.numeric(x$joint_null_max_abs) ||
      length(x$joint_null_max_abs) != nrow(x$null_replicate_status)) {
    .lt_abort("Joint-null streaming summary is malformed.")
  }
  if ((any(calibrated) && any(!is.finite(x$joint_null_max_abs))) ||
      (!any(calibrated) && any(!is.na(x$joint_null_max_abs)))) {
    .lt_abort("Joint-null summary availability is inconsistent with calibrated genes.")
  }
  .lt_assert_sha256(x$domain_identity, "domain_identity")
  .lt_assert_named_list(x$provenance, "provenance")
  .lt_assert_named_list(x$metadata, "metadata")
  .lt_assert_sha256(x$calibration_identity, "calibration_identity")
  if (!identical(x$calibration_identity,
                 .lt_hash(.lt_calibration_payload(x)))) {
    .lt_abort("Calibration semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

print.lt_calibration_spec <- function(x, ...) {
  cat("<lt_calibration_spec>", x$statistic, x$alternative, "\n")
  cat("  domain:", x$domain_mode, " adjustment request:",
      x$multiple_testing, "(not automatically applied)\n")
  invisible(x)
}

print.lt_calibration <- function(x, ...) {
  cat("<lt_calibration>", x$calibration_id, "\n")
  cat("  statistic:", x$statistic, " alternative:", x$alternative, "\n")
  cat("  calibrated:", sum(x$calibration_status == "CALIBRATED"), "/",
      length(x$gene_ids), " genes\n")
  cat("  null generator scientific validity: external / not adjudicated\n")
  invisible(x)
}
