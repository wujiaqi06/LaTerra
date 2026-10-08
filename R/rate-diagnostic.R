# Trait-free, threshold-free diagnostics for materialized rate objects.

.lt_rate_diag_probs <- function() {
  c(0, .001, .01, .05, .25, .5, .75, .95, .99, .999, 1)
}

.lt_rate_diag_qnames <- function() {
  c("min", "q001", "q01", "q05", "q25", "median", "q75",
    "q95", "q99", "q999", "max")
}

.lt_rate_neutral_value <- function(representation) {
  if (representation %in% c("GBI", "GBI_C2")) 1 else 0
}

.lt_rate_diag_safe_quantile <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(stats::setNames(rep(NA_real_, length(.lt_rate_diag_probs())),
                           .lt_rate_diag_qnames()))
  }
  stats::setNames(as.double(stats::quantile(
    x, .lt_rate_diag_probs(), names = FALSE, type = 7
  )), .lt_rate_diag_qnames())
}

.lt_rate_diag_concentration <- function(departure) {
  departure <- departure[is.finite(departure) & departure >= 0]
  total <- sum(departure)
  if (!length(departure) || !is.finite(total) || total == 0) {
    return(c(top_single = NA_real_, top_0_1pct = NA_real_,
             top_1pct = NA_real_, top_5pct = NA_real_))
  }
  departure <- sort(departure, decreasing = TRUE, method = "radix")
  share <- function(fraction) {
    sum(departure[seq_len(max(1L, ceiling(length(departure) * fraction)))]) /
      total
  }
  c(top_single = departure[[1L]] / total,
    top_0_1pct = share(.001), top_1pct = share(.01), top_5pct = share(.05))
}

.lt_rate_diag_vector <- function(value, neutral) {
  value <- value[is.finite(value)]
  departure <- abs(value - neutral)
  q <- .lt_rate_diag_safe_quantile(value)
  dq <- .lt_rate_diag_safe_quantile(departure)
  cc <- .lt_rate_diag_concentration(departure)
  c(
    available_count = length(value), zero_count = sum(value == 0),
    positive_count = sum(value > 0), negative_count = sum(value < 0),
    q,
    mean = if (length(value)) mean(value) else NA_real_,
    sd = if (length(value) > 1L) stats::sd(value) else NA_real_,
    iqr = if (length(value)) stats::IQR(value, type = 7) else NA_real_,
    mad = if (length(value)) stats::mad(value) else NA_real_,
    departure_median = unname(dq[["median"]]),
    departure_q90 = if (length(departure)) as.double(stats::quantile(
      departure, .9, names = FALSE, type = 7
    )) else NA_real_,
    departure_q95 = if (length(departure)) as.double(stats::quantile(
      departure, .95, names = FALSE, type = 7
    )) else NA_real_,
    departure_q99 = if (length(departure)) as.double(stats::quantile(
      departure, .99, names = FALSE, type = 7
    )) else NA_real_,
    departure_q999 = unname(dq[["q999"]]),
    departure_max = unname(dq[["max"]]),
    departure_total = sum(departure),
    cc
  )
}

.lt_rate_diag_raw_vector <- function(value) {
  value <- value[is.finite(value)]
  c(
    raw_observed_count = length(value),
    raw_sum = if (length(value)) sum(value) else NA_real_,
    raw_mean = if (length(value)) mean(value) else NA_real_,
    raw_median = if (length(value)) stats::median(value) else NA_real_,
    raw_zero_count = sum(value == 0),
    raw_positive_count = sum(value > 0),
    raw_zero_fraction = if (length(value)) mean(value == 0) else NA_real_
  )
}

.lt_rate_diag_group <- function(values, neutral, keys, margin,
                                raw_values = NULL, baseline_mu = NULL,
                                zero_encoded = NULL) {
  n <- if (margin == 1L) nrow(values) else ncol(values)
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    value <- if (margin == 1L) values[i, ] else values[, i]
    z <- .lt_rate_diag_vector(value, neutral)
    raw <- if (is.null(raw_values)) {
      .lt_rate_diag_raw_vector(numeric())
    } else {
      .lt_rate_diag_raw_vector(if (margin == 1L) raw_values[i, ] else raw_values[, i])
    }
    mu <- if (is.null(baseline_mu)) numeric() else {
      if (margin == 1L) baseline_mu[i, ] else baseline_mu[, i]
    }
    mu <- mu[is.finite(mu)]
    encoded <- if (is.null(zero_encoded)) logical() else {
      if (margin == 1L) zero_encoded[i, ] else zero_encoded[, i]
    }
    encoded <- encoded %in% TRUE
    rows[[i]] <- c(
      z, raw,
      zero_encoded_count = sum(encoded),
      zero_encoded_fraction_available = if (z[["available_count"]] > 0)
        sum(encoded) / z[["available_count"]] else NA_real_,
      baseline_mu_count = length(mu),
      baseline_mu_min = if (length(mu)) min(mu) else NA_real_,
      baseline_mu_median = if (length(mu)) stats::median(mu) else NA_real_,
      baseline_mu_mean = if (length(mu)) mean(mu) else NA_real_,
      baseline_mu_max = if (length(mu)) max(mu) else NA_real_
    )
  }
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(out) <- make.unique(names(out))
  out <- cbind(scientific_key = keys, out, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}

.lt_rate_diag_reason_counts <- function(rate) {
  reason <- rate$representation_reason_levels[rate$representation_reason]
  tab <- table(factor(reason, levels = rate$representation_reason_levels))
  data.frame(
    representation_reason = names(tab), count = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

.lt_rate_diag_exact_repetition <- function(values) {
  value <- values[is.finite(values)]
  if (!length(value)) {
    return(data.frame(
      extreme = c("minimum", "maximum"), value = NA_real_,
      exact_count = 0L, exact_fraction_available = NA_real_,
      interpretation = "exact_extreme_repetition_no_upstream_cause_inferred",
      stringsAsFactors = FALSE
    ))
  }
  extrema <- c(min(value), max(value))
  data.frame(
    extreme = c("minimum", "maximum"), value = extrema,
    exact_count = as.integer(vapply(extrema, function(z) sum(value == z), integer(1L))),
    exact_fraction_available = vapply(extrema, function(z) mean(value == z), numeric(1L)),
    interpretation = "exact_extreme_repetition_no_upstream_cause_inferred",
    stringsAsFactors = FALSE
  )
}

.lt_rate_diag_extremes <- function(rate, neutral, n) {
  available <- which(rate$representation_availability)
  if (!length(available) || n == 0L) {
    return(data.frame(
      rank_type = character(), rank = integer(), linear_index = integer(),
      gene_id = character(), branch_id = character(), value = double(),
      departure = double(), direction = character(), stringsAsFactors = FALSE
    ))
  }
  value <- rate$values[available]
  departure <- abs(value - neutral)
  nr <- nrow(rate$values)
  take <- function(score, rank_type, decreasing) {
    ord <- order(score, available, decreasing = decreasing, method = "radix")
    linear <- available[ord[seq_len(min(n, length(ord)))]]
    v <- rate$values[linear]
    data.frame(
      rank_type = rank_type, rank = seq_along(linear),
      linear_index = linear,
      gene_id = rate$ordered_gene_ledger[((linear - 1L) %% nr) + 1L],
      branch_id = rate$ordered_branch_ledger[((linear - 1L) %/% nr) + 1L],
      value = v, departure = abs(v - neutral),
      direction = ifelse(v > neutral, "above_neutral",
                         ifelse(v < neutral, "below_neutral", "at_neutral")),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    take(value, "highest_value", TRUE),
    take(value, "lowest_value", FALSE),
    take(departure, "largest_absolute_departure", TRUE)
  )
}

.lt_rate_diag_decomposition <- function(rate, source_matrix, baseline_fit,
                                        measurement_fit_eligibility) {
  if (is.null(source_matrix) && is.null(baseline_fit)) {
    return(list(
      compatibility = list(state = "NOT_REQUESTED", current = NA,
                           reasons = character(), repair_attempted = FALSE),
      raw_values = NULL, baseline_mu = NULL
    ))
  }
  if (is.null(source_matrix) || is.null(baseline_fit)) {
    .lt_abort("`source_matrix` and `baseline_fit` must be supplied together.")
  }
  if (is.null(measurement_fit_eligibility)) {
    measurement_fit_eligibility <- baseline_fit$domain$eligibility
  }
  status <- lt_rate_dependency_status(
    rate, source_matrix, baseline_fit, measurement_fit_eligibility
  )
  if (!isTRUE(status$current)) {
    .lt_abort(paste0(
      "STALE_RECOMPUTE_REQUIRED: diagnostic decomposition refused; ",
      paste(status$reasons, collapse = "; "), "; no automatic repair was attempted."
    ))
  }
  mu <- matrix(NA_real_, nrow(source_matrix$values), ncol(source_matrix$values),
               dimnames = dimnames(source_matrix$values))
  if (inherits(baseline_fit, "lt_c3_fit")) {
    mu[baseline_fit$edge_linear] <- baseline_fit$mu
  } else if (inherits(baseline_fit, "lt_c2_fit")) {
    mu[baseline_fit$evaluation_linear] <- baseline_fit$evaluation_mu
  } else {
    .lt_abort("`baseline_fit` must be lt_c3_fit or lt_c2_fit.")
  }
  raw <- source_matrix$values
  raw[measurement_fit_eligibility$status != "eligible"] <- NA_real_
  list(
    compatibility = unclass(status), raw_values = raw, baseline_mu = mu
  )
}

.lt_rate_diagnostic_semantic_payload <- function(x) {
  x[c(
    "schema_version", "representation", "baseline_estimand", "neutral_value",
    "ordered_gene_ledger", "ordered_branch_ledger", "source_rate_sha256",
    "native_summary", "reason_counts", "sign_counts", "exact_repetition",
    "departure_concentration", "zero_loss", "gene_summary", "branch_summary",
    "extreme_ledger", "compatibility", "provenance"
  )]
}

#' Diagnose trait-free rate geometry without changing any cell
#'
#' This threshold-free operator summarizes a materialized [lt_rate()] object.
#' It does not trim, winsorize, exclude, classify outliers, inspect traits, or
#' infer an upstream cause for repeated extrema. Optional source/baseline
#' decomposition is allowed only when both objects are supplied and the rate's
#' scientific dependencies are still current. Historical finite-sentinel
#' `logGBI` diagnostics are descriptive/read-migrate-only evidence and never
#' authorize ordinary finite-sentinel inference.
#'
#' @param rate A valid `lt_rate`.
#' @param source_matrix Optional compatible `lt_matrix`; must be paired with
#'   `baseline_fit`.
#' @param baseline_fit Optional compatible `lt_c3_fit` or `lt_c2_fit`.
#' @param measurement_fit_eligibility Optional exact Layer-2 declaration.
#' @param extreme_n Nonnegative number of cells retained for each bounded
#'   descriptive extreme ranking.
#' @return An `lt_rate_diagnostic` object.
#' @export
lt_diagnose_rate <- function(
    rate, source_matrix = NULL, baseline_fit = NULL,
    measurement_fit_eligibility = NULL, extreme_n = 20L) {
  validate_lt_rate(rate)
  if (!is.numeric(extreme_n) || length(extreme_n) != 1L || is.na(extreme_n) ||
      extreme_n < 0 || extreme_n != as.integer(extreme_n)) {
    .lt_abort("`extreme_n` must be one nonnegative integer.")
  }
  extreme_n <- as.integer(extreme_n)
  neutral <- .lt_rate_neutral_value(rate$representation)
  decomposition <- .lt_rate_diag_decomposition(
    rate, source_matrix, baseline_fit, measurement_fit_eligibility
  )
  value <- rate$values[rate$representation_availability]
  departure <- abs(value - neutral)
  native <- as.data.frame(as.list(.lt_rate_diag_vector(value, neutral)),
                          stringsAsFactors = FALSE)
  reasons <- .lt_rate_diag_reason_counts(rate)
  zero_input <- reasons$count[reasons$representation_reason ==
                                "input_zero_log_undefined"]
  boundary_zero <- reasons$count[reasons$representation_reason == "baseline_zero"]
  out <- structure(list(
    schema_version = "LT000A_T1B_LOGZERO_PROD001_v2",
    representation = rate$representation,
    baseline_estimand = rate$baseline_estimand,
    neutral_value = neutral,
    ordered_gene_ledger = rate$ordered_gene_ledger,
    ordered_branch_ledger = rate$ordered_branch_ledger,
    source_rate_sha256 = rate$output_hashes$semantic_object_sha256,
    native_summary = native,
    reason_counts = reasons,
    sign_counts = data.frame(
      sign = c("below_neutral", "at_neutral", "above_neutral"),
      count = c(sum(value < neutral), sum(value == neutral), sum(value > neutral)),
      stringsAsFactors = FALSE
    ),
    exact_repetition = .lt_rate_diag_exact_repetition(rate$values),
    departure_concentration = as.data.frame(
      as.list(.lt_rate_diag_concentration(departure)), stringsAsFactors = FALSE
    ),
    zero_loss = data.frame(
      available_count = length(value),
      zero_encoded_count = if (is.null(rate$zero_encoded)) 0L else
        as.integer(sum(rate$zero_encoded, na.rm = TRUE)),
      zero_encoded_fraction_available = if (length(value)) {
        if (is.null(rate$zero_encoded)) 0 else
          sum(rate$zero_encoded, na.rm = TRUE) / length(value)
      } else NA_real_,
      input_zero_log_undefined_count = zero_input,
      baseline_zero_count = boundary_zero,
      unavailable_count = length(rate$values) - length(value),
      interpretation = "availability_loss_is_not_bad_data_removal",
      stringsAsFactors = FALSE
    ),
    gene_summary = .lt_rate_diag_group(
      rate$values, neutral, rate$ordered_gene_ledger, 1L,
      decomposition$raw_values, decomposition$baseline_mu,
      rate$zero_encoded %||% NULL
    ),
    branch_summary = .lt_rate_diag_group(
      rate$values, neutral, rate$ordered_branch_ledger, 2L,
      decomposition$raw_values, decomposition$baseline_mu,
      rate$zero_encoded %||% NULL
    ),
    extreme_ledger = .lt_rate_diag_extremes(rate, neutral, extreme_n),
    compatibility = decomposition$compatibility,
    provenance = list(
      operation = "trait_free_threshold_free_rate_diagnosis",
      diagnose_only = TRUE, exclusion_count = 0L,
      trimming = FALSE, winsorization = FALSE,
      outlier_classification = FALSE, normality_rule = FALSE,
      trait_consulted = FALSE, rer_consulted = FALSE,
      upstream_cause_inferred = FALSE,
      zero_state_is_pseudocount = FALSE,
      zero_state_inferred_upstream_cause = FALSE,
      zero_mask_sha256 = if (is.null(rate$zero_encoded)) NA_character_ else
        .lt_hash(rate$zero_encoded),
      unusual_is_error = FALSE, extreme_is_artifact = FALSE,
      long_branch_is_anomalous = FALSE,
      bounded_extreme_n = extreme_n,
      performance_strategy = paste(
        "matrix-native summaries; margin-wise bounded loops;",
        "no full long-table materialization"
      )
    ),
    metadata = list(scientific_winner = NA_character_),
    semantic_diagnostic_sha256 = NULL
  ), class = "lt_rate_diagnostic")
  out$semantic_diagnostic_sha256 <- .lt_hash(
    .lt_rate_diagnostic_semantic_payload(out)
  )
  validate_lt_rate_diagnostic(out)
  out
}

#' Validate a trait-free rate diagnostic
#'
#' Validation is structural/scientific and does not require a registered file
#' SHA or frozen review-package identity. Mutable ordinary metadata is outside
#' the semantic diagnostic hash.
#'
#' @param x An `lt_rate_diagnostic`.
#' @return `x`, invisibly.
#' @export
validate_lt_rate_diagnostic <- function(x) {
  if (!inherits(x, "lt_rate_diagnostic")) {
    .lt_abort("`x` must inherit from lt_rate_diagnostic.")
  }
  if (!x$schema_version %in% c(
    "LT000A_T1C_RATEDIAG001_v1", "LT000A_T1B_LOGZERO_PROD001_v2"
  )) {
    .lt_abort("Unknown rate-diagnostic schema version.")
  }
  .lt_assert_scalar_character(x$representation, "representation")
  .lt_assert_scalar_character(x$baseline_estimand, "baseline_estimand")
  .lt_assert_unique_ids(x$ordered_gene_ledger, "ordered_gene_ledger")
  .lt_assert_unique_ids(x$ordered_branch_ledger, "ordered_branch_ledger")
  .lt_assert_sha256(x$source_rate_sha256, "source_rate_sha256")
  .lt_assert_sha256(x$semantic_diagnostic_sha256,
                    "semantic_diagnostic_sha256")
  if (!is.numeric(x$neutral_value) || length(x$neutral_value) != 1L ||
      is.na(x$neutral_value) ||
      !identical(x$neutral_value, .lt_rate_neutral_value(x$representation))) {
    .lt_abort("Rate diagnostic neutral value is inconsistent with representation.")
  }
  if (!is.data.frame(x$gene_summary) || !is.data.frame(x$branch_summary) ||
      !identical(x$gene_summary$scientific_key, x$ordered_gene_ledger) ||
      !identical(x$branch_summary$scientific_key, x$ordered_branch_ledger)) {
    .lt_abort("Rate diagnostic margin summaries must preserve exact scientific-key order.")
  }
  required_extreme <- c(
    "rank_type", "rank", "linear_index", "gene_id", "branch_id", "value",
    "departure", "direction"
  )
  if (!is.data.frame(x$extreme_ledger) ||
      any(!required_extreme %in% names(x$extreme_ledger)) ||
      anyDuplicated(x$extreme_ledger[c("rank_type", "rank")])) {
    .lt_abort("Rate diagnostic extreme ledger is not uniquely keyed.")
  }
  if (nrow(x$extreme_ledger) &&
      (any(!x$extreme_ledger$gene_id %in% x$ordered_gene_ledger) ||
       any(!x$extreme_ledger$branch_id %in% x$ordered_branch_ledger))) {
    .lt_abort("Rate diagnostic extreme ledger contains unknown scientific keys.")
  }
  if (!isTRUE(x$provenance$diagnose_only) ||
      !identical(x$provenance$exclusion_count, 0L) ||
      isTRUE(x$provenance$trimming) || isTRUE(x$provenance$winsorization) ||
      isTRUE(x$provenance$outlier_classification) ||
      isTRUE(x$provenance$normality_rule) ||
      isTRUE(x$provenance$trait_consulted) ||
      isTRUE(x$provenance$rer_consulted) ||
      isTRUE(x$provenance$upstream_cause_inferred)) {
    .lt_abort("Rate diagnostic must remain threshold-free diagnose-only evidence.")
  }
  if (identical(x$schema_version, "LT000A_T1B_LOGZERO_PROD001_v2")) {
    required_zero <- c(
      "zero_encoded_count", "zero_encoded_fraction_available"
    )
    if (any(!required_zero %in% names(x$zero_loss)) ||
        any(!required_zero %in% names(x$gene_summary)) ||
        any(!required_zero %in% names(x$branch_summary)) ||
        !identical(x$provenance$zero_state_is_pseudocount, FALSE) ||
        !identical(x$provenance$zero_state_inferred_upstream_cause, FALSE)) {
      .lt_abort("Production zero-state diagnostic fields are incomplete.")
    }
    encoded_count <- x$zero_loss$zero_encoded_count[[1L]]
    available_count <- x$zero_loss$available_count[[1L]]
    expected_fraction <- if (available_count > 0)
      encoded_count / available_count else NA_real_
    if (!identical(x$zero_loss$zero_encoded_fraction_available[[1L]],
                   expected_fraction) ||
        sum(x$gene_summary$zero_encoded_count) != encoded_count ||
        sum(x$branch_summary$zero_encoded_count) != encoded_count) {
      .lt_abort("Production zero-state diagnostic counts are inconsistent.")
    }
  }
  if (!x$compatibility$state %in% c("NOT_REQUESTED", "CURRENT") ||
      isTRUE(x$compatibility$repair_attempted)) {
    .lt_abort("Rate diagnostic compatibility state is invalid or implies repair.")
  }
  if (!identical(x$semantic_diagnostic_sha256,
                 .lt_hash(.lt_rate_diagnostic_semantic_payload(x)))) {
    .lt_abort("Rate diagnostic semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

#' @export
print.lt_rate_diagnostic <- function(x, ...) {
  validate_lt_rate_diagnostic(x)
  cat("<lt_rate_diagnostic> ", x$representation, "\n", sep = "")
  cat("  available cells: ", x$native_summary$available_count, "\n", sep = "")
  cat("  genes x branches: ", length(x$ordered_gene_ledger), " x ",
      length(x$ordered_branch_ledger), "\n", sep = "")
  cat("  mode: trait-free, threshold-free, diagnose-only\n")
  invisible(x)
}
