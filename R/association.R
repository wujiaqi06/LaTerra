# Trait-free-of-rate univariate branch-wise association.

.lt_association_reason_levels <- function() {
  c("eligible", "rate_unavailable", "trait_unavailable",
    "branch_user_excluded")
}

.lt_association_status_levels <- function() {
  c(
    "ESTIMABLE", "NO_ELIGIBLE_BRANCHES", "INSUFFICIENT_BRANCHES",
    "NO_TRAIT_VARIATION", "BINARY_REFERENCE_ABSENT",
    "BINARY_FOCAL_ABSENT", "NUMERICALLY_INDETERMINATE"
  )
}

.lt_rate_association_identity <- function(rate) {
  .lt_refuse_legacy_encoded_inference(
    rate, ".lt_rate_association_identity"
  )
  .lt_refuse_historical_logGBI_inference(
    rate, ".lt_rate_association_identity"
  )
  validate_lt_rate(rate)
  data.frame(
    representation = rate$representation,
    semantic_object_sha256 = rate$output_hashes$semantic_object_sha256,
    values_sha256 = rate$output_hashes$values_sha256,
    available_cell_sha256 = rate$ordered_cell_identity$available_cell_sha256,
    stringsAsFactors = FALSE
  )
}

.lt_association_domain_payload <- function(x) {
  list(
    domain_id = x$domain_id,
    gene_ids = x$gene_ids,
    branch_ids = x$branch_ids,
    domain_mode = x$domain_mode,
    rate_identities = x$rate_identities,
    branch_state_identity = x$branch_state_identity,
    eligible = x$eligible,
    reason = x$reason,
    reason_levels = x$reason_levels,
    state_index = x$state_index,
    user_branch_allowed = x$user_branch_allowed,
    ordered_gene_identity = x$ordered_gene_identity,
    ordered_branch_identity = x$ordered_branch_identity,
    cell_domain_identity = x$cell_domain_identity,
    reason_counts = x$reason_counts,
    rate_branch_count = x$rate_branch_count,
    state_branch_count = x$state_branch_count,
    matched_branch_count = x$matched_branch_count,
    rate_only_branch_ids = x$rate_only_branch_ids,
    state_only_branch_ids = x$state_only_branch_ids,
    match_fraction = x$match_fraction,
    provenance = x$provenance
  )
}

.lt_domain_rates <- function(rates, branch_state, branch_ids, mismatch,
                             domain_id, mode) {
  mismatch <- match.arg(mismatch, c("warn", "error"))
  if (is.list(rates) && length(rates)) {
    invisible(lapply(rates, .lt_refuse_legacy_encoded_inference,
                     route = ".lt_domain_rates"))
  }
  if (!is.list(rates) || !length(rates) ||
      !all(vapply(rates, inherits, logical(1L), "lt_rate"))) {
    .lt_abort("Association-domain construction requires lt_rate objects.")
  }
  invisible(lapply(rates, .lt_refuse_historical_logGBI_inference,
                   route = ".lt_domain_rates"))
  invisible(lapply(rates, validate_lt_rate))
  validate_lt_branch_state(branch_state)
  representations <- vapply(rates, `[[`, character(1L), "representation")
  if (anyDuplicated(representations)) {
    .lt_abort("Association-domain rate representations must be unique.")
  }
  if (length(rates) > 1L) {
    .lt_rate_compatible_axes(rates, require_same_baseline = TRUE)
  }
  first <- rates[[1L]]
  gene_ids <- first$ordered_gene_ledger
  rate_branch_ids <- first$ordered_branch_ledger
  rate_only <- setdiff(rate_branch_ids, branch_state$branch_ids)
  state_only <- setdiff(branch_state$branch_ids, rate_branch_ids)
  if (length(rate_only) || length(state_only)) {
    message <- paste0(
      "Exact branch-key mismatch: ", length(rate_only), " rate-only and ",
      length(state_only), " state-only branches; proceeding on exact intersection."
    )
    if (mismatch == "error") .lt_abort(message)
    warning(message, call. = FALSE)
  }
  state_index <- match(rate_branch_ids, branch_state$branch_ids)
  matched <- !is.na(state_index)
  state_available <- rep(FALSE, length(rate_branch_ids))
  state_available[matched] <-
    branch_state$availability[state_index[matched]] == "available"

  if (is.null(branch_ids)) {
    user_allowed <- rep(TRUE, length(rate_branch_ids))
    branch_ids_requested <- rate_branch_ids
  } else {
    .lt_assert_unique_ids(branch_ids, "branch_ids")
    unknown <- setdiff(branch_ids, rate_branch_ids)
    if (length(unknown)) {
      .lt_abort(paste0(
        "User analysis-domain branch IDs are not rate scientific keys: ",
        paste(utils::head(unknown, 5L), collapse = ", ")
      ))
    }
    user_allowed <- rate_branch_ids %in% branch_ids
    branch_ids_requested <- branch_ids
  }

  rate_available <- Reduce(`&`, lapply(rates, `[[`,
                                        "representation_availability"))
  state_matrix <- matrix(
    state_available, nrow = length(gene_ids), ncol = length(rate_branch_ids),
    byrow = TRUE, dimnames = list(gene_ids, rate_branch_ids)
  )
  user_matrix <- matrix(
    user_allowed, nrow = length(gene_ids), ncol = length(rate_branch_ids),
    byrow = TRUE, dimnames = list(gene_ids, rate_branch_ids)
  )
  eligible <- rate_available & state_matrix & user_matrix
  reason_label <- matrix(
    "rate_unavailable", nrow = length(gene_ids), ncol = length(rate_branch_ids),
    dimnames = list(gene_ids, rate_branch_ids)
  )
  reason_label[rate_available & !state_matrix] <- "trait_unavailable"
  reason_label[rate_available & state_matrix & !user_matrix] <-
    "branch_user_excluded"
  reason_label[eligible] <- "eligible"
  reason_levels <- .lt_association_reason_levels()
  reason <- matrix(
    as.integer(match(reason_label, reason_levels)),
    nrow = nrow(reason_label), ncol = ncol(reason_label),
    dimnames = dimnames(reason_label)
  )
  reason_counts <- data.frame(
    reason = reason_levels,
    count = tabulate(reason, nbins = length(reason_levels)),
    stringsAsFactors = FALSE
  )
  reason_counts <- reason_counts[reason_counts$count > 0L, , drop = FALSE]
  rate_identities <- do.call(rbind, lapply(rates, .lt_rate_association_identity))
  rownames(rate_identities) <- NULL
  provenance <- list(
    operator = "LaTerra_T2A_exact_keyed_layer5_domain",
    rate_construction_trait_independent = TRUE,
    exact_key_matching = TRUE,
    position_matching = FALSE,
    alias_guessing = FALSE,
    mismatch_policy = mismatch,
    user_branch_restriction_supplied = !is.null(branch_ids),
    user_branch_ids_requested = branch_ids_requested
  )
  core <- list(
    domain_id = domain_id,
    gene_ids = gene_ids,
    branch_ids = rate_branch_ids,
    domain_mode = mode,
    rate_identities = rate_identities,
    branch_state_identity = branch_state$branch_state_identity,
    eligible = eligible,
    reason = reason,
    reason_levels = reason_levels,
    state_index = as.integer(state_index),
    user_branch_allowed = user_allowed,
    ordered_gene_identity = .lt_hash(gene_ids),
    ordered_branch_identity = .lt_hash(rate_branch_ids),
    cell_domain_identity = .lt_rate_keyed_mask_hash(
      eligible, gene_ids, rate_branch_ids
    ),
    reason_counts = reason_counts,
    rate_branch_count = length(rate_branch_ids),
    state_branch_count = length(branch_state$branch_ids),
    matched_branch_count = sum(matched),
    rate_only_branch_ids = rate_only,
    state_only_branch_ids = state_only,
    match_fraction = sum(matched) / length(rate_branch_ids),
    provenance = provenance
  )
  out <- structure(c(core, list(
    association_domain_identity = .lt_hash(core),
    metadata = list()
  )), class = "lt_association_domain")
  validate_lt_association_domain(out)
  out
}

#' Construct a native association domain
#'
#' Historical finite-sentinel `lt_rate(representation = "logGBI")` inputs are
#' read/migrate-only and fail closed before domain construction.
#'
#' @param rate One `lt_rate` representation.
#' @param branch_state An `lt_branch_state`.
#' @param branch_ids Optional exact scientific-key restriction.
#' @param mismatch `warn` to use the exact intersection or `error` to require
#'   equal branch-key sets.
#' @param domain_id Stable domain identifier.
#' @return An `lt_association_domain` labelled `NATIVE_ASSOCIATION_DOMAIN`.
#' @export
lt_association_domain <- function(rate, branch_state, branch_ids = NULL,
                                  mismatch = c("warn", "error"),
                                  domain_id = paste0(
                                    "native_", rate$representation, "_",
                                    branch_state$state_id)) {
  .lt_domain_rates(
    list(rate), branch_state, branch_ids, mismatch, domain_id,
    "NATIVE_ASSOCIATION_DOMAIN"
  )
}

#' Construct an exact common association domain
#'
#' A historical finite-sentinel `logGBI` member makes the whole requested
#' common domain fail closed. It is never inferred to be an explicit replay
#' view from its values.
#'
#' @param rates A list of at least two compatible `lt_rate` representations.
#' @inheritParams lt_association_domain
#' @return An `lt_association_domain` labelled `COMMON_ASSOCIATION_DOMAIN`.
#' @export
lt_common_association_domain <- function(rates, branch_state,
                                         branch_ids = NULL,
                                         mismatch = c("warn", "error"),
                                         domain_id = paste0(
                                           "common_", branch_state$state_id)) {
  if (!is.list(rates) || length(rates) < 2L) {
    .lt_abort("Common association domain requires a list of at least two rates.")
  }
  .lt_domain_rates(
    rates, branch_state, branch_ids, mismatch, domain_id,
    "COMMON_ASSOCIATION_DOMAIN"
  )
}

#' Validate a Layer-5 association domain
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_association_domain <- function(x) {
  if (!inherits(x, "lt_association_domain")) {
    .lt_abort("`x` must inherit from `lt_association_domain`.")
  }
  .lt_assert_scalar_character(x$domain_id, "domain_id")
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  .lt_assert_unique_ids(x$branch_ids, "branch_ids")
  if (!x$domain_mode %in% c("NATIVE_ASSOCIATION_DOMAIN",
                            "COMMON_ASSOCIATION_DOMAIN")) {
    .lt_abort("Unknown association domain mode.")
  }
  if (!is.data.frame(x$rate_identities) || !identical(
    names(x$rate_identities),
    c("representation", "semantic_object_sha256", "values_sha256",
      "available_cell_sha256")
  ) || anyDuplicated(x$rate_identities$representation)) {
    .lt_abort("Association-domain rate identities are malformed or ambiguous.")
  }
  hash_columns <- c("semantic_object_sha256", "values_sha256",
                    "available_cell_sha256")
  if (any(!vapply(x$rate_identities[hash_columns], function(z)
    all(grepl("^[0-9a-f]{64}$", z)), logical(1L)))) {
    .lt_abort("Association-domain rate identity hashes are malformed.")
  }
  dn <- list(x$gene_ids, x$branch_ids)
  if (!is.logical(x$eligible) || anyNA(x$eligible) ||
      !identical(dimnames(x$eligible), dn) ||
      !is.integer(x$reason) || anyNA(x$reason) ||
      !identical(dimnames(x$reason), dn) ||
      any(!x$reason %in% seq_along(x$reason_levels)) ||
      !identical(x$reason_levels, .lt_association_reason_levels())) {
    .lt_abort("Association-domain eligibility/reason layers are malformed.")
  }
  labels <- x$reason_levels[x$reason]
  if (any(x$eligible & labels != "eligible") ||
      any(!x$eligible & labels == "eligible")) {
    .lt_abort("Association-domain eligibility and reason codes disagree.")
  }
  if (!is.integer(x$state_index) || length(x$state_index) != length(x$branch_ids) ||
      any(x$state_index < 1L, na.rm = TRUE) ||
      !is.logical(x$user_branch_allowed) || anyNA(x$user_branch_allowed) ||
      length(x$user_branch_allowed) != length(x$branch_ids)) {
    .lt_abort("Association-domain keyed state/user-branch mapping is malformed.")
  }
  if (!identical(x$ordered_gene_identity, .lt_hash(x$gene_ids)) ||
      !identical(x$ordered_branch_identity, .lt_hash(x$branch_ids)) ||
      !identical(x$cell_domain_identity, .lt_rate_keyed_mask_hash(
        x$eligible, x$gene_ids, x$branch_ids
      ))) {
    .lt_abort("Association-domain ordered or cell identity is stale.")
  }
  expected_counts <- data.frame(
    reason = x$reason_levels,
    count = tabulate(x$reason, nbins = length(x$reason_levels)),
    stringsAsFactors = FALSE
  )
  expected_counts <- expected_counts[expected_counts$count > 0L, , drop = FALSE]
  rownames(expected_counts) <- NULL
  observed_counts <- x$reason_counts
  rownames(observed_counts) <- NULL
  if (!identical(observed_counts, expected_counts)) {
    .lt_abort("Association-domain reason counts are stale.")
  }
  .lt_assert_sha256(x$branch_state_identity, "branch_state_identity")
  .lt_assert_named_list(x$provenance, "provenance")
  .lt_assert_named_list(x$metadata, "metadata")
  expected_identity <- .lt_hash(.lt_association_domain_payload(x))
  .lt_assert_sha256(x$association_domain_identity,
                    "association_domain_identity")
  if (!identical(x$association_domain_identity, expected_identity)) {
    .lt_abort("Association-domain semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

.lt_association_payload <- function(x) {
  x[c(
    "association_id", "rate_representation", "rate_identity",
    "branch_state_identity", "association_domain_identity", "trait_id",
    "branch_state_id", "trait_type", "estimand", "effect_direction",
    "gene_ids", "status", "reason", "estimates", "descriptive_scores",
    "domain_counts", "inference_status", "significance_calibration",
    "rate_neutral_point", "rate_units", "provenance"
  )]
}

.lt_safe_cor <- function(x, y, method) {
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  suppressWarnings(stats::cor(x, y, method = method))
}

#' Compute per-gene univariate branch-wise association point estimates
#'
#' Rate is the response and branch state is the predictor. The operator is
#' unweighted and exposes no p-value, confidence interval, FDR, or prediction
#' claim. Rate construction remains trait independent. Historical finite
#' `lt_rate(representation = "logGBI")` inputs are forbidden. The separately
#' typed finite encoded view is diagnostic/replay-only and is also refused by
#' ordinary inference. Use GBI for La Terra V1 inference.
#'
#' @param rate One certified `lt_rate`.
#' @param branch_state One declared `lt_branch_state`.
#' @param domain Optional compatible `lt_association_domain`; a native domain
#'   is constructed when omitted.
#' @param mismatch Passed to native domain construction when `domain` is `NULL`.
#' @param association_id Stable result identifier.
#' @param encoded_mode Deprecated compatibility argument. It cannot authorize
#'   finite-encoded inference.
#' @param acknowledge_encoded_zero Deprecated compatibility argument. Its
#'   value has no effect on the mandatory refusal.
#' @return An `lt_association` containing every gene, including explicit
#'   non-estimability states.
#' @export
lt_associate_univariate <- function(
    rate, branch_state, domain = NULL, mismatch = c("warn", "error"),
    association_id = NULL,
    encoded_mode = NULL, acknowledge_encoded_zero = FALSE) {
  .lt_refuse_legacy_encoded_inference(rate, "lt_associate_univariate")
  if (is.null(association_id)) {
    label <- if (inherits(rate, "lt_logGBI_encoded")) rate$mode else
      if (inherits(rate, "lt_log_relative")) rate$representation_id else
        rate$representation
    association_id <- paste0(label, "_", branch_state$state_id,
                             "_univariate")
  }
  if (inherits(rate, "lt_log_relative")) {
    .lt_abort_class(
      "Use `lt_associate_log_relative()` and select separate ZERO_MASS/POSITIVE_LOG components.",
      "lt_error_log_relative_requires_component_association"
    )
  }
  .lt_refuse_historical_logGBI_inference(rate, "lt_associate_univariate")
  validate_lt_rate(rate)
  validate_lt_branch_state(branch_state)
  if (!branch_state$type %in% c("binary", "continuous")) {
    .lt_abort("T2A univariate association supports only binary or continuous branch states.")
  }
  if (is.null(domain)) {
    domain <- lt_association_domain(rate, branch_state, mismatch = mismatch)
  }
  validate_lt_association_domain(domain)
  .lt_assert_scalar_character(association_id, "association_id")
  if (!identical(rate$ordered_gene_ledger, domain$gene_ids) ||
      !identical(rate$ordered_branch_ledger, domain$branch_ids)) {
    .lt_abort("Association rate and domain ordered scientific keys differ.")
  }
  identity <- .lt_rate_association_identity(rate)
  match_identity <- domain$rate_identities$representation == rate$representation &
    domain$rate_identities$semantic_object_sha256 == identity$semantic_object_sha256
  if (sum(match_identity) != 1L) {
    .lt_abort("Association rate is not an exact scientific dependency of the domain.")
  }
  if (!identical(branch_state$branch_state_identity,
                 domain$branch_state_identity)) {
    .lt_abort("Association branch state differs from the domain dependency.")
  }
  state_values <- rep(NA, length(domain$branch_ids))
  matched <- !is.na(domain$state_index)
  state_values[matched] <- branch_state$values[domain$state_index[matched]]
  n_gene <- length(domain$gene_ids)
  status <- character(n_gene)
  reason <- character(n_gene)
  intercept_value <- beta_value <- reference_mean_value <-
    focal_mean_value <- pearson_value <- spearman_value <-
    r_squared_value <- trait_mean_value <- trait_sd_value <-
    rate_mean_value <- rate_sd_value <- rep(NA_real_, n_gene)
  n_reference_value <- n_focal_value <- rep(NA_integer_, n_gene)
  n_value <- integer(n_gene)
  small_domain_value <- logical(n_gene)
  ref_key <- focal_key <- NA_character_
  binary_stats <- NULL
  if (branch_state$type == "binary") {
    ref_key <- .lt_level_key(branch_state$coding$reference_level)
    focal_key <- .lt_level_key(branch_state$coding$focal_level)
    state_keys <- rep(NA_character_, length(state_values))
    present_state <- !is.na(state_values)
    state_keys[present_state] <- vapply(
      state_values[present_state], .lt_level_key, character(1L)
    )
    summarize_group <- function(columns) {
      group_values <- rate$values[, columns, drop = FALSE]
      group_eligible <- domain$eligible[, columns, drop = FALSE]
      group_n <- rowSums(group_eligible)
      group_values[!group_eligible] <- 0
      group_sum <- rowSums(group_values)
      group_sum_squares <- rowSums(group_values * group_values)
      list(n = as.integer(group_n), sum = group_sum,
           sum_squares = group_sum_squares)
    }
    reference <- summarize_group(which(state_keys == ref_key))
    focal <- summarize_group(which(state_keys == focal_key))
    total_n <- reference$n + focal$n
    total_sum <- reference$sum + focal$sum
    total_sum_squares <- reference$sum_squares + focal$sum_squares
    total_ss <- total_sum_squares - total_sum^2 / pmax(total_n, 1L)
    total_ss[total_ss < 0 & total_ss > -1e-12 *
      pmax(total_sum_squares, 1)] <- 0
    total_ss[total_ss < 0] <- NA_real_
    reference_mean_all <- reference$sum / reference$n
    focal_mean_all <- focal$sum / focal$n
    beta_all <- focal_mean_all - reference_mean_all
    rate_mean_all <- total_sum / total_n
    rate_sd_all <- sqrt(total_ss / pmax(total_n - 1L, 1L))
    point_biserial_all <- beta_all * sqrt(reference$n * focal$n) /
      sqrt(total_n * total_ss)
    point_biserial_all[!is.finite(point_biserial_all)] <- NA_real_
    binary_stats <- list(
      n_reference = reference$n, n_focal = focal$n, n = total_n,
      reference_mean = reference_mean_all, focal_mean = focal_mean_all,
      beta = beta_all, rate_mean = rate_mean_all, rate_sd = rate_sd_all,
      point_biserial = point_biserial_all
    )
  }
  for (i in seq_len(n_gene)) {
    n_reference <- n_focal <- NA_integer_
    alpha <- beta <- pearson <- spearman <- r_squared <- NA_real_
    reference_mean <- focal_mean <- trait_mean <- trait_sd <- NA_real_
    rate_mean <- rate_sd <- NA_real_
    small <- FALSE
    if (branch_state$type == "continuous") {
      keep <- domain$eligible[i, ]
      y <- rate$values[i, keep]
      x <- state_values[keep]
      n <- length(y)
    } else {
      n_reference <- binary_stats$n_reference[[i]]
      n_focal <- binary_stats$n_focal[[i]]
      n <- binary_stats$n[[i]]
    }
    if (!n) {
      state <- "NO_ELIGIBLE_BRANCHES"
      why <- "no_cell_passes_layer5_eligibility"
    } else if (branch_state$type == "continuous") {
      x <- as.double(x); y <- as.double(y)
      trait_mean <- mean(x); trait_sd <- stats::sd(x)
      rate_mean <- mean(y); rate_sd <- stats::sd(y)
      if (n < 2L) {
        state <- "INSUFFICIENT_BRANCHES"
        why <- "continuous_slope_requires_at_least_two_observations"
      } else if (length(unique(x)) < 2L) {
        state <- "NO_TRAIT_VARIATION"
        why <- "continuous_predictor_has_no_variation"
      } else {
        denominator <- sum((x - trait_mean)^2)
        if (!is.finite(denominator) || denominator <= 0) {
          state <- "NO_TRAIT_VARIATION"
          why <- "continuous_predictor_denominator_not_positive"
        } else {
          beta <- sum((x - trait_mean) * (y - rate_mean)) / denominator
          alpha <- rate_mean - beta * trait_mean
          pearson <- .lt_safe_cor(x, y, "pearson")
          spearman <- .lt_safe_cor(x, y, "spearman")
          r_squared <- pearson^2
          if (!is.finite(beta) || !is.finite(alpha)) {
            state <- "NUMERICALLY_INDETERMINATE"
            why <- "nonfinite_continuous_point_estimate"
          } else {
            state <- "ESTIMABLE"
            why <- "point_estimate_available"
            small <- n < 3L
          }
        }
      }
    } else {
      rate_mean <- binary_stats$rate_mean[[i]]
      rate_sd <- binary_stats$rate_sd[[i]]
      if (n_reference == 0L) {
        state <- "BINARY_REFERENCE_ABSENT"
        why <- "no_eligible_reference_branch"
      } else if (n_focal == 0L) {
        state <- "BINARY_FOCAL_ABSENT"
        why <- "no_eligible_focal_branch"
      } else {
        reference_mean <- binary_stats$reference_mean[[i]]
        focal_mean <- binary_stats$focal_mean[[i]]
        alpha <- reference_mean
        beta <- binary_stats$beta[[i]]
        pearson <- binary_stats$point_biserial[[i]]
        r_squared <- pearson^2
        if (!is.finite(beta) || !is.finite(alpha)) {
          state <- "NUMERICALLY_INDETERMINATE"
          why <- "nonfinite_binary_point_estimate"
        } else {
          state <- "ESTIMABLE"
          why <- "point_estimate_available"
          small <- n_reference == 1L || n_focal == 1L
        }
      }
    }
    status[[i]] <- state; reason[[i]] <- why
    intercept_value[[i]] <- alpha
    beta_value[[i]] <- beta
    reference_mean_value[[i]] <- reference_mean
    focal_mean_value[[i]] <- focal_mean
    pearson_value[[i]] <- pearson
    spearman_value[[i]] <- spearman
    r_squared_value[[i]] <- r_squared
    trait_mean_value[[i]] <- trait_mean
    trait_sd_value[[i]] <- trait_sd
    rate_mean_value[[i]] <- rate_mean
    rate_sd_value[[i]] <- rate_sd
    n_reference_value[[i]] <- n_reference
    n_focal_value[[i]] <- n_focal
    n_value[[i]] <- n
    small_domain_value[[i]] <- small
  }
  reference_level <- focal_level <- rep(NA_character_, n_gene)
  raw_mean_difference <- rep(NA_real_, n_gene)
  if (branch_state$type == "binary") {
    reference_level[] <- as.character(branch_state$coding$reference_level)
    focal_level[] <- as.character(branch_state$coding$focal_level)
    raw_mean_difference <- beta_value
  }
  estimates <- data.frame(
    gene_id = domain$gene_ids, intercept = intercept_value,
    beta = beta_value, reference_level = reference_level,
    focal_level = focal_level, reference_mean = reference_mean_value,
    focal_mean = focal_mean_value, raw_mean_difference = raw_mean_difference,
    n_reference = n_reference_value, n_focal = n_focal_value, n = n_value,
    small_domain_warning = small_domain_value, stringsAsFactors = FALSE
  )
  descriptive_scores <- data.frame(
    gene_id = domain$gene_ids, pearson_r = pearson_value,
    spearman_rho = if (branch_state$type == "continuous")
      spearman_value else rep(NA_real_, n_gene),
    point_biserial_r = if (branch_state$type == "binary")
      pearson_value else rep(NA_real_, n_gene),
    r_squared = r_squared_value, trait_mean = trait_mean_value,
    trait_sd = trait_sd_value, rate_mean = rate_mean_value,
    rate_sd = rate_sd_value, stringsAsFactors = FALSE
  )
  domain_counts <- data.frame(
    gene_id = domain$gene_ids, eligible_count = n_value,
    eligible_terminal_count = rep(NA_integer_, n_gene),
    eligible_internal_count = rep(NA_integer_, n_gene),
    stringsAsFactors = FALSE
  )
  estimand <- if (branch_state$type == "continuous")
    "unweighted_OLS_slope_with_intercept" else
    "unweighted_focal_minus_reference_mean_difference"
  effect_direction <- if (branch_state$type == "continuous")
    "rate_response_per_declared_trait_unit" else paste0(
      as.character(branch_state$coding$focal_level), " - ",
      as.character(branch_state$coding$reference_level)
    )
  neutral <- if (rate$representation %in% c("GBI", "GBI_C2")) 1 else 0
  units <- rate$representation_provenance$units %||%
    rate$metadata$units %||% "unspecified"
  provenance <- list(
    operator = "LaTerra_T2A_univariate_branchwise_point_estimate",
    rate_is_response = TRUE,
    branch_state_is_predictor = TRUE,
    unweighted = TRUE,
    rate_transformation = "none",
    trait_transformation = "none",
    rate_construction_trait_independent = TRUE,
    association_not_prediction = TRUE,
    p_value = "not_computed",
    fdr = "not_computed"
  )
  core <- list(
    association_id = association_id,
    rate_representation = rate$representation,
    rate_identity = identity$semantic_object_sha256,
    branch_state_identity = branch_state$branch_state_identity,
    association_domain_identity = domain$association_domain_identity,
    trait_id = branch_state$trait_id,
    branch_state_id = branch_state$state_id,
    trait_type = branch_state$type,
    estimand = estimand,
    effect_direction = effect_direction,
    gene_ids = domain$gene_ids,
    status = status,
    reason = reason,
    estimates = estimates,
    descriptive_scores = descriptive_scores,
    domain_counts = domain_counts,
    inference_status = "POINT_ESTIMATE_ONLY",
    significance_calibration = "NOT_PERFORMED",
    rate_neutral_point = neutral,
    rate_units = units,
    provenance = provenance
  )
  out <- structure(c(core, list(
    association_identity = .lt_hash(core),
    metadata = list()
  )), class = "lt_association")
  validate_lt_association(out)
  out
}

#' Validate an association point-estimate result
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_association <- function(x) {
  if (!inherits(x, "lt_association")) {
    .lt_abort("`x` must inherit from `lt_association`.")
  }
  .lt_assert_scalar_character(x$association_id, "association_id")
  .lt_assert_scalar_character(x$rate_representation, "rate_representation")
  .lt_assert_sha256(x$rate_identity, "rate_identity")
  .lt_assert_sha256(x$branch_state_identity, "branch_state_identity")
  .lt_assert_sha256(x$association_domain_identity,
                    "association_domain_identity")
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  n <- length(x$gene_ids)
  if (!is.character(x$status) || length(x$status) != n || anyNA(x$status) ||
      any(!x$status %in% .lt_association_status_levels()) ||
      !is.character(x$reason) || length(x$reason) != n || anyNA(x$reason) ||
      any(!nzchar(x$reason))) {
    .lt_abort("Association gene status/reason ledger is malformed.")
  }
  for (name in c("estimates", "descriptive_scores", "domain_counts")) {
    z <- x[[name]]
    if (!is.data.frame(z) || nrow(z) != n ||
        !identical(z$gene_id, x$gene_ids)) {
      .lt_abort(paste0("Association `", name,
                       "` must retain the exact ordered gene ledger."))
    }
  }
  forbidden <- c("p_value", "pvalue", "q_value", "fdr", "confidence_interval")
  released_names <- tolower(unlist(lapply(
    x[c("estimates", "descriptive_scores", "domain_counts")], names
  )))
  if (any(released_names %in% forbidden)) {
    .lt_abort("Association result cannot expose uncalibrated inferential fields.")
  }
  if (!identical(x$inference_status, "POINT_ESTIMATE_ONLY") ||
      !identical(x$significance_calibration, "NOT_PERFORMED")) {
    .lt_abort("T2A association cannot masquerade as calibrated inference.")
  }
  .lt_assert_named_list(x$provenance, "provenance")
  .lt_assert_named_list(x$metadata, "metadata")
  expected_identity <- .lt_hash(.lt_association_payload(x))
  .lt_assert_sha256(x$association_identity, "association_identity")
  if (!identical(x$association_identity, expected_identity)) {
    .lt_abort("Association semantic identity is stale or inconsistent.")
  }
  invisible(x)
}
