# Exact keyed common domains and descriptive C3-vs-C2 sensitivity summaries.

.lt_rate_reason_counts <- function(rate) {
  reason <- rate$representation_reason_levels[rate$representation_reason]
  z <- as.data.frame(table(reason), stringsAsFactors = FALSE)
  names(z) <- c("reason", "count")
  z$count <- as.integer(z$count)
  z <- z[z$count > 0L, , drop = FALSE]
  z$representation <- rate$representation
  z[, c("representation", "reason", "count")]
}

.lt_rate_compatible_axes <- function(rates, require_same_baseline = TRUE) {
  invisible(lapply(rates, validate_lt_rate))
  first <- rates[[1L]]
  for (z in rates[-1L]) {
    if (!identical(z$ordered_gene_ledger, first$ordered_gene_ledger) ||
        !identical(z$ordered_branch_ledger, first$ordered_branch_ledger) ||
        !identical(z$coord_state, first$coord_state) ||
        !identical(z$measurement_fit_eligibility$status_sha256,
                   first$measurement_fit_eligibility$status_sha256) ||
        !identical(z$measurement_fit_eligibility$reason_sha256,
                   first$measurement_fit_eligibility$reason_sha256) ||
        !identical(z$measurement_fit_eligibility$provenance_sha256,
                   first$measurement_fit_eligibility$provenance_sha256) ||
        !identical(z$representation_provenance$source_matrix_identity_sha256,
                   first$representation_provenance$source_matrix_identity_sha256)) {
      .lt_abort("Rate objects do not share the exact keyed source/coordinate/Layer-2 identity.")
    }
    if (require_same_baseline &&
        (!identical(z$baseline_estimand, first$baseline_estimand) ||
         !identical(z$baseline_fit_id, first$baseline_fit_id) ||
         !identical(z$baseline_mu_sha256, first$baseline_mu_sha256))) {
      .lt_abort("Production common-domain construction requires one exact baseline fit identity.")
    }
  }
  invisible(TRUE)
}

.lt_common_domain_semantic_payload <- function(x) {
  list(
    representations = x$representations,
    source_domain_hashes_sha256 = x$source_domain_hashes_sha256,
    common_cell_count = x$common_cell_count,
    common_cell_sha256 = x$common_cell_sha256,
    common_domain_rule = x$common_domain_rule,
    reason_counts_sha256 = x$reason_counts_sha256,
    ordered_gene_ledger_sha256 = x$ordered_gene_ledger_sha256,
    ordered_branch_ledger_sha256 = x$ordered_branch_ledger_sha256,
    domain_loss_fraction = x$domain_loss_fraction,
    domain_loss_denominator = x$domain_loss_denominator,
    zero_cells_lost_solely_log_undefined =
      x$zero_cells_lost_solely_log_undefined,
    boundary_zero_cells_lost_ratio_denominator_zero =
      x$boundary_zero_cells_lost_ratio_denominator_zero,
    full_domain_label = x$full_domain_label,
    comparison_scope = x$comparison_scope
  )
}

#' Construct the exact keyed common domain of rate representations
#'
#' The compact result stores the common source-order linear ledger and a
#' key-order-invariant hash. Full-domain summaries are labelled
#' `UNEQUAL_DOMAIN_DESCRIPTIVE` whenever source availability differs; such a
#' difference is not estimator-performance evidence.
#'
#' @param ... `lt_rate` objects or one list of them.
#' @return An `lt_rate_common_domain`.
#' @export
lt_rate_common_domain <- function(...) {
  rates <- list(...)
  if (length(rates) == 1L && is.list(rates[[1L]]) &&
      !inherits(rates[[1L]], "lt_rate")) rates <- rates[[1L]]
  if (length(rates) < 2L || !all(vapply(rates, inherits, logical(1L), "lt_rate"))) {
    .lt_abort("Common-domain construction requires at least two lt_rate objects.")
  }
  representations <- vapply(rates, `[[`, character(1L), "representation")
  if (anyDuplicated(representations)) {
    .lt_abort("Requested common-domain representations must be unique.")
  }
  .lt_rate_compatible_axes(rates, require_same_baseline = TRUE)
  first <- rates[[1L]]
  common <- Reduce(`&`, lapply(rates, `[[`, "representation_availability"))
  union_domain <- Reduce(`|`, lapply(rates, `[[`, "representation_availability"))
  source_domains <- do.call(rbind, lapply(rates, function(z) data.frame(
    representation = z$representation,
    inclusion_count = sum(z$representation_availability),
    exclusion_count = length(z$representation_availability) -
      sum(z$representation_availability),
    available_cell_sha256 = z$ordered_cell_identity$available_cell_sha256,
    output_sha256 = z$output_hashes$output_sha256,
    stringsAsFactors = FALSE
  )))
  reason_counts <- do.call(rbind, lapply(rates, .lt_rate_reason_counts))
  reason_matrices <- lapply(rates, function(z) {
    matrix(z$representation_reason_levels[z$representation_reason],
           nrow = nrow(z$values), ncol = ncol(z$values))
  })
  names(reason_matrices) <- representations
  add_name <- intersect(representations, c("ADD_LT", "ADD_C2"))
  gbi_name <- intersect(representations, c("GBI", "GBI_C2"))
  log_name <- intersect(
    representations, c("logGBI", "logGBI_strict", "logGBI_C2")
  )
  zero_log_loss <- 0L; boundary_zero_loss <- 0L
  if (length(add_name) == 1L && length(gbi_name) == 1L && length(log_name) == 1L) {
    add <- rates[[match(add_name, representations)]]$representation_availability
    gbi <- rates[[match(gbi_name, representations)]]$representation_availability
    log_ok <- rates[[match(log_name, representations)]]$representation_availability
    log_reason <- reason_matrices[[log_name]]
    zero_log_loss <- sum(add & gbi & !log_ok &
                           log_reason == "input_zero_log_undefined")
    boundary_zero_loss <- sum(add & !gbi & !log_ok &
                                log_reason == "baseline_zero")
  }
  out <- structure(list(
    representations = representations,
    ordered_gene_ledger = first$ordered_gene_ledger,
    ordered_branch_ledger = first$ordered_branch_ledger,
    gene_display_labels = first$gene_display_labels %||%
      first$ordered_gene_ledger,
    branch_display_labels = first$branch_display_labels %||%
      first$ordered_branch_ledger,
    ordered_gene_ledger_sha256 = .lt_hash(first$ordered_gene_ledger),
    ordered_branch_ledger_sha256 = .lt_hash(first$ordered_branch_ledger),
    source_domains = source_domains,
    source_domain_hashes_sha256 = .lt_hash(source_domains),
    common_linear = which(common),
    common_cell_count = sum(common),
    common_cell_sha256 = .lt_rate_keyed_mask_hash(
      common, first$ordered_gene_ledger, first$ordered_branch_ledger
    ),
    common_domain_rule = "exact keyed intersection of requested representation availability",
    reason_counts = reason_counts,
    reason_counts_sha256 = .lt_hash(reason_counts),
    domain_loss_fraction = if (any(union_domain))
      1 - sum(common) / sum(union_domain) else NA_real_,
    domain_loss_denominator = "union_of_requested_available_domains",
    zero_cells_lost_solely_log_undefined = as.integer(zero_log_loss),
    boundary_zero_cells_lost_ratio_denominator_zero = as.integer(boundary_zero_loss),
    full_domain_label = if (all(vapply(rates[-1L], function(z) identical(
      z$representation_availability, first$representation_availability
    ), logical(1L)))) "EQUAL_DOMAIN_DESCRIPTIVE" else "UNEQUAL_DOMAIN_DESCRIPTIVE",
    comparison_scope = "EXACT_KEYED_COMMON_DOMAIN",
    scientific_pass = FALSE,
    interpretation = "domain difference is not estimator-performance evidence"
  ), class = "lt_rate_common_domain")
  out$semantic_common_domain_sha256 <- .lt_hash(
    .lt_common_domain_semantic_payload(out)
  )
  validate_lt_rate_common_domain(out)
  out
}

#' Validate a rate common-domain object
#'
#' Validation reconstructs the compact keyed mask from `common_linear` and the
#' stored scientific axes. Mutable display labels are not part of its semantic
#' identity.
#'
#' @param x An `lt_rate_common_domain`.
#' @return `x`, invisibly.
#' @export
validate_lt_rate_common_domain <- function(x) {
  if (!inherits(x, "lt_rate_common_domain")) {
    .lt_abort("`x` must inherit from lt_rate_common_domain.")
  }
  .lt_assert_unique_ids(x$representations, "representations")
  .lt_assert_unique_ids(x$ordered_gene_ledger, "ordered_gene_ledger")
  .lt_assert_unique_ids(x$ordered_branch_ledger, "ordered_branch_ledger")
  gene_display <- x$gene_display_labels %||% x$ordered_gene_ledger
  branch_display <- x$branch_display_labels %||% x$ordered_branch_ledger
  if (!is.character(gene_display) ||
      length(gene_display) != length(x$ordered_gene_ledger) || anyNA(gene_display) ||
      any(!nzchar(gene_display)) ||
      !is.character(branch_display) ||
      length(branch_display) != length(x$ordered_branch_ledger) ||
      anyNA(branch_display) || any(!nzchar(branch_display))) {
    .lt_abort("Common-domain display labels must align to scientific keys.")
  }
  total_cells <- length(x$ordered_gene_ledger) *
    length(x$ordered_branch_ledger)
  if (!is.numeric(x$common_linear) || anyNA(x$common_linear) ||
      is.unsorted(x$common_linear, strictly = TRUE) ||
      any(x$common_linear < 1L | x$common_linear > total_cells) ||
      !identical(x$common_cell_count, length(x$common_linear)) ||
      !identical(x$source_domain_hashes_sha256, .lt_hash(x$source_domains)) ||
      !identical(x$reason_counts_sha256, .lt_hash(x$reason_counts))) {
    .lt_abort("Rate common-domain ledger/count/hash contract failed.")
  }
  common <- matrix(FALSE, length(x$ordered_gene_ledger),
                   length(x$ordered_branch_ledger), dimnames = list(
                     x$ordered_gene_ledger, x$ordered_branch_ledger
                   ))
  common[x$common_linear] <- TRUE
  if (!identical(x$ordered_gene_ledger_sha256,
                 .lt_hash(x$ordered_gene_ledger)) ||
      !identical(x$ordered_branch_ledger_sha256,
                 .lt_hash(x$ordered_branch_ledger)) ||
      !identical(x$common_cell_sha256, .lt_rate_keyed_mask_hash(
        common, x$ordered_gene_ledger, x$ordered_branch_ledger
      )) ||
      !identical(x$semantic_common_domain_sha256,
                 .lt_hash(.lt_common_domain_semantic_payload(x)))) {
    .lt_abort("Rate common-domain keyed semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

.lt_paired_summary <- function(first, second, comparison, ratio = FALSE) {
  keep <- is.finite(first) & is.finite(second)
  first <- first[keep]; second <- second[keep]
  difference <- second - first
  q <- if (length(difference)) as.double(stats::quantile(
    difference, c(.05, .25, .5, .75, .95), names = FALSE, type = 7
  )) else rep(NA_real_, 5L)
  ratio_keep <- ratio & first > 0 & second > 0
  lr <- if (any(ratio_keep)) log(second[ratio_keep] / first[ratio_keep]) else numeric()
  lrq <- if (length(lr)) as.double(stats::quantile(
    lr, c(.05, .5, .95), names = FALSE, type = 7
  )) else rep(NA_real_, 3L)
  data.frame(
    comparison = comparison,
    paired_count = length(first),
    pearson_correlation = if (length(first) > 1L) stats::cor(first, second) else NA_real_,
    spearman_correlation = if (length(first) > 1L)
      stats::cor(first, second, method = "spearman") else NA_real_,
    median_paired_difference_second_minus_first = if (length(difference))
      stats::median(difference) else NA_real_,
    median_absolute_paired_difference = if (length(difference))
      stats::median(abs(difference)) else NA_real_,
    difference_q05 = q[[1L]], difference_q25 = q[[2L]],
    difference_q50 = q[[3L]], difference_q75 = q[[4L]],
    difference_q95 = q[[5L]],
    maximum_finite_absolute_difference = if (length(difference))
      max(abs(difference)) else NA_real_,
    positive_ratio_count = length(lr),
    log_ratio_q05 = lrq[[1L]], median_log_ratio = lrq[[2L]],
    log_ratio_q95 = lrq[[3L]],
    difference_direction = "second_minus_first",
    comparison_scope = "EXACT_KEYED_COMMON_DOMAIN",
    stringsAsFactors = FALSE
  )
}

#' Compare two rate representations on their exact keyed common domain
#'
#' @param first First `lt_rate`.
#' @param second Second `lt_rate`.
#' @param positive_ratio Whether positive ratios are meaningful.
#' @param comparison Stable comparison label.
#' @return One-row descriptive paired summary with domain identities.
#' @export
lt_compare_rate_pair <- function(first, second, positive_ratio = FALSE,
                                 comparison = paste(first$representation,
                                                    second$representation,
                                                    sep = "_vs_")) {
  .lt_rate_compatible_axes(list(first, second), require_same_baseline = FALSE)
  if (!is.logical(positive_ratio) || length(positive_ratio) != 1L ||
      is.na(positive_ratio)) .lt_abort("`positive_ratio` must be TRUE or FALSE.")
  .lt_assert_scalar_character(comparison, "comparison")
  common <- first$representation_availability & second$representation_availability
  out <- .lt_paired_summary(first$values[common], second$values[common],
                            comparison, ratio = positive_ratio)
  out$first_full_count <- sum(first$representation_availability)
  out$second_full_count <- sum(second$representation_availability)
  out$common_cell_sha256 <- .lt_rate_keyed_mask_hash(
    common, first$ordered_gene_ledger, first$ordered_branch_ledger
  )
  out$full_domain_label <- if (identical(
    first$representation_availability, second$representation_availability
  )) "EQUAL_DOMAIN_DESCRIPTIVE" else "UNEQUAL_DOMAIN_DESCRIPTIVE"
  out$scientific_winner <- NA_character_
  out
}

#' Exact common-positive C3-vs-C2 sensitivity comparison
#'
#' The primary paired domain requires positive Y and positive, available C3 and
#' C2 baselines. The result is descriptive sensitivity evidence only and never
#' selects a winner.
#'
#' @param c3_fit Certified production `lt_c3_fit`.
#' @param c2_fit Certified sensitivity `lt_c2_fit`.
#' @return An `lt_c3_c2_comparison`.
#' @export
lt_c3_c2_compare <- function(c3_fit, c2_fit) {
  validate_lt_c3_fit(c3_fit); validate_lt_c2_fit(c2_fit)
  if (identical(c3_fit$state, "NUMERICALLY_INDETERMINATE") ||
      identical(c2_fit$state, "NUMERICALLY_INDETERMINATE")) {
    .lt_abort("C3-vs-C2 comparison requires two certified completed fits.")
  }
  c3d <- c3_fit$domain; c2d <- c2_fit$domain
  if (!identical(c3d$matrix$gene_ids, c2d$matrix$gene_ids) ||
      !identical(c3d$matrix$branch_ids, c2d$matrix$branch_ids) ||
      !identical(c3d$hashes$semantic_values_sha256,
                 c2d$hashes$semantic_values_sha256) ||
      !identical(c3d$hashes$semantic_coordinate_state_sha256,
                 c2d$hashes$semantic_coordinate_state_sha256) ||
      !identical(c3d$hashes$measurement_eligibility_sha256,
                 c2d$hashes$measurement_eligibility_sha256) ||
      !identical(c3d$hashes$measurement_reason_sha256,
                 c2d$hashes$measurement_reason_sha256)) {
    .lt_abort("C3 and C2 fits do not share one exact keyed matrix/Layer-2 identity.")
  }
  c3_positive <- c3_fit$observed_value > 0 & c3_fit$mu > 0
  c3_linear <- c3_fit$edge_linear[c3_positive]
  c2_y <- c2d$matrix$values[c2_fit$evaluation_linear]
  c2_positive <- c2_fit$baseline_availability & c2_y > 0 &
    c2_fit$evaluation_mu > 0
  c2_linear <- c2_fit$evaluation_linear[c2_positive]
  common_linear <- sort(intersect(c3_linear, c2_linear), method = "radix")
  c3_pos <- match(common_linear, c3_fit$edge_linear)
  c2_pos <- match(common_linear, c2_fit$evaluation_linear)
  if (anyNA(c3_pos) || anyNA(c2_pos)) {
    .lt_abort("Exact keyed C3/C2 common-positive mapping failed.")
  }
  y <- c3_fit$observed_value[c3_pos]
  mu3 <- c3_fit$mu[c3_pos]
  mu2 <- c2_fit$evaluation_mu[c2_pos]
  nr <- nrow(c3d$matrix$values)
  common_mask <- matrix(FALSE, nr, ncol(c3d$matrix$values),
                        dimnames = dimnames(c3d$matrix$values))
  common_mask[common_linear] <- TRUE
  baseline <- .lt_paired_summary(mu3, mu2, "mu_C3_vs_mu_C2", ratio = TRUE)
  representations <- rbind(
    .lt_paired_summary(y - mu3, y - mu2, "ADD_LT_C3_vs_ADD_C2", ratio = FALSE),
    .lt_paired_summary(y / mu3, y / mu2, "GBI_C3_vs_GBI_C2", ratio = TRUE),
    .lt_paired_summary(log(y / mu3), log(y / mu2),
                       "logGBI_C3_vs_logGBI_C2", ratio = FALSE)
  )
  baseline_reason <- data.frame(
    reason = c2_fit$baseline_reason_levels,
    count = tabulate(c2_fit$baseline_reason,
                     nbins = length(c2_fit$baseline_reason_levels)),
    stringsAsFactors = FALSE
  )
  baseline_reason <- baseline_reason[baseline_reason$count > 0L, , drop = FALSE]
  out <- structure(list(
    c3_declared_fit_domain_count = length(c3d$edge_linear),
    c3_declared_fit_domain_sha256 = c3d$hashes$ordered_cell_ledger_sha256,
    c2_declared_fit_domain_count = length(c2d$fit_linear),
    c2_declared_fit_domain_sha256 = c2d$hashes$fit_domain_sha256,
    common_positive_cell_count = length(common_linear),
    common_positive_linear = common_linear,
    common_positive_cell_sha256 = .lt_rate_keyed_mask_hash(
      common_mask, c3d$matrix$gene_ids, c3d$matrix$branch_ids
    ),
    c2_zero_cell_fit_loss_count = sum(c2d$evaluation_mask &
                                        c2d$matrix$values == 0),
    c2_zero_cell_fit_loss_fraction = sum(c2d$evaluation_mask &
                                           c2d$matrix$values == 0) /
      sum(c2d$evaluation_mask),
    unsupported_vertex_ledger = c2d$unsupported_vertices,
    unsupported_gene_ledger = c2d$unsupported_vertices[
      c2d$unsupported_vertices$vertex_type == "gene", , drop = FALSE
    ],
    unsupported_branch_ledger = c2d$unsupported_vertices[
      c2d$unsupported_vertices$vertex_type == "branch", , drop = FALSE
    ],
    component_ledger = c2d$components$summary,
    baseline_availability_reason_counts = baseline_reason,
    baseline_paired_summary = baseline,
    representation_paired_summary = representations,
    comparison_scope = "EXACT_KEYED_COMMON_POSITIVE_DOMAIN",
    full_domain_label = "UNEQUAL_DOMAIN_DESCRIPTIVE",
    sensitivity_result = "DIVERGENCE_DESCRIBED_NO_WINNER_SELECTED",
    scientific_pass = FALSE
  ), class = "lt_c3_c2_comparison")
  out$output_sha256 <- .lt_hash(out)
  out
}
