# Matrix-QC specification, diagnostics, and optional Layer-2 admission.

.lt_qc_severities <- function() c("INFO", "WARNING", "HOLD")
.lt_qc_methods <- function() c("none", "quantile", "iqr", "mad", "absolute", "custom_mask")
.lt_qc_scopes <- function() c("whole_matrix", "within_gene", "within_branch")
.lt_qc_directions <- function() c("upper", "lower", "two_sided")
.lt_qc_quantile_tie_semantics <- function() c("strict", "inclusive", "exact_rank")

.lt_qc_normalize_rule <- function(rule, index) {
  if (!is.list(rule) || is.null(names(rule))) {
    .lt_abort("Every QC outlier rule must be a named list.")
  }
  method <- rule$method %||% "none"
  scope <- rule$scope %||% "whole_matrix"
  direction <- rule$direction %||% "two_sided"
  severity <- rule$severity %||% "WARNING"
  reason <- rule$reason %||% paste0("generic_", method, "_diagnostic")
  rule_id <- rule$rule_id %||% sprintf("rule_%03d_%s", index, method)
  tie_semantics <- rule$tie_semantics %||%
    if (identical(method, "quantile")) "strict" else "not_applicable"
  for (z in list(method, scope, direction, severity, reason, rule_id, tie_semantics)) {
    if (!is.character(z) || length(z) != 1L || is.na(z) || !nzchar(z)) {
      .lt_abort("QC rule identities and labels must be non-empty character scalars.")
    }
  }
  if (!method %in% .lt_qc_methods()) .lt_abort("Unknown QC outlier-rule method.")
  if (!scope %in% .lt_qc_scopes()) .lt_abort("Unknown QC outlier-rule scope.")
  if (!direction %in% .lt_qc_directions()) .lt_abort("Unknown QC outlier-rule direction.")
  if (method == "quantile") {
    if (!tie_semantics %in% .lt_qc_quantile_tie_semantics()) {
      .lt_abort("Quantile `tie_semantics` must be strict, inclusive, or exact_rank.")
    }
    if (tie_semantics == "exact_rank") {
      .lt_abort("Quantile `exact_rank` tie semantics are parked and not implemented in QC1-FIX001.")
    }
  } else if (tie_semantics != "not_applicable") {
    .lt_abort("`tie_semantics` applies only to quantile QC rules.")
  }
  if (!severity %in% c("INFO", "WARNING")) {
    .lt_abort("Outlier rules may emit INFO or WARNING, never HOLD.")
  }
  num <- function(name, default = NA_real_) {
    value <- rule[[name]] %||% default
    if (!is.numeric(value) || length(value) != 1L || is.nan(value)) {
      .lt_abort(sprintf("QC rule `%s` must be a numeric scalar.", name))
    }
    as.double(value)
  }
  multiplier <- num("multiplier", if (method %in% c("iqr", "mad")) 1.5 else NA_real_)
  lower_probability <- num("lower_probability", if (method == "quantile") 0.01 else NA_real_)
  upper_probability <- num("upper_probability", if (method == "quantile") 0.99 else NA_real_)
  lower <- num("lower")
  upper <- num("upper")
  if (method %in% c("iqr", "mad") && (!is.finite(multiplier) || multiplier < 0)) {
    .lt_abort("IQR/MAD multipliers must be finite and nonnegative.")
  }
  if (method == "quantile" &&
      (!is.finite(lower_probability) || !is.finite(upper_probability) ||
       lower_probability < 0 || upper_probability > 1 ||
       lower_probability > upper_probability)) {
    .lt_abort("Quantile probabilities must satisfy 0 <= lower <= upper <= 1.")
  }
  if (method == "absolute") {
    if (direction %in% c("lower", "two_sided") && !is.finite(lower)) {
      .lt_abort("Absolute lower/two-sided rules require a finite `lower` threshold.")
    }
    if (direction %in% c("upper", "two_sided") && !is.finite(upper)) {
      .lt_abort("Absolute upper/two-sided rules require a finite `upper` threshold.")
    }
    if (direction == "two_sided" && lower > upper) {
      .lt_abort("Absolute lower threshold must not exceed upper threshold.")
    }
  }
  mask <- rule$mask %||% NULL
  if (method == "custom_mask" && (is.null(mask) || !is.matrix(mask) || !is.logical(mask))) {
    .lt_abort("Custom-mask QC rules require a logical matrix `mask`.")
  }
  if (method != "custom_mask" && !is.null(mask)) {
    .lt_abort("Only a custom-mask QC rule may carry `mask`.")
  }
  list(
    rule_id = rule_id, method = method, scope = scope, direction = direction,
    tie_semantics = tie_semantics,
    multiplier = multiplier, lower_probability = lower_probability,
    upper_probability = upper_probability, lower = lower, upper = upper,
    severity = severity, reason = reason, mask = mask
  )
}

#' Specify a generic matrix-QC diagnostic procedure
#'
#' The standard profile is threshold-free and diagnoses continuous structure
#' without creating binary outlier flags. Optional rules are explicit procedure
#' parameters, not universal scientific truths.
#'
#' @param profile Named diagnostic profile; currently `standard` or `custom`.
#' @param value_domain Declared/automatic value domain.
#' @param outlier_rules Optional list of generic diagnostic rules. An empty list
#'   disables binary outlier flagging.
#' @param top_k Positive integer leverage cutoffs.
#' @param top_fractions Fractions of units used for leverage summaries.
#' @param tail_fractions Fractions used for lower/upper tail summaries.
#' @param point_mass_top_n Maximum exact repeated-value concentrations retained.
#' @return An `lt_qc_spec` object.
#' @export
lt_qc_spec <- function(profile = c("standard", "custom"),
                       value_domain = c("auto", "signed", "nonnegative", "positive"),
                       outlier_rules = NULL,
                       top_k = c(1L, 5L, 10L),
                       top_fractions = c(0.01, 0.05),
                       tail_fractions = c(0.01, 0.05),
                       point_mass_top_n = 10L) {
  profile <- match.arg(profile)
  value_domain <- match.arg(value_domain)
  outlier_rules_supplied <- !is.null(outlier_rules)
  if (is.null(outlier_rules)) {
    outlier_rules <- list()
  }
  if (!is.list(outlier_rules)) .lt_abort("`outlier_rules` must be a list.")
  if (profile == "standard" && outlier_rules_supplied && length(outlier_rules)) {
    .lt_abort("The `standard` QC profile is threshold-free; use `profile = \"custom\"` for explicit optional rules.")
  }
  rules <- lapply(seq_along(outlier_rules), function(i) {
    .lt_qc_normalize_rule(outlier_rules[[i]], i)
  })
  rule_ids <- vapply(rules, `[[`, character(1), "rule_id")
  if (anyDuplicated(rule_ids)) .lt_abort("QC rule IDs must be unique.")
  if (!is.numeric(top_k) || anyNA(top_k) || any(top_k < 1) ||
      any(top_k != as.integer(top_k))) .lt_abort("`top_k` must contain positive integers.")
  validate_fraction <- function(x, name) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x <= 0 | x > 1)) {
      .lt_abort(sprintf("`%s` must contain finite fractions in (0, 1].", name))
    }
  }
  validate_fraction(top_fractions, "top_fractions")
  validate_fraction(tail_fractions, "tail_fractions")
  if (!is.numeric(point_mass_top_n) || length(point_mass_top_n) != 1L ||
      is.na(point_mass_top_n) || point_mass_top_n < 1 ||
      point_mass_top_n != as.integer(point_mass_top_n)) {
    .lt_abort("`point_mass_top_n` must be one positive integer.")
  }
  core <- list(
    profile = profile, behavior = "DIAGNOSE_ONLY", value_domain = value_domain,
    outlier_rules = rules, top_k = sort(unique(as.integer(top_k))),
    top_fractions = sort(unique(as.double(top_fractions))),
    tail_fractions = sort(unique(as.double(tail_fractions))),
    point_mass_top_n = as.integer(point_mass_top_n),
    zero_policy = "exact_zero_valid_no_pseudocount",
    flag_contract = "flagged_not_ineligible",
    threshold_authority = "explicit_optional_rule_not_universal_truth",
    standard_profile_contract = "threshold_free_diagnose_only",
    semantic_contract = c(
      "outlier_not_error", "long_branch_not_anomalous_cell", "trim_not_diagnosis"
    ),
    future_diagnostic_placeholders = list(
      two_axis_empirical_extremeness = list(
        status = "PARKED_NOT_IMPLEMENTED",
        fields = c("within_gene_percentile", "within_branch_percentile"),
        admission_effect = "none"
      )
    )
  )
  structure(c(core, list(diagnostic_spec_id = .lt_hash(core))), class = "lt_qc_spec")
}

#' Validate a matrix-QC specification
#'
#' @param x An `lt_qc_spec`.
#' @return `x`, invisibly.
#' @export
validate_lt_qc_spec <- function(x) {
  if (!inherits(x, "lt_qc_spec")) .lt_abort("`x` must inherit from lt_qc_spec.")
  rebuilt <- lt_qc_spec(
    profile = x$profile, value_domain = x$value_domain,
    outlier_rules = x$outlier_rules, top_k = x$top_k,
    top_fractions = x$top_fractions, tail_fractions = x$tail_fractions,
    point_mass_top_n = x$point_mass_top_n
  )
  if (!identical(x$diagnostic_spec_id, rebuilt$diagnostic_spec_id)) {
    .lt_abort("QC specification hash does not match its normalized fields.")
  }
  invisible(x)
}

.lt_qc_matrix_identity <- function(x) {
  list(
    values_sha256 = .lt_hash(x$values),
    coord_state_sha256 = .lt_hash(x$coord_state),
    value_reason_sha256 = .lt_hash(x$value_reason),
    ordered_gene_ids_sha256 = .lt_hash(x$gene_ids),
    ordered_branch_ids_sha256 = .lt_hash(x$branch_ids),
    combined_sha256 = .lt_hash(list(
      x$values, x$coord_state, x$value_reason, x$gene_ids, x$branch_ids
    ))
  )
}

.lt_qc_value_summary <- function(x) {
  x <- as.double(x)
  n <- length(x)
  if (!n) {
    return(c(
      available_count = 0, exact_zero_count = 0, positive_count = 0,
      negative_count = 0, min = NA, q01 = NA, q05 = NA, q25 = NA,
      median = NA, q75 = NA, q95 = NA, q99 = NA, max = NA, mean = NA,
      sd = NA, mad_raw = NA, sum = 0, sum_abs = 0, max_abs = NA,
      positive_log_count = 0, positive_log_domain_loss = 0,
      positive_log_min = NA, positive_log_median = NA, positive_log_max = NA
    ))
  }
  q <- as.double(stats::quantile(
    x, probs = c(.01, .05, .25, .5, .75, .95, .99),
    names = FALSE, type = 7
  ))
  center <- q[[4L]]
  mad_raw <- stats::median(abs(x - center))
  positive <- x > 0
  lx <- log(x[positive])
  lq <- if (length(lx)) as.double(stats::quantile(
    lx, probs = c(0, .5, 1), names = FALSE, type = 7
  )) else rep(NA_real_, 3L)
  c(
    available_count = n, exact_zero_count = sum(x == 0),
    positive_count = sum(positive), negative_count = sum(x < 0),
    min = min(x), q01 = q[[1L]], q05 = q[[2L]], q25 = q[[3L]],
    median = q[[4L]], q75 = q[[5L]], q95 = q[[6L]], q99 = q[[7L]],
    max = max(x), mean = mean(x), sd = if (n > 1L) stats::sd(x) else NA_real_,
    mad_raw = mad_raw, sum = .lt_pairwise_sum(x),
    sum_abs = .lt_pairwise_sum(abs(x)), max_abs = max(abs(x)),
    positive_log_count = sum(positive), positive_log_domain_loss = sum(!positive),
    positive_log_min = lq[[1L]], positive_log_median = lq[[2L]],
    positive_log_max = lq[[3L]]
  )
}

.lt_qc_axis_summary <- function(x, margin = c("gene", "branch")) {
  margin <- match.arg(margin)
  nr <- nrow(x$values); nc <- ncol(x$values)
  by_gene <- margin == "gene"
  n <- if (by_gene) nr else nc
  ids <- if (by_gene) x$gene_ids else x$branch_ids
  total_other <- if (by_gene) nc else nr
  rows <- lapply(seq_len(n), function(i) {
    state <- if (by_gene) x$coord_state[i, ] else x$coord_state[, i]
    values <- if (by_gene) x$values[i, ] else x$values[, i]
    observed <- state == "observed"
    available <- observed & is.finite(values) & !is.na(values)
    s <- .lt_qc_value_summary(values[available])
    data.frame(
      id = ids[[i]], total_coordinates = total_other,
      coordinate_observed_count = sum(observed),
      coordinate_observed_fraction = sum(observed) / total_other,
      payload_available_count = sum(available),
      payload_available_fraction_of_observed = if (any(observed)) sum(available) / sum(observed) else NA_real_,
      payload_unavailable_observed_count = sum(observed & !available),
      t(s), stringsAsFactors = FALSE, check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out)[[1L]] <- if (by_gene) "gene_id" else "branch_id"
  rownames(out) <- NULL
  out
}

.lt_qc_point_masses <- function(values, top_n) {
  if (!length(values)) return(data.frame(
    rank = integer(), value = numeric(), count = integer(),
    fraction_available = numeric(), interpretation = character(),
    stringsAsFactors = FALSE
  ))
  runs <- rle(sort(as.double(values), method = "radix"))
  keep <- which(runs$lengths > 1L)
  if (!length(keep)) return(data.frame(
    rank = integer(), value = numeric(), count = integer(),
    fraction_available = numeric(), interpretation = character(),
    stringsAsFactors = FALSE
  ))
  ord <- order(-runs$lengths[keep], runs$values[keep], method = "radix")
  keep <- keep[utils::head(ord, top_n)]
  data.frame(
    rank = seq_along(keep), value = runs$values[keep],
    count = as.integer(runs$lengths[keep]),
    fraction_available = runs$lengths[keep] / length(values),
    interpretation = "exact_numeric_repetition_no_upstream_cause_inferred",
    stringsAsFactors = FALSE
  )
}

.lt_qc_extreme_repetition <- function(values) {
  values <- as.double(values)
  interpretation <- "exact_extreme_repetition_no_upstream_cause_inferred"
  if (!length(values)) return(data.frame(
    minimum_finite_value = NA_real_, minimum_exact_count = 0L,
    minimum_exact_fraction_available = NA_real_,
    maximum_finite_value = NA_real_, maximum_exact_count = 0L,
    maximum_exact_fraction_available = NA_real_,
    interpretation = interpretation, stringsAsFactors = FALSE
  ))
  minimum <- min(values)
  maximum <- max(values)
  minimum_count <- sum(values == minimum)
  maximum_count <- sum(values == maximum)
  data.frame(
    minimum_finite_value = minimum,
    minimum_exact_count = as.integer(minimum_count),
    minimum_exact_fraction_available = minimum_count / length(values),
    maximum_finite_value = maximum,
    maximum_exact_count = as.integer(maximum_count),
    maximum_exact_fraction_available = maximum_count / length(values),
    interpretation = interpretation, stringsAsFactors = FALSE
  )
}

.lt_qc_leverage <- function(values, gene_summary, branch_summary, spec, resolved_domain) {
  contribution <- if (resolved_domain == "signed") abs(values) else values
  contribution_name <- if (resolved_domain == "signed") "absolute_value" else "nonnegative_value"
  make <- function(target, x, ids = NULL, measure = contribution_name) {
    x <- as.double(x)
    total <- .lt_pairwise_sum(x)
    ord <- order(-x, if (is.null(ids)) seq_along(x) else ids, method = "radix")
    rows <- list(); k <- 0L
    for (top in spec$top_k) {
      n <- min(length(x), top); k <- k + 1L
      rows[[k]] <- data.frame(
        target = target, measure = measure, selection = "top_k",
        parameter = as.double(top), selected_count = n,
        contribution_fraction = if (total > 0) .lt_pairwise_sum(x[ord[seq_len(n)]]) / total else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    for (fraction in spec$top_fractions) {
      n <- min(length(x), max(1L, ceiling(length(x) * fraction))); k <- k + 1L
      rows[[k]] <- data.frame(
        target = target, measure = measure, selection = "top_fraction_of_units",
        parameter = fraction, selected_count = n,
        contribution_fraction = if (total > 0) .lt_pairwise_sum(x[ord[seq_len(n)]]) / total else NA_real_,
        stringsAsFactors = FALSE
      )
    }
    do.call(rbind, rows)
  }
  gene_measure <- if (resolved_domain == "signed") abs(gene_summary$sum) else gene_summary$sum
  branch_measure <- if (resolved_domain == "signed") abs(branch_summary$sum) else branch_summary$sum
  rbind(
    make("cell", contribution, measure = contribution_name),
    make("gene_margin", gene_measure, gene_summary$gene_id,
         if (resolved_domain == "signed") "absolute_gene_margin" else "nonnegative_gene_margin"),
    make("branch_margin", branch_measure, branch_summary$branch_id,
         if (resolved_domain == "signed") "absolute_branch_margin" else "nonnegative_branch_margin")
  )
}

.lt_qc_tail_summary <- function(values, spec, resolved_domain) {
  if (!length(values)) return(data.frame())
  mass <- if (resolved_domain == "signed") abs(values) else values
  mass_name <- if (resolved_domain == "signed") "absolute_value" else "nonnegative_value"
  total <- .lt_pairwise_sum(mass)
  rows <- list(); k <- 0L
  for (fraction in spec$tail_fractions) {
    thresholds <- as.double(stats::quantile(values, c(fraction, 1 - fraction),
                                             names = FALSE, type = 7))
    for (side in c("lower", "upper")) {
      ii <- if (side == "lower") values <= thresholds[[1L]] else values >= thresholds[[2L]]
      k <- k + 1L
      rows[[k]] <- data.frame(
        side = side, tail_fraction_parameter = fraction,
        threshold = thresholds[[if (side == "lower") 1L else 2L]],
        cell_count = sum(ii), cell_fraction = mean(ii), mass = mass_name,
        contribution_fraction = if (total > 0) .lt_pairwise_sum(mass[ii]) / total else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

.lt_qc_rule_groups <- function(x, scope) {
  if (scope == "whole_matrix") return(list(whole_matrix = seq_along(x$values)))
  if (scope == "within_gene") {
    return(stats::setNames(lapply(seq_len(nrow(x$values)), function(i) {
      i + (seq_len(ncol(x$values)) - 1L) * nrow(x$values)
    }), x$gene_ids))
  }
  stats::setNames(lapply(seq_len(ncol(x$values)), function(j) {
    seq_len(nrow(x$values)) + (j - 1L) * nrow(x$values)
  }), x$branch_ids)
}

.lt_qc_rule_flags <- function(x, rule, spec_id) {
  observed_available <- x$coord_state == "observed" & is.finite(x$values) & !is.na(x$values)
  if (rule$method == "none") return(list(flags = data.frame(), warnings = data.frame()))
  if (rule$method == "custom_mask") {
    mask <- rule$mask
    if (!identical(dim(mask), dim(x$values)) ||
        !identical(dimnames(mask), dimnames(x$values)) || anyNA(mask)) {
      .lt_abort("Custom QC mask must exactly match matrix dimensions and axes without NA.")
    }
    if (any(mask & !observed_available)) {
      .lt_abort("Custom QC mask may flag only observed coordinates with finite payloads.")
    }
    groups <- list(custom = which(mask))
  } else {
    groups <- .lt_qc_rule_groups(x, rule$scope)
  }
  flags <- list(); warnings <- list(); nf <- 0L; nw <- 0L
  for (group_id in names(groups)) {
    idx <- groups[[group_id]]
    idx <- idx[observed_available[idx]]
    if (!length(idx)) next
    v <- x$values[idx]
    lower <- upper <- NA_real_
    diagnostic <- rule$method
    robust <- rep(NA_real_, length(v))
    if (rule$method == "custom_mask") {
      selected <- rep(TRUE, length(idx))
      diagnostic <- "user_supplied_flag_mask"
    } else if (rule$method == "quantile") {
      threshold <- as.double(stats::quantile(
        v, c(rule$lower_probability, rule$upper_probability),
        names = FALSE, type = 7
      ))
      lower <- threshold[[1L]]; upper <- threshold[[2L]]
      selected <- if (rule$tie_semantics == "strict") {
        switch(rule$direction, upper = v > upper,
               lower = v < lower, two_sided = v < lower | v > upper)
      } else {
        switch(rule$direction, upper = v >= upper,
               lower = v <= lower, two_sided = v <= lower | v >= upper)
      }
    } else if (rule$method == "iqr") {
      q <- as.double(stats::quantile(v, c(.25, .5, .75), names = FALSE, type = 7))
      spread <- q[[3L]] - q[[1L]]
      lower <- q[[1L]] - rule$multiplier * spread
      upper <- q[[3L]] + rule$multiplier * spread
      robust <- if (spread > 0) (v - q[[2L]]) / spread else rep(NA_real_, length(v))
      selected <- switch(rule$direction, upper = v > upper,
                         lower = v < lower, two_sided = v < lower | v > upper)
    } else if (rule$method == "mad") {
      center <- stats::median(v)
      spread <- stats::median(abs(v - center))
      if (spread == 0) {
        selected <- rep(FALSE, length(v))
        nw <- nw + 1L
        warnings[[nw]] <- data.frame(
          severity = "WARNING", warning_code = "MAD_ZERO_UNDEFINED",
          scope = rule$scope, scope_id = group_id,
          detail = "MAD is exactly zero; robust scores and MAD flags are undefined, not Inf/NaN.",
          source_spec = spec_id, rule_id = rule$rule_id,
          stringsAsFactors = FALSE
        )
      } else {
        robust <- (v - center) / spread
        lower <- center - rule$multiplier * spread
        upper <- center + rule$multiplier * spread
        selected <- switch(rule$direction, upper = robust > rule$multiplier,
                           lower = robust < -rule$multiplier,
                           two_sided = abs(robust) > rule$multiplier)
      }
    } else if (rule$method == "absolute") {
      lower <- rule$lower; upper <- rule$upper
      selected <- switch(rule$direction, upper = v > upper,
                         lower = v < lower, two_sided = v < lower | v > upper)
    }
    if (!any(selected)) next
    hit <- idx[selected]
    gi <- ((hit - 1L) %% nrow(x$values)) + 1L
    bi <- ((hit - 1L) %/% nrow(x$values)) + 1L
    threshold_text <- sprintf(
      "method=%s;scope=%s;direction=%s;tie_semantics=%s;lower=%s;upper=%s;multiplier=%s;prob_lower=%s;prob_upper=%s",
      rule$method, rule$scope, rule$direction, rule$tie_semantics,
      .lt_hex(lower), .lt_hex(upper),
      .lt_hex(rule$multiplier), .lt_hex(rule$lower_probability),
      .lt_hex(rule$upper_probability)
    )
    nf <- nf + 1L
    flags[[nf]] <- data.frame(
      gene_id = x$gene_ids[gi], branch_id = x$branch_ids[bi],
      gene_index = gi, branch_index = bi,
      diagnostic_name = diagnostic, diagnostic_value = v[selected],
      robust_score = robust[selected], threshold_rule = threshold_text,
      severity = rule$severity, reason = rule$reason,
      source_spec = spec_id, rule_id = rule$rule_id,
      stringsAsFactors = FALSE
    )
  }
  flag_df <- if (length(flags)) do.call(rbind, flags) else data.frame()
  warning_df <- if (length(warnings)) do.call(rbind, warnings) else data.frame()
  list(flags = flag_df, warnings = warning_df)
}

.lt_qc_empty_flags <- function() data.frame(
  flag_id = character(), gene_id = character(), branch_id = character(),
  gene_index = integer(), branch_index = integer(), diagnostic_name = character(),
  diagnostic_value = numeric(), robust_score = numeric(), threshold_rule = character(),
  severity = character(), reason = character(), source_spec = character(),
  rule_id = character(), stringsAsFactors = FALSE
)

.lt_qc_flag_identity <- function(flags) {
  flags[, setdiff(names(flags), c("gene_index", "branch_index")), drop = FALSE]
}

.lt_qc_empty_messages <- function(kind = c("warning", "hold")) {
  kind <- match.arg(kind)
  if (kind == "warning") data.frame(
    severity = character(), warning_code = character(), scope = character(),
    scope_id = character(), detail = character(), source_spec = character(),
    rule_id = character(), stringsAsFactors = FALSE
  ) else data.frame(
    severity = character(), hold_code = character(), scope = character(),
    detail = character(), source_spec = character(), stringsAsFactors = FALSE
  )
}

.lt_qc_threshold_config <- function(spec) {
  if (!length(spec$outlier_rules)) return(data.frame(
    rule_id = character(), method = character(), scope = character(),
    direction = character(), tie_semantics = character(),
    multiplier = numeric(), lower_probability = numeric(),
    upper_probability = numeric(), lower = numeric(), upper = numeric(),
    severity = character(), reason = character(), custom_mask_sha256 = character(),
    stringsAsFactors = FALSE
  ))
  do.call(rbind, lapply(spec$outlier_rules, function(rule) data.frame(
    rule_id = rule$rule_id, method = rule$method, scope = rule$scope,
    direction = rule$direction, tie_semantics = rule$tie_semantics,
    multiplier = rule$multiplier,
    lower_probability = rule$lower_probability,
    upper_probability = rule$upper_probability, lower = rule$lower,
    upper = rule$upper, severity = rule$severity, reason = rule$reason,
    custom_mask_sha256 = if (is.null(rule$mask)) NA_character_ else .lt_hash(rule$mask),
    stringsAsFactors = FALSE
  )))
}

.lt_qc_future_diagnostics <- function(spec) {
  list(
    two_axis_empirical_extremeness = list(
      status = "PARKED_NOT_IMPLEMENTED",
      diagnostic_only = TRUE,
      fields = data.frame(
        field = c("within_gene_percentile", "within_branch_percentile"),
        value = c(NA_real_, NA_real_),
        status = c("not_computed", "not_computed"),
        stringsAsFactors = FALSE
      ),
      threshold = NULL,
      admission_effect = "none",
      source_spec = spec$diagnostic_spec_id
    )
  )
}

.lt_qc_top_cells <- function(x, available_linear, resolved_domain, n = 10L) {
  if (!length(available_linear)) return(data.frame(
    rank = integer(), gene_id = character(), branch_id = character(),
    value = numeric(), contribution_measure = character(),
    contribution_value = numeric(), contribution_fraction = numeric(),
    upper_tail_rank_fraction = numeric(), robust_score = numeric(),
    stringsAsFactors = FALSE
  ))
  values <- x$values[available_linear]
  contribution <- if (resolved_domain == "signed") abs(values) else values
  measure <- if (resolved_domain == "signed") "absolute_value" else "nonnegative_value"
  total <- .lt_pairwise_sum(contribution)
  gi <- ((available_linear - 1L) %% nrow(x$values)) + 1L
  bi <- ((available_linear - 1L) %/% nrow(x$values)) + 1L
  order_all <- order(-contribution, x$gene_ids[gi], x$branch_ids[bi], method = "radix")
  keep <- order_all[seq_len(min(n, length(order_all)))]
  sorted_values <- sort(values, method = "radix")
  center <- stats::median(values)
  spread <- stats::median(abs(values - center))
  data.frame(
    rank = seq_along(keep), gene_id = x$gene_ids[gi[keep]],
    branch_id = x$branch_ids[bi[keep]], value = values[keep],
    contribution_measure = measure, contribution_value = contribution[keep],
    contribution_fraction = if (total > 0) contribution[keep] / total else NA_real_,
    upper_tail_rank_fraction = (length(values) - findInterval(values[keep], sorted_values,
      left.open = FALSE, rightmost.closed = TRUE) + 1L) / length(values),
    robust_score = if (spread > 0) (values[keep] - center) / spread else NA_real_,
    stringsAsFactors = FALSE
  )
}

#' Diagnose a branch-associated matrix without changing it
#'
#' @param x An `lt_matrix`.
#' @param spec Optional `lt_qc_spec`; when omitted, `profile` constructs one.
#' @param profile Named profile used only when `spec` is omitted.
#' @return An immutable-style `lt_diagnostic` object.
#' @export
lt_qc_matrix <- function(x, spec = NULL, profile = "standard") {
  validate_lt_matrix(x)
  input_identity <- .lt_qc_matrix_identity(x)
  if (is.null(spec)) spec <- lt_qc_spec(profile = profile)
  validate_lt_qc_spec(spec)
  observed <- x$coord_state == "observed"
  available <- observed & is.finite(x$values) & !is.na(x$values)
  available_linear <- which(available)
  values <- x$values[available_linear]
  resolved_domain <- if (spec$value_domain == "auto") {
    if (any(values < 0)) "signed" else "nonnegative"
  } else spec$value_domain

  holds <- .lt_qc_empty_messages("hold")
  if (spec$value_domain == "nonnegative" && any(values < 0)) {
    holds <- rbind(holds, data.frame(
      severity = "HOLD", hold_code = "DECLARED_NONNEGATIVE_DOMAIN_VIOLATION",
      scope = "whole_matrix",
      detail = "Finite negative values contradict the explicitly declared nonnegative QC domain.",
      source_spec = spec$diagnostic_spec_id, stringsAsFactors = FALSE
    ))
  }
  if (spec$value_domain == "positive" && any(values <= 0)) {
    holds <- rbind(holds, data.frame(
      severity = "HOLD", hold_code = "DECLARED_POSITIVE_DOMAIN_VIOLATION",
      scope = "whole_matrix",
      detail = "Finite zero/negative values contradict the explicitly declared positive QC domain.",
      source_spec = spec$diagnostic_spec_id, stringsAsFactors = FALSE
    ))
  }

  payload_state <- .lt_payload_state(x$values, x$coord_state)
  state_key <- paste(as.vector(x$coord_state), as.vector(payload_state), sep = "\r")
  state_counts <- as.data.frame(table(state_key), stringsAsFactors = FALSE)
  split_key <- strsplit(as.character(state_counts$state_key), "\r", fixed = TRUE)
  cell_state_counts <- data.frame(
    coord_state = vapply(split_key, `[[`, character(1), 1L),
    payload_state = vapply(split_key, `[[`, character(1), 2L),
    count = as.integer(state_counts$Freq), stringsAsFactors = FALSE
  )
  cell_state_counts <- cell_state_counts[cell_state_counts$count > 0L, , drop = FALSE]
  cell_state_counts$fraction_all_cells <- cell_state_counts$count / length(x$values)
  cell_state_counts <- cell_state_counts[order(
    match(cell_state_counts$coord_state, lt_allowed_states()),
    cell_state_counts$payload_state, method = "radix"
  ), , drop = FALSE]
  rownames(cell_state_counts) <- NULL

  gene_summary <- .lt_qc_axis_summary(x, "gene")
  branch_summary <- .lt_qc_axis_summary(x, "branch")
  overall <- .lt_qc_value_summary(values)
  total_abs <- .lt_pairwise_sum(abs(values))
  total_nonnegative <- if (resolved_domain == "signed") NA_real_ else .lt_pairwise_sum(values)
  gene_abs_margin_total <- .lt_pairwise_sum(abs(gene_summary$sum))
  branch_abs_margin_total <- .lt_pairwise_sum(abs(branch_summary$sum))
  gene_summary$cell_absolute_mass_fraction <- if (total_abs > 0) gene_summary$sum_abs / total_abs else NA_real_
  branch_summary$cell_absolute_mass_fraction <- if (total_abs > 0) branch_summary$sum_abs / total_abs else NA_real_
  gene_summary$absolute_margin_fraction <- if (gene_abs_margin_total > 0) abs(gene_summary$sum) / gene_abs_margin_total else NA_real_
  branch_summary$absolute_margin_fraction <- if (branch_abs_margin_total > 0) abs(branch_summary$sum) / branch_abs_margin_total else NA_real_
  gene_summary$concentration_definition <- if (resolved_domain == "signed") {
    "absolute_value_and_absolute_gene_margin"
  } else "nonnegative_value_and_gene_margin"
  branch_summary$concentration_definition <- if (resolved_domain == "signed") {
    "absolute_value_and_absolute_branch_margin"
  } else "nonnegative_value_and_branch_margin"

  point_masses <- .lt_qc_point_masses(values, spec$point_mass_top_n)
  extreme_exact_repetition <- .lt_qc_extreme_repetition(values)
  tail_summary <- .lt_qc_tail_summary(values, spec, resolved_domain)
  leverage <- .lt_qc_leverage(values, gene_summary, branch_summary, spec, resolved_domain)
  metric_domains <- data.frame(
    metric = c("raw_distribution", "contribution_concentration", "positive_log_distribution"),
    mathematical_domain = c(
      "finite_signed", if (resolved_domain == "signed") "absolute_finite_value" else "finite_nonnegative_value",
      "strictly_positive_finite_value"
    ),
    admitted_count = c(length(values), length(values), sum(values > 0)),
    domain_loss_count = c(0L, 0L, sum(values <= 0)),
    zero_handling = c("retained", "retained", "outside_metric_domain_no_pseudocount"),
    stringsAsFactors = FALSE
  )

  rule_results <- lapply(spec$outlier_rules, function(rule) {
    .lt_qc_rule_flags(x, rule, spec$diagnostic_spec_id)
  })
  flags_raw <- lapply(rule_results, `[[`, "flags")
  flags_raw <- flags_raw[vapply(flags_raw, nrow, integer(1)) > 0L]
  flags <- if (length(flags_raw)) do.call(rbind, flags_raw) else .lt_qc_empty_flags()[, -1L, drop = FALSE]
  if (nrow(flags)) {
    flags <- flags[order(flags$gene_id, flags$branch_id, flags$rule_id,
                         flags$diagnostic_name, method = "radix"), , drop = FALSE]
    identity_columns <- setdiff(names(flags), c("gene_index", "branch_index"))
    flags$flag_id <- vapply(seq_len(nrow(flags)), function(i) paste0(
      "flag_", substr(.lt_hash(as.list(flags[i, identity_columns, drop = FALSE])), 1L, 20L)
    ), character(1))
    flags <- flags[, names(.lt_qc_empty_flags()), drop = FALSE]
    rownames(flags) <- NULL
  } else flags <- .lt_qc_empty_flags()

  warnings_raw <- lapply(rule_results, `[[`, "warnings")
  warnings_raw <- warnings_raw[vapply(warnings_raw, nrow, integer(1)) > 0L]
  warnings <- if (length(warnings_raw)) do.call(rbind, warnings_raw) else .lt_qc_empty_messages("warning")
  flag_counts <- if (nrow(flags)) aggregate(
    flag_id ~ rule_id + severity, flags, length
  ) else data.frame(rule_id = character(), severity = character(), flag_id = integer())
  if (nrow(flag_counts)) {
    for (i in seq_len(nrow(flag_counts))) warnings <- rbind(warnings, data.frame(
      severity = flag_counts$severity[[i]], warning_code = "DIAGNOSTIC_FLAGS_PRESENT",
      scope = "rule", scope_id = flag_counts$rule_id[[i]],
      detail = sprintf("%d diagnostic flags; no Layer-2 exclusion was applied.", flag_counts$flag_id[[i]]),
      source_spec = spec$diagnostic_spec_id, rule_id = flag_counts$rule_id[[i]],
      stringsAsFactors = FALSE
    ))
  }
  if (nrow(point_masses)) warnings <- rbind(warnings, data.frame(
    severity = "INFO", warning_code = "EXACT_POINT_MASS_PRESENT",
    scope = "whole_matrix", scope_id = "whole_matrix",
    detail = sprintf("%d repeated-value concentrations retained; no upstream cause inferred.", nrow(point_masses)),
    source_spec = spec$diagnostic_spec_id, rule_id = NA_character_,
    stringsAsFactors = FALSE
  ))
  if (any(values <= 0)) warnings <- rbind(warnings, data.frame(
    severity = "INFO", warning_code = "POSITIVE_ONLY_METRIC_DOMAIN_LOSS",
    scope = "whole_matrix", scope_id = "whole_matrix",
    detail = sprintf("%d zero/negative finite cells are outside positive-log metrics; values were unchanged.", sum(values <= 0)),
    source_spec = spec$diagnostic_spec_id, rule_id = NA_character_,
    stringsAsFactors = FALSE
  ))
  warnings <- warnings[order(match(warnings$severity, .lt_qc_severities()),
                             warnings$warning_code, warnings$scope_id,
                             method = "radix"), , drop = FALSE]
  rownames(warnings) <- NULL

  matrix_overall <- data.frame(
    gene_count = nrow(x$values), branch_count = ncol(x$values),
    cell_count = length(x$values), coordinate_observed_count = sum(observed),
    coordinate_observed_fraction = mean(observed),
    coordinate_unavailable_count = sum(!observed),
    payload_available_count = length(values),
    payload_available_fraction_of_observed = if (any(observed)) length(values) / sum(observed) else NA_real_,
    observed_payload_unavailable_count = sum(observed & !available),
    exact_zero_count = sum(values == 0),
    exact_zero_fraction_of_available = if (length(values)) mean(values == 0) else NA_real_,
    positive_count = sum(values > 0), negative_count = sum(values < 0),
    requested_value_domain = spec$value_domain, resolved_value_domain = resolved_domain,
    concentration_measure = if (resolved_domain == "signed") "absolute_value" else "nonnegative_value",
    total_absolute_value = total_abs, total_nonnegative_value = total_nonnegative,
    stringsAsFactors = FALSE
  )
  distribution <- as.data.frame(as.list(overall), stringsAsFactors = FALSE)
  cell_summary <- list(
    state_counts = cell_state_counts,
    top_leverage_cells = .lt_qc_top_cells(
      x, available_linear, resolved_domain, max(spec$top_k)
    ),
    metric_domains = metric_domains,
    future_diagnostics = .lt_qc_future_diagnostics(spec)
  )
  matrix_summary <- list(
    overall = matrix_overall, value_distribution = distribution,
    tail_concentration = tail_summary, leverage = leverage,
    point_masses = point_masses,
    extreme_exact_repetition = extreme_exact_repetition,
    metric_domains = metric_domains
  )
  threshold_config <- .lt_qc_threshold_config(spec)
  matrix_identity <- input_identity
  payload_identity <- list(payload = x$payload, sha256 = .lt_hash(x$payload))
  output_hashes <- list(
    cell_summary_sha256 = .lt_hash(cell_summary),
    gene_summary_sha256 = .lt_hash(gene_summary),
    branch_summary_sha256 = .lt_hash(branch_summary),
    matrix_summary_sha256 = .lt_hash(matrix_summary),
    flag_ledger_sha256 = .lt_hash(.lt_qc_flag_identity(flags)),
    warnings_sha256 = .lt_hash(warnings),
    holds_sha256 = .lt_hash(holds), threshold_config_sha256 = .lt_hash(threshold_config)
  )
  output_hashes$diagnostic_outputs_sha256 <- .lt_hash(output_hashes)
  diagnostic_provenance <- list(
    engine = "LaTerra_generic_matrix_qc_v1",
    package_version = as.character(utils::packageVersion("LaTerra")),
    behavior = "DIAGNOSE_ONLY", trait_independent = TRUE,
    standard_profile_threshold_free = identical(spec$profile, "standard") &&
      !length(spec$outlier_rules),
    semantic_contract = spec$semantic_contract,
    numeric_trim_performed = FALSE,
    raw_matrix_modified = FALSE, eligibility_modified = FALSE,
    deterministic_ordering = "identity_keyed_flags_and_axis_ledgers",
    matrix_identity_sha256 = matrix_identity$combined_sha256,
    payload_identity_sha256 = payload_identity$sha256
  )
  out <- structure(list(
    diagnostic_spec_id = spec$diagnostic_spec_id,
    matrix_identity = matrix_identity, payload_identity = payload_identity,
    cell_summary = cell_summary, gene_summary = gene_summary,
    branch_summary = branch_summary, matrix_summary = matrix_summary,
    flags = flags, warnings = warnings, holds = holds,
    threshold_config = threshold_config,
    diagnostic_provenance = diagnostic_provenance,
    output_hashes = output_hashes, spec = spec
  ), class = "lt_diagnostic")
  validate_lt_diagnostic(out)
  if (!identical(input_identity, .lt_qc_matrix_identity(x))) {
    .lt_abort("QC diagnosis mutated raw matrix identity.")
  }
  out
}

#' Validate a matrix-QC diagnostic object
#'
#' @param x An `lt_diagnostic`.
#' @return `x`, invisibly.
#' @export
validate_lt_diagnostic <- function(x) {
  if (!inherits(x, "lt_diagnostic")) .lt_abort("`x` must inherit from lt_diagnostic.")
  validate_lt_qc_spec(x$spec)
  if (!identical(x$diagnostic_spec_id, x$spec$diagnostic_spec_id)) {
    .lt_abort("Diagnostic/spec identities differ.")
  }
  required_hashes <- c(
    "cell_summary_sha256", "gene_summary_sha256", "branch_summary_sha256",
    "matrix_summary_sha256", "flag_ledger_sha256", "warnings_sha256",
    "holds_sha256", "threshold_config_sha256", "diagnostic_outputs_sha256"
  )
  if (!identical(names(x$output_hashes), required_hashes)) {
    .lt_abort("Diagnostic output-hash schema mismatch.")
  }
  observed_hashes <- list(
    cell_summary_sha256 = .lt_hash(x$cell_summary),
    gene_summary_sha256 = .lt_hash(x$gene_summary),
    branch_summary_sha256 = .lt_hash(x$branch_summary),
    matrix_summary_sha256 = .lt_hash(x$matrix_summary),
    flag_ledger_sha256 = .lt_hash(.lt_qc_flag_identity(x$flags)),
    warnings_sha256 = .lt_hash(x$warnings),
    holds_sha256 = .lt_hash(x$holds),
    threshold_config_sha256 = .lt_hash(x$threshold_config)
  )
  observed_hashes$diagnostic_outputs_sha256 <- .lt_hash(observed_hashes)
  if (!identical(x$output_hashes, observed_hashes)) {
    .lt_abort("Diagnostic summary/ledger hashes do not match contents.")
  }
  if (nrow(x$flags)) {
    required_flags <- names(.lt_qc_empty_flags())
    if (!identical(names(x$flags), required_flags) ||
        any(!x$flags$severity %in% c("INFO", "WARNING")) ||
        anyDuplicated(x$flags$flag_id)) {
      .lt_abort("Diagnostic flag ledger violates its frozen schema.")
    }
  }
  if (nrow(x$holds) && any(x$holds$severity != "HOLD")) {
    .lt_abort("Only structurally/mathematically established issues may enter the HOLD ledger.")
  }
  invisible(x)
}

.lt_qc_decision_mask <- function(x, diagnostic, exclude) {
  if (!is.list(exclude) || is.null(names(exclude)) || any(!nzchar(names(exclude)))) {
    .lt_abort("`exclude` must be an explicit named decision list.")
  }
  allowed <- c("flag_ids", "rule_ids", "mask")
  unknown <- setdiff(names(exclude), allowed)
  if (length(unknown)) .lt_abort("Unknown QC exclusion decision field.")
  if (!length(intersect(names(exclude), allowed))) {
    .lt_abort("QC application requires an explicit flag/rule/mask decision.")
  }
  selected <- matrix(FALSE, nrow(x$values), ncol(x$values), dimnames = dimnames(x$values))
  chosen_flags <- character()
  if ("flag_ids" %in% names(exclude)) {
    ids <- exclude$flag_ids
    if (!is.character(ids) || anyNA(ids)) .lt_abort("`flag_ids` must be character.")
    unknown_ids <- setdiff(ids, diagnostic$flags$flag_id)
    if (length(unknown_ids)) .lt_abort("Unknown diagnostic flag ID in QC application.")
    chosen_flags <- union(chosen_flags, ids)
  }
  if ("rule_ids" %in% names(exclude)) {
    ids <- exclude$rule_ids
    if (!is.character(ids) || anyNA(ids)) .lt_abort("`rule_ids` must be character.")
    known <- diagnostic$threshold_config$rule_id
    if (length(setdiff(ids, known))) .lt_abort("Unknown diagnostic rule ID in QC application.")
    chosen_flags <- union(chosen_flags, diagnostic$flags$flag_id[
      diagnostic$flags$rule_id %in% ids
    ])
  }
  if (length(chosen_flags)) {
    z <- diagnostic$flags[match(chosen_flags, diagnostic$flags$flag_id), , drop = FALSE]
    selected[cbind(z$gene_index, z$branch_index)] <- TRUE
  }
  if ("mask" %in% names(exclude)) {
    mask <- exclude$mask
    if (!is.matrix(mask) || !is.logical(mask) || anyNA(mask) ||
        !identical(dim(mask), dim(x$values)) ||
        !identical(dimnames(mask), dimnames(x$values))) {
      .lt_abort("Explicit QC exclusion mask must exactly match matrix axes without NA.")
    }
    selected <- selected | mask
  }
  if (any(selected & x$coord_state != "observed")) {
    .lt_abort("QC application cannot select Layer-2 `not_applicable` / non-observed coordinates.")
  }
  list(mask = selected, chosen_flag_ids = sort(unique(chosen_flags)),
       decision = exclude)
}

#' Explicitly apply selected QC decisions to Layer 2
#'
#' @param x The unchanged source `lt_matrix`.
#' @param diagnostic An aligned `lt_diagnostic`.
#' @param exclude Explicit named list containing any of `flag_ids`, `rule_ids`,
#'   or an aligned logical `mask`.
#' @param measurement_fit_eligibility Optional existing Layer-2 declaration.
#' @param decision_id Stable decision identifier.
#' @param authority Human/machine authority for this explicit application.
#' @return An `lt_measurement_eligibility`; matrix values are never rewritten.
#' @export
lt_apply_qc <- function(x, diagnostic, exclude,
                        measurement_fit_eligibility = NULL,
                        decision_id, authority) {
  validate_lt_matrix(x)
  validate_lt_diagnostic(diagnostic)
  .lt_assert_scalar_character(decision_id, "decision_id")
  .lt_assert_scalar_character(authority, "authority")
  before <- .lt_qc_matrix_identity(x)
  if (!identical(before, diagnostic$matrix_identity)) {
    .lt_abort("QC diagnostic does not belong to this exact raw matrix identity.")
  }
  if (nrow(diagnostic$holds)) {
    .lt_abort("QC application is held by a frozen structural/mathematical diagnostic failure.")
  }
  base <- measurement_fit_eligibility %||% lt_measurement_eligibility(x)
  validate_lt_measurement_eligibility(base, x)
  selected <- .lt_qc_decision_mask(x, diagnostic, exclude)
  if (any(base$status[selected$mask] == "not_applicable")) {
    .lt_abort("QC application cannot select a Layer-2 `not_applicable` cell.")
  }
  eligible_selected <- selected$mask & base$status == "eligible"
  already_ineligible_selected <- selected$mask & base$status == "ineligible"
  unresolved_selected <- selected$mask & base$status == "unresolved"
  if (sum(selected$mask) != sum(eligible_selected) +
      sum(already_ineligible_selected) + sum(unresolved_selected)) {
    .lt_abort("Selected Layer-2 cells do not form the frozen eligible/ineligible/unresolved transition partition.")
  }
  status <- base$status
  reason <- base$reason
  newly_ineligible <- eligible_selected | unresolved_selected
  status[newly_ineligible] <- "ineligible"
  reason[newly_ineligible] <- "measurement_qc_excluded"
  selected_linear <- which(selected$mask)
  selected_gene <- ((selected_linear - 1L) %% nrow(x$values)) + 1L
  selected_branch <- ((selected_linear - 1L) %/% nrow(x$values)) + 1L
  selected_cell_keys <- paste(
    x$gene_ids[selected_gene], x$branch_ids[selected_branch], sep = "\r"
  )
  prior_state_ledger <- data.frame(
    cell_key = selected_cell_keys,
    gene_id = x$gene_ids[selected_gene],
    branch_id = x$branch_ids[selected_branch],
    prior_status = base$status[selected_linear],
    prior_reason = base$reason[selected_linear],
    resulting_status = status[selected_linear],
    resulting_reason = reason[selected_linear],
    transition = ifelse(
      base$status[selected_linear] == "eligible", "eligible_to_ineligible",
      ifelse(base$status[selected_linear] == "unresolved",
             "unresolved_to_ineligible", "ineligible_preserved")
    ),
    stringsAsFactors = FALSE
  )
  prior_state_ledger <- prior_state_ledger[
    order(prior_state_ledger$cell_key, method = "radix"), , drop = FALSE
  ]
  rownames(prior_state_ledger) <- NULL
  counts <- list(
    selected_cell_count = as.integer(sum(selected$mask)),
    newly_excluded_count = as.integer(sum(eligible_selected)),
    already_ineligible_selected_count = as.integer(sum(already_ineligible_selected)),
    unresolved_resolved_to_ineligible_count = as.integer(sum(unresolved_selected))
  )
  decision_payload <- list(
    decision_id = decision_id, authority = authority,
    selected_flag_ids = selected$chosen_flag_ids,
    selected_cell_keys = sort(selected_cell_keys, method = "radix"),
    transition_counts = counts,
    prior_status_sha256 = base$status_sha256,
    prior_reason_sha256 = base$reason_sha256,
    prior_state_ledger_sha256 = .lt_hash(prior_state_ledger),
    diagnostic_spec_id = diagnostic$diagnostic_spec_id,
    diagnostic_outputs_sha256 = diagnostic$output_hashes$diagnostic_outputs_sha256,
    flag_ledger_sha256 = diagnostic$output_hashes$flag_ledger_sha256,
    threshold_config_sha256 = diagnostic$output_hashes$threshold_config_sha256
  )
  application <- list(
    decision_id = decision_id, authority = authority,
    explicit_application = TRUE,
    selected_cell_count = counts$selected_cell_count,
    newly_excluded_count = counts$newly_excluded_count,
    already_ineligible_selected_count = counts$already_ineligible_selected_count,
    unresolved_resolved_to_ineligible_count = counts$unresolved_resolved_to_ineligible_count,
    selected_flag_ids = selected$chosen_flag_ids,
    selected_cell_keys = sort(selected_cell_keys, method = "radix"),
    selected_cell_prior_state = prior_state_ledger,
    selected_cell_prior_state_sha256 = .lt_hash(prior_state_ledger),
    decision_sha256 = .lt_hash(decision_payload),
    source_diagnostic_spec_id = diagnostic$diagnostic_spec_id,
    source_diagnostic_outputs_sha256 = diagnostic$output_hashes$diagnostic_outputs_sha256,
    source_flag_ledger_sha256 = diagnostic$output_hashes$flag_ledger_sha256,
    source_threshold_config_sha256 = diagnostic$output_hashes$threshold_config_sha256,
    prior_status_sha256 = base$status_sha256,
    prior_reason_sha256 = base$reason_sha256,
    raw_matrix_modified = FALSE,
    value_operation = "none_no_rewrite_winsorize_smooth_clip_replace"
  )
  provenance <- base$provenance
  provenance$dataset_specific_exclusion_ledger_supplied <- TRUE
  if (!is.null(provenance$qc_application)) {
    history <- provenance$qc_application_history %||% list()
    history[[length(history) + 1L]] <- provenance$qc_application
    provenance$qc_application_history <- history
  }
  provenance$qc_application <- application
  out <- lt_measurement_eligibility(x, status, reason, provenance)
  if (!identical(before, .lt_qc_matrix_identity(x))) {
    .lt_abort("QC application mutated raw matrix identity.")
  }
  out
}

#' @export
print.lt_qc_spec <- function(x, ...) {
  validate_lt_qc_spec(x)
  cat("<lt_qc_spec> ", x$diagnostic_spec_id, "\n", sep = "")
  cat("  profile: ", x$profile, "\n", sep = "")
  cat("  behavior: DIAGNOSE_ONLY\n")
  cat("  outlier rules: ", length(x$outlier_rules), "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_diagnostic <- function(x, ...) {
  validate_lt_diagnostic(x)
  overall <- x$matrix_summary$overall
  cat("<lt_diagnostic> ", x$diagnostic_spec_id, "\n", sep = "")
  cat("  matrix: ", overall$gene_count, " genes x ", overall$branch_count, " branches\n", sep = "")
  cat("  available payloads: ", overall$payload_available_count, "\n", sep = "")
  cat("  flags: ", nrow(x$flags), " (diagnostic only)\n", sep = "")
  cat("  holds: ", nrow(x$holds), "\n", sep = "")
  invisible(x)
}
