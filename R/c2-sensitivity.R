# Positive-cell two-way log-OLS sensitivity estimand.

.lt_c2_baseline_reason_levels <- function() c(
  "available", "unsupported_gene", "unsupported_branch",
  "unsupported_gene_and_branch", "cross_component_baseline_unavailable",
  "numeric_indeterminate"
)

#' Construct the C2 positive-cell fit and evaluation domains
#'
#' C2 is a required alternative estimand for sensitivity. Its fit domain is
#' strictly positive; exact zero is never changed or given a pseudocount. Its
#' baseline-evaluation domain additionally includes eligible observed zeros,
#' but only when gene and branch effects are identified in the same positive
#' component.
#'
#' @param x An `lt_matrix`.
#' @param measurement_fit_eligibility Exact Layer-2 declaration.
#' @return An `lt_c2_domain`.
#' @export
lt_c2_domain <- function(x, measurement_fit_eligibility = NULL) {
  validate_lt_matrix(x)
  eligibility <- measurement_fit_eligibility %||%
    lt_measurement_eligibility(x)
  validate_lt_measurement_eligibility(eligibility, x)
  observed <- x$coord_state == "observed"
  admitted <- observed & eligibility$status == "eligible"
  if (any(admitted & !is.na(x$values) & x$values < 0)) {
    .lt_abort("C2 cannot admit a negative payload.")
  }
  evaluation_mask <- admitted & !is.na(x$values) & is.finite(x$values) &
    x$values >= 0
  fit_mask <- evaluation_mask & x$values > 0
  components <- .lt_components_mask(fit_mask, x$gene_ids, x$branch_ids)
  fit_linear <- which(fit_mask)
  nr <- nrow(x$values)
  fit_gene_index <- ((fit_linear - 1L) %% nr) + 1L
  fit_branch_index <- ((fit_linear - 1L) %/% nr) + 1L
  fit_component <- components$gene_component[fit_gene_index]
  if (any(fit_component == 0L) ||
      any(fit_component != components$branch_component[fit_branch_index])) {
    .lt_abort("C2 positive-edge component assignment is inconsistent.")
  }
  unsupported_gene_index <- which(!components$active_gene)
  unsupported_branch_index <- which(!components$active_branch)
  unsupported <- rbind(
    data.frame(
      vertex_type = rep("gene", length(unsupported_gene_index)),
      vertex_index = unsupported_gene_index,
      vertex_id = x$gene_ids[!components$active_gene],
      state = rep("UNSUPPORTED_VERTEX", length(unsupported_gene_index)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      vertex_type = rep("branch", length(unsupported_branch_index)),
      vertex_index = unsupported_branch_index,
      vertex_id = x$branch_ids[!components$active_branch],
      state = rep("UNSUPPORTED_VERTEX", length(unsupported_branch_index)),
      stringsAsFactors = FALSE
    )
  )
  hashes <- list(
    semantic_values_sha256 = .lt_hash(x$values),
    semantic_coordinate_state_sha256 = .lt_hash(x$coord_state),
    measurement_eligibility_sha256 = eligibility$status_sha256,
    measurement_reason_sha256 = eligibility$reason_sha256,
    measurement_provenance_sha256 = .lt_hash(eligibility$provenance),
    ordered_gene_ledger_sha256 = .lt_hash(x$gene_ids),
    ordered_branch_ledger_sha256 = .lt_hash(x$branch_ids),
    fit_domain_sha256 = .lt_rate_keyed_mask_hash(fit_mask, x$gene_ids, x$branch_ids),
    evaluation_domain_sha256 = .lt_rate_keyed_mask_hash(
      evaluation_mask, x$gene_ids, x$branch_ids
    ),
    positive_component_ledger_sha256 = .lt_hash(components$summary),
    unsupported_vertex_ledger_sha256 = .lt_hash(unsupported),
    coordinate_provenance_sha256 = .lt_hash(
      .lt_coordinate_provenance_semantic_payload(x$coordinate_provenance)
    )
  )
  hashes$domain_sha256 <- .lt_hash(hashes)
  out <- structure(list(
    matrix = x,
    eligibility = eligibility,
    fit_mask = fit_mask,
    evaluation_mask = evaluation_mask,
    fit_linear = fit_linear,
    fit_gene_index = fit_gene_index,
    fit_branch_index = fit_branch_index,
    fit_component = fit_component,
    fit_value = x$values[fit_linear],
    fit_log_value = log(x$values[fit_linear]),
    components = components,
    unsupported_vertices = unsupported,
    hashes = hashes,
    fit_domain_rule = paste(
      "coord_state observed AND Layer-2 eligible AND Y finite AND Y > 0"
    ),
    evaluation_domain_rule = paste(
      "coord_state observed AND Layer-2 eligible AND Y finite AND Y >= 0;",
      "baseline only within one identified positive component"
    ),
    estimand = "C2_positive_cell_log_OLS",
    pseudocount = "none"
  ), class = "lt_c2_domain")
  validate_lt_c2_domain(out)
  out
}

#' Validate a C2 domain
#'
#' @param x An `lt_c2_domain`.
#' @return `x`, invisibly.
#' @export
validate_lt_c2_domain <- function(x) {
  if (!inherits(x, "lt_c2_domain")) .lt_abort("`x` must inherit from lt_c2_domain.")
  validate_lt_matrix(x$matrix)
  validate_lt_measurement_eligibility(x$eligibility, x$matrix)
  if (!is.logical(x$fit_mask) || !is.logical(x$evaluation_mask) ||
      !identical(dim(x$fit_mask), dim(x$matrix$values)) ||
      !identical(dim(x$evaluation_mask), dim(x$matrix$values)) ||
      anyNA(x$fit_mask) || anyNA(x$evaluation_mask)) {
    .lt_abort("C2 fit/evaluation masks must be complete logical matrices.")
  }
  if (any(x$fit_mask & !x$evaluation_mask) ||
      any(x$matrix$values[x$fit_mask] <= 0) ||
      !identical(which(x$fit_mask), x$fit_linear) ||
      !identical(x$matrix$values[x$fit_linear], x$fit_value) ||
      !identical(log(x$fit_value), x$fit_log_value)) {
    .lt_abort("C2 positive fit ledger differs from its frozen domain rule.")
  }
  observed_hashes <- list(
    semantic_values_sha256 = .lt_hash(x$matrix$values),
    semantic_coordinate_state_sha256 = .lt_hash(x$matrix$coord_state),
    measurement_eligibility_sha256 = x$eligibility$status_sha256,
    measurement_reason_sha256 = x$eligibility$reason_sha256,
    measurement_provenance_sha256 = .lt_hash(x$eligibility$provenance),
    ordered_gene_ledger_sha256 = .lt_hash(x$matrix$gene_ids),
    ordered_branch_ledger_sha256 = .lt_hash(x$matrix$branch_ids),
    fit_domain_sha256 = .lt_rate_keyed_mask_hash(
      x$fit_mask, x$matrix$gene_ids, x$matrix$branch_ids
    ),
    evaluation_domain_sha256 = .lt_rate_keyed_mask_hash(
      x$evaluation_mask, x$matrix$gene_ids, x$matrix$branch_ids
    ),
    positive_component_ledger_sha256 = .lt_hash(x$components$summary),
    unsupported_vertex_ledger_sha256 = .lt_hash(x$unsupported_vertices),
    coordinate_provenance_sha256 = .lt_hash(
      .lt_coordinate_provenance_semantic_payload(x$matrix$coordinate_provenance)
    )
  )
  observed_hashes$domain_sha256 <- .lt_hash(observed_hashes)
  if (!identical(x$hashes, observed_hashes)) {
    .lt_abort("C2 domain hashes do not match its exact matrix/Layer-2 ledgers.")
  }
  invisible(x)
}

#' Deterministic C2 numerical control
#'
#' @param tolerance Maximum alternating-effect change.
#' @param max_iterations Maximum alternating sweeps.
#' @param normal_equation_tolerance Maximum absolute mean score residual.
#' @param replay Require an independent same-environment solver replay.
#' @return A named control list.
#' @export
lt_c2_control <- function(tolerance = 1e-10, max_iterations = 1000L,
                          normal_equation_tolerance = 4e-10,
                          replay = TRUE) {
  for (z in list(tolerance, normal_equation_tolerance)) {
    if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0) {
      .lt_abort("C2 tolerances must be positive finite scalars.")
    }
  }
  if (!is.numeric(max_iterations) || length(max_iterations) != 1L ||
      is.na(max_iterations) || max_iterations < 1 ||
      max_iterations != as.integer(max_iterations)) {
    .lt_abort("`max_iterations` must be a positive integer.")
  }
  if (!is.logical(replay) || length(replay) != 1L || is.na(replay) || !replay) {
    .lt_abort("C2 sensitivity requires deterministic replay.")
  }
  list(
    tolerance = as.double(tolerance),
    max_iterations = as.integer(max_iterations),
    normal_equation_tolerance = as.double(normal_equation_tolerance),
    replay = TRUE,
    solver = "deterministic_componentwise_alternating_conditional_means_v1",
    gauge = "mean_branch_log_effect_zero_within_component",
    objective = "positive_cell_log_scale_residual_sum_of_squares",
    pseudocount = "none"
  )
}

.lt_c2_group_sum <- function(index, weights, size) {
  grouped <- rowsum(
    matrix(weights, ncol = 1L), group = index,
    reorder = FALSE, na.rm = FALSE
  )
  out <- numeric(size)
  out[as.integer(rownames(grouped))] <- grouped[, 1L]
  out
}

.lt_c2_component_identity_sha256 <- function(domain, component_index) {
  comp <- domain$components$components[[component_index]]
  component_id <- domain$components$summary$component_id[[component_index]]
  .lt_hash(list(
    component_id = component_id,
    gene_ids = domain$matrix$gene_ids[comp$genes],
    branch_ids = domain$matrix$branch_ids[comp$branches],
    fit_domain_sha256 = .lt_rate_keyed_mask_hash(
      domain$fit_mask[comp$genes, comp$branches, drop = FALSE],
      domain$matrix$gene_ids[comp$genes],
      domain$matrix$branch_ids[comp$branches]
    )
  ))
}

.lt_c2_fit_semantic_payload <- function(x) {
  list(
    run_id = x$run_id,
    state = x$state,
    sensitivity_estimand = x$sensitivity_estimand,
    domain_sha256 = x$domain$hashes$domain_sha256,
    control_sha256 = .lt_hash(x$control),
    gene_effect_keys_sha256 = .lt_hash(x$gene_effect_keys),
    branch_effect_keys_sha256 = .lt_hash(x$branch_effect_keys),
    gene_log_effect_sha256 = .lt_hash(x$gene_log_effect),
    branch_log_effect_sha256 = .lt_hash(x$branch_log_effect),
    gene_factor_sha256 = .lt_hash(x$gene_factor),
    branch_factor_sha256 = .lt_hash(x$branch_factor),
    fit_mu_sha256 = .lt_hash(x$fit_mu),
    evaluation_linear_sha256 = .lt_hash(x$evaluation_linear),
    evaluation_mu_sha256 = .lt_hash(x$evaluation_mu),
    baseline_availability_sha256 = .lt_hash(x$baseline_availability),
    baseline_reason_sha256 = .lt_hash(x$baseline_reason),
    baseline_reason_vocabulary_sha256 = .lt_hash(x$baseline_reason_levels),
    certificates_sha256 = .lt_hash(x$certificates),
    deterministic_replay_pass = x$deterministic_replay_pass
  )
}

.lt_c2_solve_once <- function(domain, control) {
  ng <- nrow(domain$matrix$values)
  nb <- ncol(domain$matrix$values)
  alpha <- rep(NA_real_, ng)
  beta <- rep(NA_real_, nb)
  mu_fit <- rep(NA_real_, length(domain$fit_linear))
  certificates <- vector("list", length(domain$components$components))
  held <- FALSE
  for (i in seq_along(domain$components$components)) {
    comp <- domain$components$components[[i]]
    idx <- if (length(domain$components$components) == 1L) {
      seq_along(domain$fit_linear)
    } else which(domain$fit_component == i)
    gi <- domain$fit_gene_index[idx]
    bi <- domain$fit_branch_index[idx]
    ly <- domain$fit_log_value[idx]
    gd <- tabulate(gi, nbins = ng)
    bd <- tabulate(bi, nbins = nb)
    a <- rep(NA_real_, ng); b <- rep(NA_real_, nb)
    a[comp$genes] <- 0; b[comp$branches] <- 0
    converged <- FALSE; delta <- Inf; iteration <- 0L
    for (iteration in seq_len(control$max_iterations)) {
      previous_a <- a[comp$genes]
      previous_b <- b[comp$branches]
      a[comp$genes] <- .lt_c2_group_sum(
        gi, ly - b[bi], ng
      )[comp$genes] / gd[comp$genes]
      b[comp$branches] <- .lt_c2_group_sum(
        bi, ly - a[gi], nb
      )[comp$branches] / bd[comp$branches]
      shift <- mean(b[comp$branches])
      b[comp$branches] <- b[comp$branches] - shift
      a[comp$genes] <- a[comp$genes] + shift
      delta <- max(abs(a[comp$genes] - previous_a),
                   abs(b[comp$branches] - previous_b))
      if (is.finite(delta) && delta <= control$tolerance) {
        converged <- TRUE
        break
      }
    }
    residual <- ly - a[gi] - b[bi]
    gene_score <- .lt_c2_group_sum(gi, residual, ng)[comp$genes] /
      gd[comp$genes]
    branch_score <- .lt_c2_group_sum(bi, residual, nb)[comp$branches] /
      bd[comp$branches]
    normal_residual <- max(abs(gene_score), abs(branch_score))
    gauge_residual <- abs(mean(b[comp$branches]))
    fitted <- exp(a[gi] + b[bi])
    objective <- .lt_pairwise_sum(residual * residual)
    finite <- all(is.finite(a[comp$genes])) && all(is.finite(b[comp$branches])) &&
      all(is.finite(fitted)) && all(fitted > 0) && is.finite(objective)
    pass <- converged && finite &&
      normal_residual <= control$normal_equation_tolerance &&
      gauge_residual <= max(control$normal_equation_tolerance,
                            64 * .Machine$double.eps)
    component_id <- domain$components$summary$component_id[[i]]
    component_identity <- list(
      component_id = component_id,
      gene_ids = domain$matrix$gene_ids[comp$genes],
      branch_ids = domain$matrix$branch_ids[comp$branches],
      fit_domain_sha256 = .lt_rate_keyed_mask_hash(
        domain$fit_mask[comp$genes, comp$branches, drop = FALSE],
        domain$matrix$gene_ids[comp$genes],
        domain$matrix$branch_ids[comp$branches]
      )
    )
    certificates[[i]] <- list(
      component_id = component_id,
      component_identity_sha256 = .lt_hash(component_identity),
      component_membership_sha256 =
        domain$components$summary$membership_sha256[[i]],
      edge_count = length(idx),
      objective = objective,
      iterations = as.integer(iteration),
      maximum_absolute_effect_change = delta,
      maximum_absolute_normal_equation_mean_residual = normal_residual,
      gauge_residual = gauge_residual,
      fitted_value_sha256 = .lt_hash(fitted),
      solver_control_sha256 = .lt_hash(control),
      converged = converged,
      finite_positive_fitted_values = finite,
      certificate_pass = pass,
      state = if (pass) "FINITE_INTERIOR" else "NUMERICALLY_INDETERMINATE"
    )
    if (!pass) held <- TRUE
    alpha[comp$genes] <- a[comp$genes]
    beta[comp$branches] <- b[comp$branches]
    mu_fit[idx] <- fitted
  }
  if (held || anyNA(mu_fit)) {
    return(list(
      state = "NUMERICALLY_INDETERMINATE", alpha = alpha, beta = beta,
      mu_fit = NULL, certificates = certificates
    ))
  }
  list(state = "FINITE_INTERIOR", alpha = alpha, beta = beta,
       mu_fit = mu_fit, certificates = certificates)
}

#' Fit the C2 positive-cell log-OLS sensitivity baseline
#'
#' Fits each connected positive component independently and evaluates a
#' baseline only for supported gene/branch pairs in that same component. No
#' cross-component outer product is materialized.
#'
#' @param domain An `lt_c2_domain`.
#' @param control A list from `lt_c2_control()`.
#' @param run_id Stable sensitivity-run identifier.
#' @return An `lt_c2_fit`.
#' @export
lt_c2_fit <- function(domain, control = lt_c2_control(), run_id = "c2_sensitivity") {
  validate_lt_c2_domain(domain)
  .lt_assert_scalar_character(run_id, "run_id")
  expected_control <- do.call(lt_c2_control, control[c(
    "tolerance", "max_iterations", "normal_equation_tolerance", "replay"
  )])
  if (!identical(control, expected_control)) {
    .lt_abort("C2 control does not match its normalized frozen schema.")
  }
  fit <- .lt_c2_solve_once(domain, control)
  replay <- .lt_c2_solve_once(domain, control)
  replay_pass <- identical(fit$state, replay$state) &&
    identical(.lt_hash(fit$alpha), .lt_hash(replay$alpha)) &&
    identical(.lt_hash(fit$beta), .lt_hash(replay$beta)) &&
    identical(.lt_hash(fit$mu_fit), .lt_hash(replay$mu_fit)) &&
    identical(.lt_hash(fit$certificates), .lt_hash(replay$certificates))
  fit$certificates <- lapply(fit$certificates, function(certificate) {
    certificate$deterministic_replay_pass <- replay_pass
    certificate
  })
  if (!replay_pass || identical(fit$state, "NUMERICALLY_INDETERMINATE")) {
    out <- structure(list(
      run_id = run_id, state = "NUMERICALLY_INDETERMINATE",
      domain = domain, control = control,
      gene_effect_keys = domain$matrix$gene_ids,
      branch_effect_keys = domain$matrix$branch_ids,
      gene_log_effect = fit$alpha, branch_log_effect = fit$beta,
      gene_factor = exp(fit$alpha), branch_factor = exp(fit$beta),
      fit_mu = NULL, evaluation_linear = which(domain$evaluation_mask),
      evaluation_mu = NULL, baseline_availability = NULL,
      baseline_reason = NULL, baseline_reason_levels = .lt_c2_baseline_reason_levels(),
      certificates = fit$certificates,
      deterministic_replay_pass = replay_pass,
      sensitivity_estimand = "C2_positive_cell_log_OLS",
      scientific_pass = FALSE
    ), class = "lt_c2_fit")
    out$gene_log_effect_sha256 <- .lt_hash(out$gene_log_effect)
    out$branch_log_effect_sha256 <- .lt_hash(out$branch_log_effect)
    out$gene_factor_sha256 <- .lt_hash(out$gene_factor)
    out$branch_factor_sha256 <- .lt_hash(out$branch_factor)
    out$fit_mu_sha256 <- .lt_hash(out$fit_mu)
    out$evaluation_mu_sha256 <- .lt_hash(out$evaluation_mu)
    out$baseline_availability_sha256 <- .lt_hash(out$baseline_availability)
    out$baseline_reason_sha256 <- .lt_hash(out$baseline_reason)
    out$semantic_fit_sha256 <- .lt_hash(.lt_c2_fit_semantic_payload(out))
    validate_lt_c2_fit(out)
    return(out)
  }
  evaluation_linear <- which(domain$evaluation_mask)
  nr <- nrow(domain$matrix$values)
  eg <- ((evaluation_linear - 1L) %% nr) + 1L
  eb <- ((evaluation_linear - 1L) %/% nr) + 1L
  gc <- domain$components$gene_component[eg]
  bc <- domain$components$branch_component[eb]
  available <- gc > 0L & bc > 0L & gc == bc
  reason <- rep("available", length(evaluation_linear))
  reason[gc == 0L & bc > 0L] <- "unsupported_gene"
  reason[gc > 0L & bc == 0L] <- "unsupported_branch"
  reason[gc == 0L & bc == 0L] <- "unsupported_gene_and_branch"
  reason[gc > 0L & bc > 0L & gc != bc] <-
    "cross_component_baseline_unavailable"
  evaluation_mu <- rep(NA_real_, length(evaluation_linear))
  evaluation_mu[available] <- exp(fit$alpha[eg[available]] + fit$beta[eb[available]])
  if (any(!is.finite(evaluation_mu[available])) ||
      any(evaluation_mu[available] <= 0)) {
    .lt_abort("C2 evaluation baseline must be finite and positive on identified cells.")
  }
  baseline_reason <- match(reason, .lt_c2_baseline_reason_levels())
  out <- structure(list(
    run_id = run_id,
    state = "FINITE_INTERIOR",
    domain = domain,
    control = control,
    gene_effect_keys = domain$matrix$gene_ids,
    branch_effect_keys = domain$matrix$branch_ids,
    gene_log_effect = fit$alpha,
    branch_log_effect = fit$beta,
    gene_factor = exp(fit$alpha),
    branch_factor = exp(fit$beta),
    fit_mu = fit$mu_fit,
    evaluation_linear = evaluation_linear,
    evaluation_mu = evaluation_mu,
    baseline_availability = available,
    baseline_reason = as.integer(baseline_reason),
    baseline_reason_levels = .lt_c2_baseline_reason_levels(),
    certificates = fit$certificates,
    deterministic_replay_pass = TRUE,
    fit_mu_sha256 = .lt_hash(fit$mu_fit),
    evaluation_mu_sha256 = .lt_hash(evaluation_mu),
    gene_factor_sha256 = .lt_hash(exp(fit$alpha)),
    branch_factor_sha256 = .lt_hash(exp(fit$beta)),
    baseline_availability_sha256 = .lt_hash(available),
    baseline_reason_sha256 = .lt_hash(as.integer(baseline_reason)),
    sensitivity_estimand = "C2_positive_cell_log_OLS",
    factor_gauge_warning = "factor_values_are_gauge_dependent_within_component",
    scientific_pass = FALSE
  ), class = "lt_c2_fit")
  out$gene_log_effect_sha256 <- .lt_hash(out$gene_log_effect)
  out$branch_log_effect_sha256 <- .lt_hash(out$branch_log_effect)
  out$semantic_fit_sha256 <- .lt_hash(.lt_c2_fit_semantic_payload(out))
  validate_lt_c2_fit(out)
  out
}

#' Validate a C2 fit
#'
#' The validator recomputes effect/factor/baseline hashes, component certificate
#' identities, and one top-level semantic fit identity. The numerical solver is
#' not rerun by validation.
#'
#' @param x An `lt_c2_fit`.
#' @return `x`, invisibly.
#' @export
validate_lt_c2_fit <- function(x) {
  if (!inherits(x, "lt_c2_fit")) .lt_abort("`x` must inherit from lt_c2_fit.")
  validate_lt_c2_domain(x$domain)
  .lt_assert_scalar_character(x$run_id, "run_id")
  if (!x$state %in% c("FINITE_INTERIOR", "NUMERICALLY_INDETERMINATE")) {
    .lt_abort("Unknown C2 fit state.")
  }
  if (!identical(x$sensitivity_estimand, "C2_positive_cell_log_OLS")) {
    .lt_abort("C2 fit sensitivity estimand identity is invalid.")
  }
  normalized_control <- do.call(lt_c2_control, x$control[c(
    "tolerance", "max_iterations", "normal_equation_tolerance", "replay"
  )])
  if (!identical(x$control, normalized_control)) {
    .lt_abort("C2 fit control is not the normalized frozen control object.")
  }
  ng <- nrow(x$domain$matrix$values); nb <- ncol(x$domain$matrix$values)
  if (!identical(x$gene_effect_keys, x$domain$matrix$gene_ids) ||
      !identical(x$branch_effect_keys, x$domain$matrix$branch_ids) ||
      !is.numeric(x$gene_log_effect) || length(x$gene_log_effect) != ng ||
      !is.numeric(x$branch_log_effect) || length(x$branch_log_effect) != nb ||
      !is.numeric(x$gene_factor) || length(x$gene_factor) != ng ||
      !is.numeric(x$branch_factor) || length(x$branch_factor) != nb) {
    .lt_abort("C2 effect/factor vectors must align to exact scientific-key axes.")
  }
  if (!identical(x$gene_factor, exp(x$gene_log_effect)) ||
      !identical(x$branch_factor, exp(x$branch_log_effect)) ||
      !identical(x$gene_log_effect_sha256, .lt_hash(x$gene_log_effect)) ||
      !identical(x$branch_log_effect_sha256, .lt_hash(x$branch_log_effect)) ||
      !identical(x$gene_factor_sha256, .lt_hash(x$gene_factor)) ||
      !identical(x$branch_factor_sha256, .lt_hash(x$branch_factor))) {
    .lt_abort("C2 effect/factor semantic identity is stale or inconsistent.")
  }
  component_count <- length(x$domain$components$components)
  if (!is.list(x$certificates) || length(x$certificates) != component_count) {
    .lt_abort("C2 certificates must be one-to-one with domain components.")
  }
  for (i in seq_len(component_count)) {
    certificate <- x$certificates[[i]]
    if (!identical(certificate$component_id,
                   x$domain$components$summary$component_id[[i]]) ||
        !identical(certificate$component_membership_sha256,
                   x$domain$components$summary$membership_sha256[[i]]) ||
        !identical(certificate$component_identity_sha256,
                   .lt_c2_component_identity_sha256(x$domain, i)) ||
        !identical(certificate$solver_control_sha256, .lt_hash(x$control)) ||
        !identical(certificate$deterministic_replay_pass,
                   x$deterministic_replay_pass)) {
      .lt_abort("C2 certificate component/control/replay identity is inconsistent.")
    }
  }
  if (!isTRUE(x$deterministic_replay_pass) && x$state != "NUMERICALLY_INDETERMINATE") {
    .lt_abort("Released C2 fit requires deterministic replay identity.")
  }
  if (x$state == "NUMERICALLY_INDETERMINATE") {
    if (!is.null(x$fit_mu) || !is.null(x$evaluation_mu)) {
      .lt_abort("Numerically indeterminate C2 fit must not release baseline values.")
    }
    if (!identical(x$fit_mu_sha256, .lt_hash(x$fit_mu)) ||
        !identical(x$evaluation_mu_sha256, .lt_hash(x$evaluation_mu)) ||
        !identical(x$baseline_availability_sha256,
                   .lt_hash(x$baseline_availability)) ||
        !identical(x$baseline_reason_sha256, .lt_hash(x$baseline_reason)) ||
        !identical(x$semantic_fit_sha256,
                   .lt_hash(.lt_c2_fit_semantic_payload(x)))) {
      .lt_abort("Indeterminate C2 fit semantic identity is stale or inconsistent.")
    }
    return(invisible(x))
  }
  if (!all(vapply(x$certificates, `[[`, logical(1L), "certificate_pass"))) {
    .lt_abort("Released C2 fit has a failed component certificate.")
  }
  if (!is.numeric(x$fit_mu) || length(x$fit_mu) != length(x$domain$fit_linear) ||
      any(!is.finite(x$fit_mu)) || any(x$fit_mu <= 0) ||
      !identical(x$fit_mu_sha256, .lt_hash(x$fit_mu))) {
    .lt_abort("C2 fit baseline does not match its positive-domain hash.")
  }
  if (!is.numeric(x$evaluation_mu) ||
      length(x$evaluation_mu) != length(x$evaluation_linear) ||
      !is.logical(x$baseline_availability) ||
      !identical(x$baseline_availability, !is.na(x$evaluation_mu)) ||
      !identical(x$evaluation_mu_sha256, .lt_hash(x$evaluation_mu)) ||
      !identical(x$baseline_availability_sha256, .lt_hash(x$baseline_availability)) ||
      !identical(x$baseline_reason_sha256, .lt_hash(x$baseline_reason))) {
    .lt_abort("C2 evaluation baseline availability/reason hashes are inconsistent.")
  }
  if (!identical(x$evaluation_linear, which(x$domain$evaluation_mask)) ||
      !identical(x$baseline_reason_levels, .lt_c2_baseline_reason_levels()) ||
      !identical(x$semantic_fit_sha256,
                 .lt_hash(.lt_c2_fit_semantic_payload(x)))) {
    .lt_abort("C2 top-level semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

.lt_c2_identity_gate <- function(x, eligibility, fit) {
  validate_lt_c2_fit(fit)
  validate_lt_matrix(x)
  validate_lt_measurement_eligibility(eligibility, x)
  d <- fit$domain
  if (!identical(x$gene_ids, d$matrix$gene_ids) ||
      !identical(x$branch_ids, d$matrix$branch_ids) ||
      !identical(.lt_hash(x$values), d$hashes$semantic_values_sha256) ||
      !identical(.lt_hash(x$coord_state), d$hashes$semantic_coordinate_state_sha256) ||
      !identical(eligibility$status_sha256,
                 d$hashes$measurement_eligibility_sha256) ||
      !identical(eligibility$reason_sha256, d$hashes$measurement_reason_sha256) ||
      !identical(.lt_hash(eligibility$provenance),
                 d$hashes$measurement_provenance_sha256) ||
      !identical(
        .lt_hash(.lt_coordinate_provenance_semantic_payload(x$coordinate_provenance)),
        d$hashes$coordinate_provenance_sha256
      )) {
    .lt_abort("C2 sensitivity input/axis/Layer-2 identity mismatch.")
  }
  invisible(TRUE)
}

#' Construct a C2 sensitivity representation
#'
#' @param x Exact input `lt_matrix`.
#' @param c2_fit Certified `lt_c2_fit`.
#' @param representation One of `ADD_C2`, `GBI_C2`, or `logGBI_C2`.
#' @param measurement_fit_eligibility Exact Layer-2 declaration used by C2.
#' @param metadata Optional named metadata.
#' @return An `lt_rate` labelled as C2 sensitivity, never production.
#' @export
lt_c2_sensitivity <- function(
    x, c2_fit, representation = c("ADD_C2", "GBI_C2", "logGBI_C2"),
    measurement_fit_eligibility = c2_fit$domain$eligibility,
    metadata = list()) {
  representation <- match.arg(representation)
  .lt_assert_named_list(metadata, "metadata")
  .lt_c2_identity_gate(x, measurement_fit_eligibility, c2_fit)
  if (identical(c2_fit$state, "NUMERICALLY_INDETERMINATE")) {
    .lt_abort("Required C2 sensitivity component is NUMERICALLY_INDETERMINATE.")
  }
  ncell <- length(x$values); dn <- dimnames(x$values)
  reason <- rep("coordinate_unavailable", ncell)
  observed <- x$coord_state == "observed"
  reason[observed & measurement_fit_eligibility$status == "ineligible"] <-
    "measurement_ineligible"
  reason[observed & measurement_fit_eligibility$status == "unresolved"] <-
    "measurement_unresolved"
  reason[observed & measurement_fit_eligibility$status == "eligible" &
           is.na(x$values)] <- "observed_payload_unavailable"
  support <- rep("BASELINE_UNAVAILABLE", ncell)
  ev <- c2_fit$evaluation_linear
  baseline_available <- c2_fit$baseline_availability
  support[ev[baseline_available]] <- "C2_SUPPORTED"
  baseline_reason <- c2_fit$baseline_reason_levels[c2_fit$baseline_reason]
  support[ev[!baseline_available & grepl("unsupported", baseline_reason)]] <-
    "C2_UNSUPPORTED_VERTEX"
  support[ev[!baseline_available & grepl("cross_component", baseline_reason)]] <-
    "C2_CROSS_COMPONENT"
  reason[ev[!baseline_available]] <- ifelse(
    baseline_reason[!baseline_available] == "cross_component_baseline_unavailable",
    "cross_component_baseline_unavailable", "baseline_unavailable"
  )
  y <- x$values[ev]
  mu <- c2_fit$evaluation_mu
  available <- baseline_available
  if (representation == "logGBI_C2") {
    available <- available & y > 0
    reason[ev[baseline_available & y == 0]] <- "input_zero_log_undefined"
  }
  reason[ev[available]] <- "available"
  values <- rep(NA_real_, ncell)
  if (representation == "ADD_C2") {
    values[ev[available]] <- y[available] - mu[available]
  } else if (representation == "GBI_C2") {
    values[ev[available]] <- y[available] / mu[available]
  } else {
    values[ev[available]] <- log(y[available] / mu[available])
  }
  if (any(is.infinite(values), na.rm = TRUE) || any(is.nan(values))) {
    .lt_abort("C2 sensitivity cannot store Inf, -Inf, or NaN.")
  }
  dim(values) <- dim(x$values); dimnames(values) <- dn
  dim(reason) <- dim(x$values); dimnames(reason) <- dn
  dim(support) <- dim(x$values); dimnames(support) <- dn
  method <- sub("_C2$", "", representation)
  units <- if (method == "ADD") x$payload$units else "dimensionless"
  sensitivity_spec <- structure(list(
    rate_id = paste0(c2_fit$run_id, "_", representation),
    method = representation,
    trait_independent = TRUE,
    production = FALSE,
    sensitivity_estimand = "C2_positive_cell_log_OLS",
    pseudocount = "none"
  ), class = "lt_c2_sensitivity_spec")
  .lt_new_rate(
    rate_spec = sensitivity_spec,
    representation = representation,
    values = values,
    x = x,
    eligibility = measurement_fit_eligibility,
    baseline_estimand = "C2_positive_cell_log_OLS",
    baseline_fit_id = c2_fit$run_id,
    baseline_mu_sha256 = c2_fit$evaluation_mu_sha256,
    baseline_mu_all_admitted_sha256 = c2_fit$fit_mu_sha256,
    baseline_support_labels = support,
    reason_labels = reason,
    representation_provenance = list(
      operator = paste0("LaTerra_T1B_sensitivity_", representation),
      sensitivity_estimand = "C2_positive_cell_log_OLS",
      production = FALSE, trait_independent = TRUE,
      pseudocount = "none", numeric_trimming = FALSE,
      units = units,
      factor_gauge_warning = c2_fit$factor_gauge_warning
    ),
    metadata = utils::modifyList(list(
      units = units, production = FALSE, sensitivity = TRUE,
      sensitivity_estimand = "C2_positive_cell_log_OLS",
      scientific_pass = FALSE
    ), metadata)
  )
}

#' @export
print.lt_c2_domain <- function(x, ...) {
  validate_lt_c2_domain(x)
  cat("<lt_c2_domain>\n")
  cat("  positive fit cells: ", length(x$fit_linear), "\n", sep = "")
  cat("  positive components: ", length(x$components$components), "\n", sep = "")
  cat("  unsupported vertices: ", nrow(x$unsupported_vertices), "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_c2_fit <- function(x, ...) {
  validate_lt_c2_fit(x)
  cat("<lt_c2_fit> ", x$run_id, "\n", sep = "")
  cat("  state: ", x$state, "\n", sep = "")
  cat("  sensitivity estimand: C2_positive_cell_log_OLS\n")
  cat("  fit cells: ", length(x$domain$fit_linear), "\n", sep = "")
  cat("  baseline-supported evaluation cells: ",
      if (is.null(x$baseline_availability)) 0 else sum(x$baseline_availability),
      "\n", sep = "")
  invisible(x)
}
