#' C3 support/existence states
#'
#' @return The frozen component, vertex, and cell state vocabulary used by the
#'   C3 baseline engine.
#' @export
lt_c3_existence_states <- function() {
  c(
    "FINITE_INTERIOR", "UNSUPPORTED_VERTEX", "ZERO_MARGIN_BOUNDARY",
    "ZERO_TOTAL_COMPONENT", "FACIAL_BOUNDARY", "NUMERICALLY_INDETERMINATE"
  )
}

.lt_eligibility_states <- function() {
  c("eligible", "ineligible", "unresolved", "not_applicable")
}
.lt_eligibility_reasons <- function() c(
  "eligible_declared", "analysis_domain_excluded", "taxon_domain_excluded",
  "measurement_qc_excluded", "user_predeclared_excluded", "authority_unresolved",
  "coordinate_not_observed"
)

#' Declare the C3 measurement-fit eligibility layer
#'
#' Eligibility is an input/declaration layer. It is not inferred from a
#' diagnostic and does not certify upstream sequence, alignment, tree, or
#' measurement quality. When `status` is omitted, every observed finite
#' nonnegative observed cell is declared eligible with reason
#' `eligible_declared`; other observed cells remain unresolved. Layer 2 is
#' `not_applicable` with reason `coordinate_not_observed` wherever coordinate
#' truth is not `observed`.
#'
#' @param x An `lt_matrix`.
#' @param status Optional character matrix containing `eligible`, `ineligible`,
#'   `unresolved`, or `not_applicable`.
#' @param reason Optional character reason matrix with the same axes.
#' @param provenance Named declaration provenance.
#' @return An `lt_measurement_eligibility` object.
#' @export
lt_measurement_eligibility <- function(x, status = NULL, reason = NULL,
                                       provenance = list()) {
  validate_lt_matrix(x)
  dn <- dimnames(x$values)
  if (is.null(status)) {
    observed <- x$coord_state == "observed"
    admissible <- observed & is.finite(x$values) &
      !is.na(x$values) & x$values >= 0
    status <- matrix("not_applicable", nrow(x$values), ncol(x$values), dimnames = dn)
    status[observed] <- "unresolved"
    status[admissible] <- "eligible"
    reason <- matrix("coordinate_not_observed", nrow(x$values), ncol(x$values),
                     dimnames = dn)
    reason[observed] <- "authority_unresolved"
    reason[admissible] <- "eligible_declared"
    default_provenance <- list(
      declaration = "observed_finite_nonnegative_default",
      upstream_qc_certified = FALSE,
      dataset_specific_exclusion_ledger_supplied = FALSE
    )
    provenance <- utils::modifyList(default_provenance, provenance)
  } else {
    if (is.null(reason)) {
      .lt_abort("An explicit eligibility `status` requires a `reason` matrix.")
    }
    default_provenance <- list(
      declaration = "caller_supplied_layer_2",
      upstream_qc_certified = FALSE,
      dataset_specific_exclusion_ledger_supplied = TRUE
    )
    provenance <- utils::modifyList(default_provenance, provenance)
  }
  .lt_assert_named_list(provenance, "provenance")
  if (is.matrix(status) && identical(dim(status), dim(x$values)) &&
      is.null(dimnames(status))) dimnames(status) <- dn
  if (is.matrix(reason) && identical(dim(reason), dim(x$values)) &&
      is.null(dimnames(reason))) dimnames(reason) <- dn
  out <- structure(
    list(
      status = status,
      reason = reason,
      provenance = provenance,
      status_sha256 = .lt_hash(status),
      reason_sha256 = .lt_hash(reason)
    ),
    class = "lt_measurement_eligibility"
  )
  validate_lt_measurement_eligibility(out, x)
  out
}

#' Validate a scoped measurement-fit eligibility layer
#'
#' @param x An `lt_measurement_eligibility` object.
#' @param matrix_object The aligned `lt_matrix` providing coordinate truth.
#' @return `x`, invisibly.
#' @export
validate_lt_measurement_eligibility <- function(x, matrix_object) {
  if (!inherits(x, "lt_measurement_eligibility")) {
    .lt_abort("`x` must inherit from lt_measurement_eligibility.")
  }
  validate_lt_matrix(matrix_object)
  status <- x$status
  reason <- x$reason
  dn <- dimnames(matrix_object$values)
  if (!is.matrix(status) || !is.character(status) ||
      !is.matrix(reason) || !is.character(reason)) {
    .lt_abort("Eligibility `status` and `reason` must be character matrices.")
  }
  if (!identical(dim(status), dim(matrix_object$values)) ||
      !identical(dim(reason), dim(matrix_object$values))) {
    .lt_abort("Matrix/state/eligibility dimensions must agree.")
  }
  if (!identical(dimnames(status), dn) || !identical(dimnames(reason), dn)) {
    .lt_abort("Eligibility axes must exactly match the ordered matrix axes.")
  }
  if (anyNA(status) || anyNA(reason) || any(!nzchar(status)) || any(!nzchar(reason))) {
    .lt_abort("Eligibility status and reason must be explicit non-missing labels.")
  }
  unknown_state <- setdiff(unique(as.vector(status)), .lt_eligibility_states())
  unknown_reason <- setdiff(unique(as.vector(reason)), .lt_eligibility_reasons())
  if (length(unknown_state)) {
    .lt_abort(sprintf("Unknown measurement eligibility labels: %s.",
                      paste(unknown_state, collapse = ", ")))
  }
  if (length(unknown_reason)) {
    .lt_abort(sprintf("Unknown measurement eligibility reasons: %s.",
                      paste(unknown_reason, collapse = ", ")))
  }
  observed <- matrix_object$coord_state == "observed"
  if (any(!observed & status != "not_applicable")) {
    .lt_abort("Non-observed coordinates require Layer-2 status `not_applicable`.")
  }
  if (any(observed & status == "not_applicable")) {
    .lt_abort("Observed coordinates cannot use Layer-2 status `not_applicable`.")
  }
  if (any(status == "not_applicable" & reason != "coordinate_not_observed")) {
    .lt_abort("Not-applicable cells require reason `coordinate_not_observed`.")
  }
  if (any(status != "not_applicable" & reason == "coordinate_not_observed")) {
    .lt_abort("Reason `coordinate_not_observed` is exclusive to not-applicable cells.")
  }
  if (any(status == "eligible" & reason != "eligible_declared")) {
    .lt_abort("Eligible cells must use reason `eligible_declared`.")
  }
  if (any(status == "unresolved" & reason != "authority_unresolved")) {
    .lt_abort("Unresolved cells must use reason `authority_unresolved`.")
  }
  if (any(status == "ineligible" & reason %in%
          c("eligible_declared", "authority_unresolved", "coordinate_not_observed"))) {
    .lt_abort("Ineligible cells require an explicit exclusion reason.")
  }
  if (!identical(x$status_sha256, .lt_hash(status)) ||
      !identical(x$reason_sha256, .lt_hash(reason))) {
    .lt_abort("Layer-2 status/reason hashes do not match their aligned matrices.")
  }
  invisible(x)
}

.lt_payload_state <- function(values, coord_state) {
  out <- matrix("coordinate_no_payload", nrow(values), ncol(values),
                dimnames = dimnames(values))
  observed <- coord_state == "observed"
  out[observed & is.na(values)] <- "nonfinite"
  out[observed & !is.na(values) & !is.finite(values)] <- "nonfinite"
  out[observed & is.finite(values) & values < 0] <- "negative"
  out[observed & is.finite(values) & values == 0] <- "exact_zero"
  out[observed & is.finite(values) & values > 0] <- "finite_positive"
  out
}

#' Construct the declared C3 fit domain
#'
#' @param x An `lt_matrix`.
#' @param measurement_fit_eligibility Optional
#'   `lt_measurement_eligibility`; a transparent declaration is constructed
#'   when omitted.
#' @param diagnostic_warnings Character diagnostics retained for audit only.
#' @return An `lt_c3_domain` object.
#' @export
lt_c3_domain <- function(x, measurement_fit_eligibility = NULL,
                         diagnostic_warnings = character()) {
  validate_lt_matrix(x)
  eligibility <- measurement_fit_eligibility %||%
    lt_measurement_eligibility(x)
  validate_lt_measurement_eligibility(eligibility, x)
  if (!is.character(diagnostic_warnings) || anyNA(diagnostic_warnings)) {
    .lt_abort("`diagnostic_warnings` must be a non-missing character vector.")
  }
  payload_state <- .lt_payload_state(x$values, x$coord_state)
  admitted_by_labels <- x$coord_state == "observed" &
    eligibility$status == "eligible"
  if (any(admitted_by_labels & payload_state == "negative")) {
    .lt_abort("Admitted negative payload terminates C3 domain validation.")
  }
  if (any(admitted_by_labels & payload_state == "nonfinite")) {
    .lt_abort("Admitted nonfinite payload terminates C3 domain validation.")
  }
  fit_mask <- admitted_by_labels &
    payload_state %in% c("finite_positive", "exact_zero")
  edge_linear <- which(fit_mask)
  nr <- nrow(x$values)
  edge_gene_index <- ((edge_linear - 1L) %% nr) + 1L
  edge_branch_index <- ((edge_linear - 1L) %/% nr) + 1L
  authority_hash <- function(candidate, fallback) {
    if (is.character(candidate) && length(candidate) == 1L && !is.na(candidate) &&
        grepl("^[0-9a-f]{64}$", candidate)) candidate else fallback
  }
  semantic_values_hash <- .lt_hash(x$values)
  semantic_state_hash <- .lt_hash(x$coord_state)
  semantic_gene_hash <- .lt_hash(x$gene_ids)
  semantic_branch_hash <- .lt_hash(x$branch_ids)
  hashes <- list(
    input_values_sha256 = authority_hash(
      x$metadata$source_values_sha256, semantic_values_hash
    ),
    coordinate_state_sha256 = authority_hash(
      x$metadata$source_state_sha256, semantic_state_hash
    ),
    measurement_eligibility_sha256 = .lt_hash(eligibility$status),
    measurement_reason_sha256 = .lt_hash(eligibility$reason),
    ordered_gene_ledger_sha256 = authority_hash(
      x$metadata$ordered_gene_ledger_sha256, semantic_gene_hash
    ),
    ordered_branch_ledger_sha256 = authority_hash(
      x$coordinate_provenance$ordered_branch_ledger_sha256, semantic_branch_hash
    ),
    semantic_values_sha256 = semantic_values_hash,
    semantic_coordinate_state_sha256 = semantic_state_hash,
    semantic_ordered_gene_ledger_sha256 = semantic_gene_hash,
    semantic_ordered_branch_ledger_sha256 = semantic_branch_hash,
    ordered_cell_ledger_sha256 = .lt_hash(list(
      gene_ids = x$gene_ids, branch_ids = x$branch_ids,
      edge_linear = edge_linear
    ))
  )
  out <- structure(
    list(
      matrix = x,
      eligibility = eligibility,
      payload_state = payload_state,
      fit_mask = fit_mask,
      edge_linear = edge_linear,
      edge_gene_index = edge_gene_index,
      edge_branch_index = edge_branch_index,
      edge_value = x$values[edge_linear],
      diagnostic_warnings = diagnostic_warnings,
      hashes = hashes,
      domain_rule = paste(
        "coord_state == observed AND measurement_fit_eligibility == eligible",
        "AND Y finite AND Y >= 0"
      )
    ),
    class = "lt_c3_domain"
  )
  validate_lt_c3_domain(out)
  out
}

#' Validate a C3 fit-domain object
#'
#' @param x An `lt_c3_domain`.
#' @return `x`, invisibly.
#' @export
validate_lt_c3_domain <- function(x) {
  if (!inherits(x, "lt_c3_domain")) .lt_abort("`x` must inherit from lt_c3_domain.")
  validate_lt_matrix(x$matrix)
  validate_lt_measurement_eligibility(x$eligibility, x$matrix)
  if (!is.logical(x$fit_mask) || !identical(dim(x$fit_mask), dim(x$matrix$values))) {
    .lt_abort("C3 fit mask must be a logical matrix aligned to the input matrix.")
  }
  if (anyNA(x$fit_mask)) .lt_abort("C3 fit mask must not contain NA.")
  if (!identical(which(x$fit_mask), x$edge_linear)) {
    .lt_abort("C3 ordered edge ledger does not match the fit mask.")
  }
  if (any(x$edge_value < 0) || any(!is.finite(x$edge_value))) {
    .lt_abort("C3 edge payloads must be finite and nonnegative.")
  }
  invisible(x)
}

#' C3 deterministic numerical control
#'
#' @param fit_rtol Relative fitted-margin tolerance.
#' @param fit_atol_scale Absolute tolerance multiplier.
#' @param max_iter Maximum IPF sweeps.
#' @param certificate_max_edges Deterministic existence-certificate resource
#'   ceiling. Exceeding it fails closed as `NUMERICALLY_INDETERMINATE`.
#' @param replay Whether to rerun the numerical solver and compare hashes.
#' @return A named control list.
#' @export
lt_c3_control <- function(fit_rtol = 1e-10, fit_atol_scale = 1e-12,
                          max_iter = 10000L,
                          certificate_max_edges = 2e7,
                          replay = TRUE) {
  if (!is.numeric(fit_rtol) || length(fit_rtol) != 1L ||
      !is.finite(fit_rtol) || fit_rtol <= 0) .lt_abort("Invalid `fit_rtol`.")
  if (!is.numeric(fit_atol_scale) || length(fit_atol_scale) != 1L ||
      !is.finite(fit_atol_scale) || fit_atol_scale <= 0) .lt_abort("Invalid `fit_atol_scale`.")
  if (!is.numeric(max_iter) || length(max_iter) != 1L || is.na(max_iter) ||
      max_iter < 1) .lt_abort("Invalid `max_iter`.")
  if (!is.numeric(certificate_max_edges) || length(certificate_max_edges) != 1L ||
      is.na(certificate_max_edges) || certificate_max_edges < 0) {
    .lt_abort("Invalid `certificate_max_edges`.")
  }
  if (!is.logical(replay) || length(replay) != 1L || is.na(replay)) {
    .lt_abort("`replay` must be TRUE or FALSE.")
  }
  if (!replay) {
    .lt_abort("Production C3 control requires deterministic replay.")
  }
  list(
    fit_rtol = as.double(fit_rtol),
    fit_atol_scale = as.double(fit_atol_scale),
    max_iter = as.integer(max_iter),
    certificate_max_edges = as.double(certificate_max_edges),
    replay = replay,
    summation_algorithm = "compensated_adjacent_pairwise_twosum_v1",
    canonical_accumulation_order = "branch_outer_gene_inner_column_major",
    existence_algorithm = "exact_residual_scc_facial_reduction_v1"
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.lt_components_mask <- function(mask, gene_ids, branch_ids) {
  active_g <- rowSums(mask) > 0
  active_b <- colSums(mask) > 0
  seen_g <- logical(nrow(mask))
  seen_b <- logical(ncol(mask))
  comps <- list()
  while (any(active_g & !seen_g) || any(active_b & !seen_b)) {
    g <- logical(nrow(mask)); b <- logical(ncol(mask))
    ug <- which(active_g & !seen_g)
    if (length(ug)) g[ug[[1L]]] <- TRUE else b[which(active_b & !seen_b)[[1L]]] <- TRUE
    repeat {
      b2 <- b
      if (any(g)) b2 <- b2 | (active_b & colSums(mask[g, , drop = FALSE]) > 0)
      g2 <- g
      if (any(b2)) g2 <- g2 | (active_g & rowSums(mask[, b2, drop = FALSE]) > 0)
      if (identical(g2, g) && identical(b2, b)) break
      g <- g2; b <- b2
    }
    comps[[length(comps) + 1L]] <- list(genes = which(g), branches = which(b))
    seen_g <- seen_g | g; seen_b <- seen_b | b
  }
  if (!length(comps)) {
    return(list(
      components = list(), gene_component = integer(nrow(mask)),
      branch_component = integer(ncol(mask)), active_gene = active_g,
      active_branch = active_b, summary = data.frame()
    ))
  }
  keys <- vapply(comps, function(z) paste(
    if (length(z$genes)) min(gene_ids[z$genes]) else "",
    if (length(z$branches)) min(branch_ids[z$branches]) else "",
    paste(gene_ids[z$genes], collapse = "\r"),
    paste(branch_ids[z$branches], collapse = "\r"), sep = "\n"
  ), character(1))
  comps <- comps[order(keys, method = "radix")]
  gcid <- integer(nrow(mask)); bcid <- integer(ncol(mask))
  summary <- lapply(seq_along(comps), function(i) {
    z <- comps[[i]]
    gcid[z$genes] <<- i; bcid[z$branches] <<- i
    data.frame(
      component_id = sprintf("component_%03d", i),
      gene_count = length(z$genes), branch_count = length(z$branches),
      edge_count = sum(mask[z$genes, z$branches, drop = FALSE]),
      smallest_gene_id = min(gene_ids[z$genes]),
      smallest_branch_id = min(branch_ids[z$branches]),
      membership_sha256 = .lt_hash(list(gene_ids[z$genes], branch_ids[z$branches],
                                        which(mask[z$genes, z$branches, drop = FALSE]))),
      stringsAsFactors = FALSE
    )
  })
  list(
    components = comps, gene_component = gcid, branch_component = bcid,
    active_gene = active_g, active_branch = active_b,
    summary = do.call(rbind, summary)
  )
}

.lt_matrix_margins <- function(values, mask) {
  gene <- vapply(seq_len(nrow(values)), function(i) {
    .lt_pairwise_sum(values[i, mask[i, ], drop = TRUE])
  }, numeric(1))
  branch <- vapply(seq_len(ncol(values)), function(j) {
    .lt_pairwise_sum(values[mask[, j], j, drop = TRUE])
  }, numeric(1))
  total <- .lt_pairwise_sum(values[which(mask)])
  n_edges <- sum(mask)
  sum_atol <- 64 * .Machine$double.eps * max(1, abs(total)) *
    max(1, ceiling(log2(max(1, n_edges))))
  gene_residual <- abs(.lt_pairwise_sum(gene) - total)
  branch_residual <- abs(.lt_pairwise_sum(branch) - total)
  list(
    gene = gene, branch = branch, total = total, sum_atol = sum_atol,
    gene_total_residual = gene_residual,
    branch_total_residual = branch_residual,
    certified = is.finite(total) && all(is.finite(gene)) && all(gene >= 0) &&
      all(is.finite(branch)) && all(branch >= 0) &&
      gene_residual <= sum_atol && branch_residual <= sum_atol
  )
}

.lt_component_margin_certificates <- function(values, mask, components, margins) {
  if (!length(components$components)) return(data.frame())
  do.call(rbind, lapply(seq_along(components$components), function(i) {
    z <- components$components[[i]]
    local <- mask[z$genes, z$branches, drop = FALSE]
    n_edges <- sum(local)
    ll <- which(local)
    lr <- ((ll - 1L) %% length(z$genes)) + 1L
    lc <- ((ll - 1L) %/% length(z$genes)) + 1L
    global_linear <- z$genes[lr] + (z$branches[lc] - 1L) * nrow(values)
    total <- .lt_pairwise_sum(values[global_linear])
    gr <- .lt_pairwise_sum(margins$gene[z$genes])
    bc <- .lt_pairwise_sum(margins$branch[z$branches])
    atol <- 64 * .Machine$double.eps * max(1, abs(total)) *
      max(1, ceiling(log2(max(1, n_edges))))
    data.frame(
      component_id = sprintf("component_%03d", i),
      edge_count = n_edges, total = total,
      gene_total_residual = abs(gr - total),
      branch_total_residual = abs(bc - total), sum_atol = atol,
      finite_nonnegative_margins = all(is.finite(margins$gene[z$genes])) &&
        all(margins$gene[z$genes] >= 0) &&
        all(is.finite(margins$branch[z$branches])) &&
        all(margins$branch[z$branches] >= 0),
      pass = abs(gr - total) <= atol && abs(bc - total) <= atol,
      stringsAsFactors = FALSE
    )
  }))
}

.lt_scc <- function(n, from, to) {
  if (n == 0L) return(integer())
  if (!length(from)) return(seq_len(n))
  pairs <- unique(data.frame(from = as.integer(from), to = as.integer(to)))
  adj <- split(pairs$to, factor(pairs$from, levels = seq_len(n)))
  radj <- split(pairs$from, factor(pairs$to, levels = seq_len(n)))
  seen <- logical(n); finish <- integer(n); nf <- 0L
  for (root in seq_len(n)) if (!seen[root]) {
    sv <- root; si <- 1L; seen[root] <- TRUE
    while (length(sv)) {
      k <- length(sv); v <- sv[[k]]; nei <- adj[[v]]; pos <- si[[k]]
      if (pos <= length(nei)) {
        w <- nei[[pos]]; si[[k]] <- pos + 1L
        if (!seen[w]) { seen[w] <- TRUE; sv <- c(sv, w); si <- c(si, 1L) }
      } else {
        nf <- nf + 1L; finish[[nf]] <- v
        sv <- sv[-k]; si <- si[-k]
      }
    }
  }
  comp <- integer(n); cid <- 0L
  for (root in rev(finish)) if (comp[root] == 0L) {
    cid <- cid + 1L; stack <- root; comp[root] <- cid
    while (length(stack)) {
      v <- stack[[length(stack)]]; stack <- stack[-length(stack)]
      for (w in radj[[v]]) if (comp[w] == 0L) {
        comp[w] <- cid; stack <- c(stack, w)
      }
    }
  }
  comp
}

.lt_condensation_rank <- function(n_scc, from_scc, to_scc) {
  if (n_scc <= 1L) return(0L)
  keep <- from_scc != to_scc
  pairs <- unique(data.frame(from = from_scc[keep], to = to_scc[keep]))
  indeg <- tabulate(pairs$to, nbins = n_scc)
  adj <- split(pairs$to, factor(pairs$from, levels = seq_len(n_scc)))
  rank <- integer(n_scc)
  queue <- which(indeg == 0L)
  done <- 0L
  while (length(queue)) {
    v <- queue[[1L]]; queue <- queue[-1L]; done <- done + 1L
    for (w in adj[[v]]) {
      rank[w] <- max(rank[w], rank[v] + 1L)
      indeg[w] <- indeg[w] - 1L
      if (indeg[w] == 0L) queue <- sort(c(queue, w))
    }
  }
  if (done != n_scc) .lt_abort("Internal error: SCC condensation is not acyclic.")
  rank
}

.lt_certify_face <- function(values, candidate_mask, gene_ids, branch_ids,
                             max_edges) {
  n_edges <- sum(candidate_mask)
  if (n_edges > max_edges) {
    return(list(
      state = "NUMERICALLY_INDETERMINATE", retained_mask = candidate_mask & FALSE,
      forced_mask = candidate_mask & FALSE,
      certificate = list(
        method = "exact_residual_scc_facial_reduction_v1",
        status = "RESOURCE_CEILING_EXCEEDED", edge_count = n_edges,
        certificate_max_edges = max_edges, exact_or_interval_postcheck = FALSE
      )
    ))
  }
  positive_mask <- candidate_mask & values > 0
  pos <- .lt_components_mask(positive_mask, gene_ids, branch_ids)
  if (!all(pos$active_gene[rowSums(candidate_mask) > 0]) ||
      !all(pos$active_branch[colSums(candidate_mask) > 0])) {
    .lt_abort("Positive margins require every candidate vertex to have a positive edge.")
  }
  n_blocks <- length(pos$components)
  zero_linear <- which(candidate_mask & values == 0)
  nr <- nrow(values)
  zg <- ((zero_linear - 1L) %% nr) + 1L
  zb <- ((zero_linear - 1L) %/% nr) + 1L
  from <- pos$gene_component[zg]
  to <- pos$branch_component[zb]
  scc <- .lt_scc(n_blocks, from, to)
  retained_zero <- if (length(zero_linear)) scc[from] == scc[to] else logical()
  retained <- positive_mask
  retained[zero_linear[retained_zero]] <- TRUE
  forced <- candidate_mask & !retained
  if (any(values[forced] > 0)) {
    .lt_abort("Fatal certificate contradiction: positive Y is off the certified face.")
  }
  n_scc <- if (length(scc)) max(scc) else 0L
  state <- if (any(forced)) "FACIAL_BOUNDARY" else "FINITE_INTERIOR"
  if (any(forced)) {
    ranks <- .lt_condensation_rank(n_scc, scc[from], scc[to])
    fg <- ((which(forced) - 1L) %% nr) + 1L
    fb <- ((which(forced) - 1L) %/% nr) + 1L
    coeff <- ranks[scc[pos$branch_component[fb]]] -
      ranks[scc[pos$gene_component[fg]]]
    if (any(coeff <= 0L)) .lt_abort("Dual exposing rank failed strict sign verification.")
    denom <- sum(coeff)
    dual <- list(
      representation = "integer_condensation_rank_over_edge_coefficient_sum",
      denominator = denom,
      gene_rank_numerator = -ranks[scc[pos$gene_component]],
      branch_rank_numerator = ranks[scc[pos$branch_component]],
      forced_edge_coefficient_numerator = coeff,
      t_dot_h_exact = "0_by_positive_edge_zero_coefficient_and_zero_payload_elsewhere",
      all_candidate_edge_coefficients_nonnegative = TRUE,
      all_forced_edge_coefficients_strictly_positive = TRUE
    )
  } else {
    dual <- list(
      representation = "no_proper_exposer_residual_graph_single_scc_per_component",
      t_dot_h_exact = "not_applicable_no_proper_face",
      all_candidate_edge_coefficients_nonnegative = TRUE,
      all_forced_edge_coefficients_strictly_positive = TRUE
    )
  }
  final <- .lt_components_mask(retained, gene_ids, branch_ids)
  certificate <- list(
    method = "exact_residual_scc_facial_reduction_v1",
    theorem = paste(
      "binary64 payloads are exact nonnegative rational flows; residual SCC",
      "membership characterizes edges positive in some feasible flow; averaging",
      "cycle augmentations gives a strictly positive symbolic primal witness"
    ),
    original_nonnegative_primal_witness = "input_Y_on_candidate_edges",
    positive_support_component_count = n_blocks,
    residual_scc_count = n_scc,
    retained_edge_count = sum(retained),
    forced_zero_edge_count = sum(forced),
    positive_y_off_face_count = sum(values[forced] > 0),
    strict_positivity_certified = TRUE,
    minimal_face_certified = TRUE,
    primal_certificate_type = "EXACT_SYMBOLIC_RESIDUAL_CYCLE_CONVEX_WITNESS",
    dual_certificate_type = if (any(forced))
      "EXACT_INTEGER_CONDENSATION_EXPOSER" else "EXACT_NO_PROPER_FACE_SCC_CERTIFICATE",
    dual = dual,
    exact_or_interval_postcheck = TRUE,
    ambiguous_sign_count = 0L,
    retained_edge_sha256 = .lt_hash(which(retained)),
    forced_edge_sha256 = .lt_hash(which(forced)),
    final_component_sha256 = .lt_hash(final$summary)
  )
  certificate$certificate_sha256 <- .lt_hash(certificate)
  list(state = state, retained_mask = retained, forced_mask = forced,
       certificate = certificate, final_components = final)
}

#' Construct graph, boundary, and exact facial-support certificates
#'
#' @param domain An `lt_c3_domain`.
#' @param control A list returned by `lt_c3_control()`.
#' @return An internal-support object consumed by `lt_c3_fit()`.
#' @export
lt_c3_support <- function(domain, control = lt_c3_control()) {
  validate_lt_c3_domain(domain)
  values <- domain$matrix$values
  mask <- domain$fit_mask
  genes <- domain$matrix$gene_ids
  branches <- domain$matrix$branch_ids
  initial <- .lt_components_mask(mask, genes, branches)
  margins <- .lt_matrix_margins(values, mask)
  initial_margin_certificates <- .lt_component_margin_certificates(
    values, mask, initial, margins
  )
  if (!margins$certified ||
      (nrow(initial_margin_certificates) && !all(initial_margin_certificates$pass))) {
    return(structure(list(
      state = "NUMERICALLY_INDETERMINATE", domain = domain,
      initial_components = initial, margins = margins,
      initial_component_margin_certificates = initial_margin_certificates,
      reason = "deterministic_margin_certificate_failed"
    ), class = "lt_c3_support"))
  }
  unsupported_ids <- c(genes[!initial$active_gene], branches[!initial$active_branch])
  unsupported <- data.frame(
    vertex_type = c(rep("gene", sum(!initial$active_gene)),
                    rep("branch", sum(!initial$active_branch))),
    vertex_id = unsupported_ids,
    state = rep("UNSUPPORTED_VERTEX", length(unsupported_ids)),
    stringsAsFactors = FALSE
  )
  candidate <- mask
  edge_state <- rep(NA_character_, length(domain$edge_linear))
  edge_pos <- integer(length(mask)); edge_pos[domain$edge_linear] <- seq_along(domain$edge_linear)
  zero_total <- list(); zero_margin <- data.frame()
  for (i in seq_along(initial$components)) {
    z <- initial$components[[i]]
    cmask <- candidate & FALSE
    cmask[z$genes, z$branches] <- candidate[z$genes, z$branches, drop = FALSE]
    total <- .lt_pairwise_sum(values[which(cmask)])
    if (identical(total, 0)) {
      idx <- edge_pos[which(cmask)]
      edge_state[idx] <- "ZERO_TOTAL_COMPONENT"
      candidate[cmask] <- FALSE
      zero_total[[length(zero_total) + 1L]] <- data.frame(
        component_id = sprintf("component_%03d", i), edge_count = sum(cmask),
        state = "ZERO_TOTAL_COMPONENT", stringsAsFactors = FALSE
      )
    }
  }
  gene_degree <- rowSums(candidate)
  branch_degree <- colSums(candidate)
  zero_g <- gene_degree > 0 & margins$gene == 0
  zero_b <- branch_degree > 0 & margins$branch == 0
  boundary_mask <- candidate & (matrix(zero_g, nrow(values), ncol(values)) |
                                  matrix(zero_b, nrow(values), ncol(values), byrow = TRUE))
  if (any(boundary_mask)) {
    edge_state[edge_pos[which(boundary_mask)]] <- "ZERO_MARGIN_BOUNDARY"
    candidate[boundary_mask] <- FALSE
  }
  zero_margin <- rbind(
    data.frame(vertex_type = rep("gene", sum(zero_g)), vertex_id = genes[zero_g],
               state = rep("ZERO_MARGIN_BOUNDARY", sum(zero_g)), stringsAsFactors = FALSE),
    data.frame(vertex_type = rep("branch", sum(zero_b)), vertex_id = branches[zero_b],
               state = rep("ZERO_MARGIN_BOUNDARY", sum(zero_b)), stringsAsFactors = FALSE)
  )
  reduced <- .lt_components_mask(candidate, genes, branches)
  reduced_margin_certificates <- .lt_component_margin_certificates(
    values, candidate, reduced, margins
  )
  if (nrow(reduced_margin_certificates) && !all(reduced_margin_certificates$pass)) {
    return(structure(list(
      state = "NUMERICALLY_INDETERMINATE", domain = domain,
      initial_components = initial, reduced_components = reduced,
      margins = margins,
      initial_component_margin_certificates = initial_margin_certificates,
      reduced_component_margin_certificates = reduced_margin_certificates,
      reason = "reduced_component_margin_certificate_failed"
    ), class = "lt_c3_support"))
  }
  retained_global <- candidate & FALSE
  forced_global <- candidate & FALSE
  existence <- list()
  indeterminate <- FALSE
  for (i in seq_along(reduced$components)) {
    z <- reduced$components[[i]]
    cmask <- candidate & FALSE
    cmask[z$genes, z$branches] <- candidate[z$genes, z$branches, drop = FALSE]
    face <- .lt_certify_face(values, cmask, genes, branches,
                             control$certificate_max_edges)
    face$component_id <- sprintf("reduced_component_%03d", i)
    existence[[i]] <- face
    if (face$state == "NUMERICALLY_INDETERMINATE") {
      indeterminate <- TRUE
    } else {
      retained_global <- retained_global | face$retained_mask
      forced_global <- forced_global | face$forced_mask
      if (any(face$forced_mask))
        edge_state[edge_pos[which(face$forced_mask)]] <- "FACIAL_BOUNDARY"
    }
  }
  edge_state[edge_pos[which(retained_global)]] <- "FINITE_INTERIOR"
  final <- .lt_components_mask(retained_global, genes, branches)
  structure(list(
    state = if (indeterminate) "NUMERICALLY_INDETERMINATE" else
      if (any(forced_global)) "FACIAL_BOUNDARY" else
        if (any(retained_global)) "FINITE_INTERIOR" else
          if (length(zero_total)) "ZERO_TOTAL_COMPONENT" else "ZERO_MARGIN_BOUNDARY",
    domain = domain, control = control, margins = margins,
    initial_components = initial, reduced_components = reduced,
    initial_component_margin_certificates = initial_margin_certificates,
    reduced_component_margin_certificates = reduced_margin_certificates,
    final_components = final, unsupported_vertices = unsupported,
    zero_total_components = if (length(zero_total)) do.call(rbind, zero_total) else data.frame(),
    zero_margin_vertices = zero_margin,
    zero_margin_edge_mask = boundary_mask,
    facial_forced_edge_mask = forced_global,
    final_support_mask = retained_global,
    edge_support_state = edge_state,
    existence_certificates = existence
  ), class = "lt_c3_support")
}

.lt_margin_errors <- function(target, fitted, atol, rtol) {
  abs_error <- abs(fitted - target)
  rel_error <- ifelse(target > 0, abs_error / target, 0)
  list(
    pass = all(abs_error <= atol + rtol * target),
    max_absolute = if (length(abs_error)) max(abs_error) else 0,
    max_relative = if (length(rel_error)) max(rel_error) else 0,
    failing = which(abs_error > atol + rtol * target)
  )
}

.lt_objective_terms <- function(y, mu) {
  out <- numeric(length(y))
  zero <- y == 0
  out[zero] <- mu[zero]
  nz <- !zero
  d <- (y[nz] - mu[nz]) / mu[nz]
  small <- abs(d) < 1e-4
  term <- numeric(length(d))
  ds <- d[small]
  term[small] <- mu[nz][small] *
    (ds^2 / 2 - ds^3 / 6 + ds^4 / 12 - ds^5 / 20 + ds^6 / 30)
  term[!small] <- y[nz][!small] * log(y[nz][!small] / mu[nz][!small]) -
    y[nz][!small] + mu[nz][!small]
  bound <- 64 * .Machine$double.eps * pmax(1, y[nz], mu[nz])
  term[term < 0 & abs(term) <= bound] <- 0
  out[nz] <- term
  out
}

.lt_ipf_component <- function(values, mask, gene_index, branch_index, target_gene,
                              target_branch, control) {
  local_mask <- mask[gene_index, branch_index, drop = FALSE]
  R <- target_gene[gene_index]
  C <- target_branch[branch_index]
  n_edges <- sum(local_mask)
  T <- .lt_pairwise_sum(R)
  fit_atol <- control$fit_atol_scale * max(1, T / max(1, n_edges))
  G <- rep(1, length(gene_index)); B <- rep(1, length(branch_index))
  converged <- FALSE; denom_fail <- 0L; overflow <- 0L; underflow <- 0L
  for (iter in seq_len(control$max_iter)) {
    dg <- as.vector(local_mask %*% B)
    denom_fail <- denom_fail + sum(!is.finite(dg) | dg <= 0)
    if (denom_fail) break
    G <- R / dg
    db <- as.vector(crossprod(local_mask, G))
    denom_fail <- denom_fail + sum(!is.finite(db) | db <= 0)
    if (denom_fail) break
    B <- C / db
    overflow <- sum(!is.finite(G)) + sum(!is.finite(B))
    underflow <- sum(G == 0) + sum(B == 0)
    if (overflow || underflow) break
    ghat <- G * as.vector(local_mask %*% B)
    bhat <- B * as.vector(crossprod(local_mask, G))
    eg <- .lt_margin_errors(R, ghat, fit_atol, control$fit_rtol)
    eb <- .lt_margin_errors(C, bhat, fit_atol, control$fit_rtol)
    if (eg$pass && eb$pass) { converged <- TRUE; break }
  }
  if (!converged || denom_fail || overflow || underflow) {
    return(list(
      state = "NUMERICALLY_INDETERMINATE", convergence_flag = converged,
      iteration_count = iter, denominator_failure_count = denom_fail,
      overflow_count = overflow, underflow_to_zero_count = underflow
    ))
  }
  alpha <- log(G); beta <- log(B)
  shift <- .lt_pairwise_sum(alpha) / length(alpha)
  alpha <- alpha - shift; beta <- beta + shift
  coords <- which(local_mask, arr.ind = TRUE)
  y <- values[cbind(gene_index[coords[, 1L]], branch_index[coords[, 2L]])]
  mu <- exp(alpha[coords[, 1L]] + beta[coords[, 2L]])
  overflow <- sum(!is.finite(mu)); underflow <- sum(mu == 0)
  if (overflow || underflow) {
    return(list(
      state = "NUMERICALLY_INDETERMINATE", convergence_flag = TRUE,
      iteration_count = iter, denominator_failure_count = denom_fail,
      overflow_count = overflow, underflow_to_zero_count = underflow
    ))
  }
  pos <- seq_along(mu)
  gp <- split(pos, factor(coords[, 1L], levels = seq_along(gene_index)))
  bp <- split(pos, factor(coords[, 2L], levels = seq_along(branch_index)))
  ghat <- vapply(gp, function(ii) .lt_pairwise_sum(mu[ii]), numeric(1))
  bhat <- vapply(bp, function(ii) .lt_pairwise_sum(mu[ii]), numeric(1))
  eg <- .lt_margin_errors(R, ghat, fit_atol, control$fit_rtol)
  eb <- .lt_margin_errors(C, bhat, fit_atol, control$fit_rtol)
  objective_terms <- .lt_objective_terms(y, mu)
  objective <- .lt_pairwise_sum(objective_terms)
  gauge_residual <- abs(.lt_pairwise_sum(alpha))
  reconstructed <- exp(alpha[coords[, 1L]] + beta[coords[, 2L]])
  reconstruction_error <- max(abs(reconstructed - mu))
  objective_terms_nonnegative <- all(objective_terms >= 0)
  gauge_tolerance <- 64 * .Machine$double.eps * max(1, length(alpha))
  gauge_certified <- is.finite(gauge_residual) && gauge_residual <= gauge_tolerance
  state <- if (eg$pass && eb$pass && is.finite(objective) && objective >= 0 &&
               objective_terms_nonnegative && gauge_certified &&
               all(mu > 0) && all(is.finite(mu))) "FINITE_INTERIOR" else
    "NUMERICALLY_INDETERMINATE"
  list(
    state = state, convergence_flag = converged, iteration_count = iter,
    gene_index = gene_index, branch_index = branch_index,
    local_edge_row = coords[, 1L], local_edge_col = coords[, 2L],
    global_edge_linear = gene_index[coords[, 1L]] +
      (branch_index[coords[, 2L]] - 1L) * nrow(values),
    mu = mu, gene_log_factor = alpha, branch_log_factor = beta,
    gene_factor = exp(alpha), branch_factor = exp(beta),
    objective = objective, objective_terms_nonnegative = objective_terms_nonnegative,
    fit_atol = fit_atol, fit_rtol = control$fit_rtol,
    gene_errors = eg, branch_errors = eb,
    denominator_failure_count = denom_fail, overflow_count = overflow,
    underflow_to_zero_count = underflow,
    positive_y_with_zero_mu_count = sum(y > 0 & mu == 0),
    gauge_shift = shift, gauge_residual = gauge_residual,
    gauge_tolerance = gauge_tolerance, gauge_certified = gauge_certified,
    reconstruction_error = reconstruction_error,
    mu_sha256 = .lt_hash(mu), gene_factor_sha256 = .lt_hash(exp(alpha)),
    branch_factor_sha256 = .lt_hash(exp(beta))
  )
}

.lt_certificate_required <- c(
  "work_order", "run_id", "component_id", "parent_component_id",
  "input_values_sha256", "coordinate_state_sha256",
  "measurement_eligibility_sha256", "ordered_gene_ledger_sha256",
  "ordered_branch_ledger_sha256", "ordered_cell_ledger_sha256",
  "configuration_sha256", "software_environment_sha256",
  "initial_gene_count", "initial_branch_count", "initial_edge_count",
  "unsupported_vertex_ledger_sha256", "zero_margin_vertex_ledger_sha256",
  "forced_zero_edge_ledger_sha256", "facial_edge_ledger_sha256",
  "final_component_membership_sha256", "gene_margin_sha256",
  "branch_margin_sha256", "component_total_hex", "summation_algorithm",
  "canonical_accumulation_order", "sum_atol_hex",
  "gene_total_residual_hex", "branch_total_residual_hex",
  "marginal_cone_method_version", "primal_certificate_type",
  "primal_certificate_sha256", "dual_certificate_type",
  "dual_certificate_sha256", "strict_positivity_certified",
  "minimal_face_certified", "exact_or_interval_postcheck",
  "ambiguous_sign_count", "solver_name_and_version",
  "deterministic_settings_sha256", "convergence_flag", "iteration_count",
  "mu_sha256", "mu_all_finite", "mu_strictly_positive_on_certified_support",
  "mu_exact_zero_off_certified_face", "positive_y_with_zero_mu_count",
  "objective_name", "objective_hex", "objective_finite_nonnegative",
  "objective_terms_nonnegative", "gene_margin_agreement_pass",
  "branch_margin_agreement_pass",
  "maximum_gene_margin_absolute_error_hex",
  "maximum_gene_margin_relative_error_hex",
  "maximum_branch_margin_absolute_error_hex",
  "maximum_branch_margin_relative_error_hex", "fit_atol_hex", "fit_rtol_hex",
  "denominator_failure_count", "overflow_count", "underflow_to_zero_count",
  "gauge_rule", "gauge_shift_hex", "gene_factor_sha256",
  "branch_factor_sha256", "gauge_residual_hex", "gauge_tolerance_hex",
  "gauge_certified", "reconstruction_error_hex", "reconstruction_certified",
  "replay_mode", "replay_mu_sha256", "replay_certificate_sha256",
  "deterministic_replay_pass", "commands_sha256"
)

#' Validate a machine-readable C3 fit certificate
#'
#' @param certificate A certificate emitted by `lt_c3_fit()`.
#' @return The certificate, invisibly, or an error.
#' @export
validate_lt_c3_certificate <- function(certificate) {
  if (!is.list(certificate) || is.null(names(certificate))) {
    .lt_abort("C3 certificate must be a named list.")
  }
  missing <- setdiff(.lt_certificate_required, names(certificate))
  if (length(missing)) {
    .lt_abort(sprintf("C3 certificate missing required fields: %s.",
                      paste(missing, collapse = ", ")))
  }
  if (!isTRUE(certificate$strict_positivity_certified) ||
      !isTRUE(certificate$minimal_face_certified) ||
      !isTRUE(certificate$exact_or_interval_postcheck) ||
      !isTRUE(certificate$mu_all_finite) ||
      !isTRUE(certificate$mu_strictly_positive_on_certified_support) ||
      !isTRUE(certificate$mu_exact_zero_off_certified_face) ||
      !isTRUE(certificate$objective_finite_nonnegative) ||
      !isTRUE(certificate$objective_terms_nonnegative) ||
      !isTRUE(certificate$gene_margin_agreement_pass) ||
      !isTRUE(certificate$branch_margin_agreement_pass) ||
      !isTRUE(certificate$convergence_flag) ||
      !isTRUE(certificate$gauge_certified) ||
      !isTRUE(certificate$reconstruction_certified) ||
      !isTRUE(certificate$deterministic_replay_pass) ||
      certificate$ambiguous_sign_count != 0L ||
      certificate$positive_y_with_zero_mu_count != 0L ||
      certificate$denominator_failure_count != 0L ||
      certificate$overflow_count != 0L || certificate$underflow_to_zero_count != 0L) {
    .lt_abort("C3 certificate contains a failed production gate.")
  }
  invisible(certificate)
}

#' Fit the certified C3 baseline working mean
#'
#' Fits only certified support edges and returns no off-domain products.
#' The normative target is a generalized-KL/Bregman margin working mean, not a
#' literal Poisson probability model.
#'
#' @param domain An `lt_c3_domain`.
#' @param control A list from `lt_c3_control()`.
#' @param run_id Stable run identifier.
#' @return An `lt_c3_fit` with edge-aligned supported means and certificates.
#' @export
lt_c3_fit <- function(domain, control = lt_c3_control(), run_id = "c3_run") {
  validate_lt_c3_domain(domain)
  .lt_assert_scalar_character(run_id, "run_id")
  support <- lt_c3_support(domain, control)
  if (support$state == "NUMERICALLY_INDETERMINATE") {
    return(structure(list(
      run_id = run_id, state = "NUMERICALLY_INDETERMINATE",
      domain = domain, support = support, mu = NULL, certificates = list(),
      scientific_pass = FALSE
    ), class = "lt_c3_fit"))
  }
  values <- domain$matrix$values
  final <- support$final_components
  mu <- rep(NA_real_, length(domain$edge_linear))
  boundary <- support$edge_support_state %in%
    c("ZERO_TOTAL_COMPONENT", "ZERO_MARGIN_BOUNDARY", "FACIAL_BOUNDARY")
  mu[boundary] <- 0
  edge_pos <- integer(length(values)); edge_pos[domain$edge_linear] <- seq_along(domain$edge_linear)
  fits <- list(); certs <- list(); held <- FALSE
  config_hash <- .lt_hash(control)
  env <- list(R = R.version.string, platform = R.version$platform,
              BLAS = extSoftVersion()[["BLAS"]], LaTerra = "0.0.0.9002")
  env_hash <- .lt_hash(env)
  for (i in seq_along(final$components)) {
    z <- final$components[[i]]
    fit <- .lt_ipf_component(
      values, support$final_support_mask, z$genes, z$branches,
      support$margins$gene, support$margins$branch, control
    )
    if (fit$state != "FINITE_INTERIOR") { held <- TRUE; fits[[i]] <- fit; next }
    replay <- if (control$replay) .lt_ipf_component(
      values, support$final_support_mask, z$genes, z$branches,
      support$margins$gene, support$margins$branch, control
    ) else fit
    replay_pass <- identical(fit$mu_sha256, replay$mu_sha256) &&
      identical(fit$iteration_count, replay$iteration_count) &&
      identical(fit$state, replay$state)
    idx <- edge_pos[fit$global_edge_linear]
    mu[idx] <- fit$mu
    parent_face <- Filter(function(x) any(x$retained_mask[fit$global_edge_linear]),
                          support$existence_certificates)
    face_cert <- parent_face[[1L]]$certificate
    comp_mask <- support$final_support_mask[z$genes, z$branches, drop = FALSE]
    comp_edges <- sum(comp_mask)
    comp_local_linear <- which(comp_mask)
    comp_local_row <- ((comp_local_linear - 1L) %% length(z$genes)) + 1L
    comp_local_col <- ((comp_local_linear - 1L) %/% length(z$genes)) + 1L
    comp_global_linear <- z$genes[comp_local_row] +
      (z$branches[comp_local_col] - 1L) * nrow(values)
    comp_total <- .lt_pairwise_sum(values[comp_global_linear])
    comp_sum_atol <- 64 * .Machine$double.eps * max(1, abs(comp_total)) *
      max(1, ceiling(log2(max(1, comp_edges))))
    comp_gene_total_residual <- abs(
      .lt_pairwise_sum(support$margins$gene[z$genes]) - comp_total
    )
    comp_branch_total_residual <- abs(
      .lt_pairwise_sum(support$margins$branch[z$branches]) - comp_total
    )
    cert <- list(
      certificate_status = if (support$state == "FACIAL_BOUNDARY")
        "PASS_FACIAL_BOUNDARY" else "PASS_FINITE_INTERIOR",
      work_order = "LT000A-T1A", run_id = run_id,
      component_id = sprintf("final_component_%03d", i),
      parent_component_id = parent_face[[1L]]$component_id,
      input_values_sha256 = domain$hashes$input_values_sha256,
      coordinate_state_sha256 = domain$hashes$coordinate_state_sha256,
      measurement_eligibility_sha256 = domain$hashes$measurement_eligibility_sha256,
      ordered_gene_ledger_sha256 = domain$hashes$ordered_gene_ledger_sha256,
      ordered_branch_ledger_sha256 = domain$hashes$ordered_branch_ledger_sha256,
      ordered_cell_ledger_sha256 = domain$hashes$ordered_cell_ledger_sha256,
      configuration_sha256 = config_hash, software_environment_sha256 = env_hash,
      initial_gene_count = nrow(values), initial_branch_count = ncol(values),
      initial_edge_count = length(domain$edge_linear),
      unsupported_vertex_ledger_sha256 = .lt_hash(support$unsupported_vertices),
      zero_margin_vertex_ledger_sha256 = .lt_hash(support$zero_margin_vertices),
      forced_zero_edge_ledger_sha256 = .lt_hash(which(support$zero_margin_edge_mask)),
      facial_edge_ledger_sha256 = .lt_hash(which(support$facial_forced_edge_mask)),
      final_component_membership_sha256 = final$summary$membership_sha256[[i]],
      gene_margin_sha256 = .lt_hash(support$margins$gene[z$genes]),
      branch_margin_sha256 = .lt_hash(support$margins$branch[z$branches]),
      component_total_hex = .lt_hex(comp_total),
      summation_algorithm = control$summation_algorithm,
      canonical_accumulation_order = control$canonical_accumulation_order,
      sum_atol_hex = .lt_hex(comp_sum_atol),
      gene_total_residual_hex = .lt_hex(comp_gene_total_residual),
      branch_total_residual_hex = .lt_hex(comp_branch_total_residual),
      marginal_cone_method_version = control$existence_algorithm,
      primal_certificate_type = face_cert$primal_certificate_type,
      primal_certificate_sha256 = face_cert$certificate_sha256,
      dual_certificate_type = face_cert$dual_certificate_type,
      dual_certificate_sha256 = .lt_hash(face_cert$dual),
      strict_positivity_certified = face_cert$strict_positivity_certified,
      minimal_face_certified = face_cert$minimal_face_certified,
      exact_or_interval_postcheck = face_cert$exact_or_interval_postcheck,
      ambiguous_sign_count = face_cert$ambiguous_sign_count,
      solver_name_and_version = "deterministic_bipartite_IPF_v1",
      deterministic_settings_sha256 = config_hash,
      convergence_flag = fit$convergence_flag, iteration_count = fit$iteration_count,
      mu_sha256 = fit$mu_sha256, mu_all_finite = all(is.finite(fit$mu)),
      mu_strictly_positive_on_certified_support = all(fit$mu > 0),
      mu_exact_zero_off_certified_face = all(mu[boundary][!is.na(mu[boundary])] == 0),
      positive_y_with_zero_mu_count = fit$positive_y_with_zero_mu_count,
      objective_name = "generalized_KL_Bregman_working_mean",
      objective_hex = .lt_hex(fit$objective),
      objective_finite_nonnegative = is.finite(fit$objective) && fit$objective >= 0,
      objective_terms_nonnegative = fit$objective_terms_nonnegative,
      gene_margin_agreement_pass = fit$gene_errors$pass,
      branch_margin_agreement_pass = fit$branch_errors$pass,
      maximum_gene_margin_absolute_error_hex = .lt_hex(fit$gene_errors$max_absolute),
      maximum_gene_margin_relative_error_hex = .lt_hex(fit$gene_errors$max_relative),
      maximum_branch_margin_absolute_error_hex = .lt_hex(fit$branch_errors$max_absolute),
      maximum_branch_margin_relative_error_hex = .lt_hex(fit$branch_errors$max_relative),
      fit_atol_hex = .lt_hex(fit$fit_atol), fit_rtol_hex = .lt_hex(fit$fit_rtol),
      denominator_failure_count = fit$denominator_failure_count,
      overflow_count = fit$overflow_count,
      underflow_to_zero_count = fit$underflow_to_zero_count,
      gauge_rule = "sum_gene_log_factors_zero_per_positive_component",
      gauge_shift_hex = .lt_hex(fit$gauge_shift),
      gene_factor_sha256 = fit$gene_factor_sha256,
      branch_factor_sha256 = fit$branch_factor_sha256,
      gauge_residual_hex = .lt_hex(fit$gauge_residual),
      gauge_tolerance_hex = .lt_hex(fit$gauge_tolerance),
      gauge_certified = fit$gauge_certified,
      reconstruction_error_hex = .lt_hex(fit$reconstruction_error),
      reconstruction_certified = is.finite(fit$reconstruction_error) &&
        fit$reconstruction_error == 0,
      replay_mode = if (control$replay) "same_environment_full_solver_replay" else "disabled",
      replay_mu_sha256 = replay$mu_sha256,
      replay_certificate_sha256 = .lt_hash(list(
        replay$state, replay$iteration_count, replay$mu_sha256,
        replay$gene_factor_sha256, replay$branch_factor_sha256
      )),
      deterministic_replay_pass = replay_pass,
      commands_sha256 = .lt_hash("lt_c3_fit(domain, control, run_id)"),
      prohibitions = c("no_off_domain_materialization", "no_cross_component_products",
                       "no_poisson_inference_claim", "no_trait_dependent_fit_domain",
                       "no_pseudocount")
    )
    validate_lt_c3_certificate(cert)
    fits[[i]] <- fit; certs[[i]] <- cert
  }
  if (any(is.na(mu[support$edge_support_state != "UNSUPPORTED_VERTEX"]))) held <- TRUE
  structure(list(
    run_id = run_id,
    state = if (held) "NUMERICALLY_INDETERMINATE" else support$state,
    domain = domain, support = support,
    edge_linear = domain$edge_linear,
    edge_gene_index = domain$edge_gene_index,
    edge_branch_index = domain$edge_branch_index,
    observed_value = domain$edge_value,
    mu = if (held) NULL else mu,
    baseline_support_state = support$edge_support_state,
    component_fits = fits, certificates = certs,
    scientific_pass = FALSE,
    interpretation = "software evidence only; independent scientific adjudication required"
  ), class = "lt_c3_fit")
}

#' Validate a C3 fit object
#'
#' @param x An `lt_c3_fit`.
#' @return `x`, invisibly.
#' @export
validate_lt_c3_fit <- function(x) {
  if (!inherits(x, "lt_c3_fit")) .lt_abort("`x` must inherit from lt_c3_fit.")
  validate_lt_c3_domain(x$domain)
  if (!x$state %in% lt_c3_existence_states()) .lt_abort("Unknown C3 fit state.")
  if (identical(x$state, "NUMERICALLY_INDETERMINATE")) {
    if (!is.null(x$mu)) .lt_abort("Held C3 fits must not release partial production mu.")
    return(invisible(x))
  }
  if (!is.numeric(x$mu) || length(x$mu) != length(x$domain$edge_linear)) {
    .lt_abort("C3 fitted means must align exactly to the admitted edge ledger.")
  }
  if (any(!is.finite(x$mu)) || any(x$mu < 0)) {
    .lt_abort("Released C3 fitted means must be finite and nonnegative.")
  }
  if (any(x$observed_value > 0 & x$mu == 0)) {
    .lt_abort("A positive observed payload cannot have fitted mean zero.")
  }
  if (!identical(x$edge_linear, x$domain$edge_linear)) {
    .lt_abort("C3 fit edge order differs from the declared domain.")
  }
  for (certificate in x$certificates) validate_lt_c3_certificate(certificate)
  invisible(x)
}

#' @export
print.lt_c3_domain <- function(x, ...) {
  validate_lt_c3_domain(x)
  cat("<lt_c3_domain>\n")
  cat("  matrix: ", nrow(x$matrix$values), " genes x ", ncol(x$matrix$values), " branches\n", sep = "")
  cat("  admitted edges: ", length(x$edge_linear), "\n", sep = "")
  cat("  exact zeros: ", sum(x$edge_value == 0), "\n", sep = "")
  cat("  upstream QC certified: ", isTRUE(x$eligibility$provenance$upstream_qc_certified), "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_c3_fit <- function(x, ...) {
  cat("<lt_c3_fit> ", x$run_id, "\n", sep = "")
  cat("  state: ", x$state, "\n", sep = "")
  cat("  fit-domain edges: ", length(x$domain$edge_linear), "\n", sep = "")
  cat("  fitted means released: ", if (is.null(x$mu)) 0 else length(x$mu), "\n", sep = "")
  cat("  scientific pass: FALSE (independent scientific adjudication required)\n")
  invisible(x)
}
