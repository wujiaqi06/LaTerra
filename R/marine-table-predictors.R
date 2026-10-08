.lt_marine_reporting_clarification <- function() {
  list(id = "TARGET001-REPORTING-CLARIFICATION001",
    sha256 = "1324bdc74bbd57c8b7a579d4585ba32daca85f9a7621b84e0b3144b63de6b1df",
    qualifier = "presentation_rule_frozen_by_main_control_original_generator_unlocated",
    fields = c("OUT_TABLE_S6/Fig5B_perm_summary!F2", "OUT_TABLE_S5/predictor_annotation!H2:H149"))
}

.lt_marine_table_percent <- function(proportion) {
  if (!is.numeric(proportion) || length(proportion) != 1L ||
      !is.finite(proportion) || proportion < 0 || proportion > 1)
    .lt_marine_abort("S6 observed percent requires one finite proportion in [0,1]; no stored-answer fallback.",
      "reporting_payload")
  if (Sys.getlocale("LC_NUMERIC") != "C")
    .lt_marine_abort("Set LC_NUMERIC=C for the frozen reporting formatter.", "reporting_locale")
  sprintf("%.1f%%", 100 * proportion)
}

.lt_marine_table_direction <- function(marine, aquatic, marine_selected, aquatic_selected) {
  if (!is.logical(marine_selected) || !is.logical(aquatic_selected) ||
      length(marine_selected) != 1L || length(aquatic_selected) != 1L ||
      is.na(marine_selected) || is.na(aquatic_selected) || !(marine_selected || aquatic_selected))
    .lt_marine_abort("S5 direction requires typed current selected-union membership.", "reporting_payload")
  phrase <- function(value, selected, axis) {
    if (!selected) return(character())
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value == 0)
      .lt_marine_abort("Selected S5 coefficient must be finite and nonzero; no epsilon or forced sign.",
        "reporting_payload")
    paste(axis, if (value > 0) "fast direction positive" else "slow direction negative")
  }
  paste(c(phrase(marine, marine_selected, "marine"),
    phrase(aquatic, aquatic_selected, "aquatic")), collapse = "; ")
}

.lt_marine_table_bool <- function(x, field) {
  if (is.logical(x) && !anyNA(x)) return(x)
  if (!is.character(x) || anyNA(x) || any(!x %in% c("True", "False")))
    .lt_marine_abort(paste0("Unrecognized typed display flag: ", field), "reporting_payload")
  x == "True"
}

.lt_marine_table_s5 <- function(full, annotation) {
  ids <- c("fix_marine_binary", "fix_aquatic_v2")
  beta <- lapply(ids, function(id) {
    model <- full$models[[id]]
    b <- model$beta
    if (is.null(model$fit) || !is.numeric(b) || !length(b) || any(!is.finite(b)))
      .lt_marine_abort("S5 requires the actual finite baseline fitted coefficients.", "reporting_dependency")
    .lt_assert_unique_ids(names(b), "S5 fitted coefficient gene keys")
    fit_beta <- as.matrix(model$fit$beta)[, 1L]
    if (!identical(names(b), names(fit_beta)) || !identical(as.numeric(b), as.numeric(fit_beta)))
      .lt_marine_abort("S5 coefficient vector differs from its fitted source support.", "reporting_dependency")
    b
  })
  selected <- lapply(beta, function(b) names(b)[b != 0])
  genes <- union(selected[[1L]], selected[[2L]])
  .lt_assert_unique_ids(annotation$gene, "S5 approved annotation gene keys")
  if (!length(genes) || !setequal(genes, annotation$gene))
    .lt_marine_abort("Fresh baseline selected union and approved annotation keys differ.", "reporting_keys")
  lookup <- match(genes, annotation$gene)
  headers <- unlist(.lt_marine_reporting_layout("OUT_TABLE_S5", "predictor_annotation")$schema$headers,
    use.names = FALSE)
  admitted <- headers[c(12:20, 24:33)]
  if (!all(admitted %in% names(annotation)))
    .lt_marine_abort("Missing approved S5 annotation/display fields.", "reporting_fields")
  out <- data.frame(gene = genes, stringsAsFactors = FALSE)
  out <- cbind(out, annotation[lookup, admitted, drop = FALSE]); rownames(out) <- NULL
  out$selected_in_marine <- genes %in% selected[[1L]]
  out$selected_in_aquatic <- genes %in% selected[[2L]]
  out$marine_coef <- out$aquatic_coef <- numeric(length(genes))
  out$marine_coef[out$selected_in_marine] <- beta[[1L]][genes[out$selected_in_marine]]
  out$aquatic_coef[out$selected_in_aquatic] <- beta[[2L]][genes[out$selected_in_aquatic]]
  out$lasso_group <- ifelse(out$selected_in_marine & out$selected_in_aquatic, "shared",
    ifelse(out$selected_in_marine, "marine_specific", "aquatic_specific"))
  out$max_abs_coef <- pmax(abs(out$marine_coef), abs(out$aquatic_coef))
  out$coefficient_direction_summary <- vapply(seq_along(genes), function(i)
    .lt_marine_table_direction(out$marine_coef[i], out$aquatic_coef[i],
      out$selected_in_marine[i], out$selected_in_aquatic[i]), "")
  out$coef_marker_marine_gt_0_1 <- ifelse(abs(out$marine_coef) > .1, "#", NA_character_)
  out$coef_marker_aquatic_gt_0_1 <- ifelse(abs(out$aquatic_coef) > .1, "*", NA_character_)
  color_coef <- ifelse(out$lasso_group == "shared", (out$marine_coef + out$aquatic_coef) / 2,
    ifelse(out$selected_in_marine, out$marine_coef, out$aquatic_coef))
  out$label_color_direction <- ifelse(color_coef >= 0, "fast_positive_red", "slow_negative_blue")
  flags <- c("counted_in_Figure4C_circle_size", "displayed_in_approved_representative_gene_table",
    "recommended_for_Figure4C_original", "recommended_for_Figure4C_display", "recommended_for_TableS5_only",
    "keep_unassigned", "reported_in_Supplementary_Table_S5")
  for (field in flags) out[[field]] <- .lt_marine_table_bool(out[[field]], field)
  counted <- out$counted_in_Figure4C_circle_size
  rep <- out$displayed_in_approved_representative_gene_table
  if (any(rep & !counted) || anyNA(out$display_module) || any(!nzchar(out$display_module)))
    .lt_marine_abort("Approved representative/module metadata is inconsistent.", "reporting_payload")
  rank_map <- c(high = 0L, medium = 1L, `medium-low` = 2L, low = 2L, unassigned = 3L)
  confidence <- unname(rank_map[tolower(trimws(out$annotation_confidence))])
  confidence[is.na(confidence)] <- 99L
  marker <- abs(out$marine_coef) > .1 | abs(out$aquatic_coef) > .1
  cell_key <- paste(out$display_module, out$lasso_group, sep = "\r")
  out$Figure4C_cell_total_genes <- NA_integer_
  out$Figure4C_cell_representative_genes <- out$Figure4C_cell_all_genes <- NA_character_
  for (cell in unique(cell_key[counted])) {
    rows <- which(counted & cell_key == cell)
    prefix <- data.frame(confidence = confidence[rows], max_abs_coef = out$max_abs_coef[rows], marker = marker[rows])
    if (anyDuplicated(prefix))
      .lt_marine_abort("Fresh S5 cell sort prefix has a tie; unadmitted priority is required. List assembly stopped.",
        "reporting_priority_tie")
    ordered <- rows[order(prefix$confidence, -prefix$max_abs_coef, -as.integer(prefix$marker), method = "radix")]
    reps <- rows[rep[rows]]
    approved_order <- out$approved_representative_display_order_within_cell[reps]
    if (!is.numeric(approved_order) || anyNA(approved_order) || any(!is.finite(approved_order)) ||
        anyDuplicated(approved_order) || !setequal(approved_order, seq_along(reps)))
      .lt_marine_abort("Missing/ambiguous approved representative order; no re-selection.", "reporting_keys")
    reps <- reps[order(approved_order)]
    out$Figure4C_cell_total_genes[rows] <- length(rows)
    out$Figure4C_cell_all_genes[rows] <- paste(out$gene[ordered], collapse = ";")
    out$Figure4C_cell_representative_genes[rows] <- paste(out$gene[reps], collapse = ";")
  }
  group_order <- match(out$lasso_group, c("shared", "marine_specific", "aquatic_specific"))
  out <- out[order(group_order, out$display_module, out$gene, method = "radix"), headers, drop = FALSE]
  modules <- unique(out$display_module)
  summary <- lapply(modules, function(module) {
    m <- out[out$display_module == module, , drop = FALSE]
    counted <- m$counted_in_Figure4C_circle_size
    rep <- m$displayed_in_approved_representative_gene_table
    only <- m$recommended_for_TableS5_only
    data.frame(display_module = module, n_total_predictors_in_TableS5 = nrow(m),
      n_counted_in_Figure4C_circle_size = sum(counted),
      n_displayed_representative_labels_in_approved_table = sum(rep), n_TableS5_only = sum(only),
      n_unassigned = sum(m$keep_unassigned), n_shared = sum(m$lasso_group == "shared"),
      n_marine_specific = sum(m$lasso_group == "marine_specific"),
      n_aquatic_specific = sum(m$lasso_group == "aquatic_specific"),
      all_predictors = paste(m$gene, collapse = ";"),
      predictors_counted_in_Figure4C = paste(m$gene[counted], collapse = ";"),
      representative_labels = paste(m$gene[rep], collapse = ";"),
      TableS5_only_predictors = paste(m$gene[only], collapse = ";"), stringsAsFactors = FALSE)
  })
  list(predictor_annotation = out, display_module_summary = do.call(rbind, summary),
    provenance = list(clarification = .lt_marine_reporting_clarification(),
      selected_membership_from_fitted_source = TRUE, priority_input_read = FALSE,
      fresh_prefix_unique_in_each_counted_cell = TRUE, representative_selection_recomputed = FALSE,
      admitted_input_fields = admitted, original_numeric_annotation_columns_consumed = FALSE))
}
