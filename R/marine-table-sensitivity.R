.lt_marine_reporting_fields <- function(source_id, sheet, current) {
  x <- .lt_marine_reporting_layout(source_id, sheet)
  dynamic <- c("EXISTING_ADMITTED_INPUT_JOIN", "RECOMPUTE_FROM_THIS_RUN", "THIS_RUN_KEY_OR_MANIFEST_JOIN")
  cols <- unlist(x$schema$headers, use.names = FALSE)[
    colSums(matrix(x$roles %in% dynamic, nrow(x$roles))) > 0L]
  required <- union(cols, unlist(x$schema$key_fields, use.names = FALSE))
  if (!is.data.frame(current) || anyDuplicated(names(current)) || !all(required %in% names(current)))
    .lt_marine_abort("Missing derived reporting fields; no archived output fallback.", "reporting_fields")
  current[required]
}

.lt_marine_table_s4_sensitivity <- function(models) {
  ids <- c("fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear", "fix_drop_sea_otter",
    "fix_drop_pinniped", "fix_whale_only", "fix_pinniped_only", "fix_aquatic_v2",
    "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia", "fix_aquatic_v2_noMarineEdge",
    "fix_aquatic_v2_noNonMarineCarnivores", "fix_aquatic_v2_noNonMarineRodents", "fix_aquatic_v2_noHippoLechwe")
  if (!identical(names(models), ids))
    .lt_marine_abort("S4 requires all fourteen current coherent global screens, not Phase13 exceptions.", "reporting_keys")
  # Same source-traced run display labels already used in the full-data output;
  # original Stage2 QC joins drop_taxa_or_group, May22 replaces baseline by baseline.
  drop <- c("baseline", "Cetacea", "Ursus_maritimus", "Enhydra_lutris_kenyoni", "Pinnipedia",
    "Cetacea only", "Pinnipedia only", "baseline", "Cetacea", "Pinnipedia", "MarineEdge",
    "NonMarineCarnivores", "NonMarineRodents", "HippoLechwe")
  rows <- lapply(seq_along(models), function(j) {
    m <- models[[j]]; s <- m$global_screen
    if (m$run_id != ids[j] || nrow(s$summary) != 1L ||
        nrow(s$tested) + sum(!s$gene_ledger$tested) != nrow(s$gene_ledger))
      .lt_marine_abort("Incomplete global screening ledger for S4.", "reporting_dependency")
    trait <- if (j <= 7L) "marine_binary" else "aquatic_v2"
    baseline <- j %in% c(1L, 8L)
    column <- if (baseline) trait else m$trait_column
    n_sig <- nrow(s$significant)
    slow <- sum(s$significant$tvalue > 0); fast <- sum(s$significant$tvalue < 0)
    data.frame(run_id = ids[j], trait = trait, trait_column = column,
      analysis_type = if (baseline) "baseline_ttest" else "drop_sensitivity_ttest", drop_rule = drop[j],
      n_input_genes = nrow(s$gene_ledger), n_tested_genes = nrow(s$tested),
      n_skipped_genes = sum(!s$gene_ledger$tested), n_branches_screened = sum(m$state$screened),
      n_significant_FDR_0_01 = n_sig, n_slow = slow, n_fast = fast,
      slow_proportion = if (n_sig) slow / n_sig else NA_real_,
      fast_proportion = if (n_sig) fast / n_sig else NA_real_, stringsAsFactors = FALSE)
  })
  list(marine_sensitivity = do.call(rbind, rows[1:7]), aquatic_sensitivity = do.call(rbind, rows[8:14]))
}

.lt_marine_phase13_report_decision <- function(auc, prop_no_selected) {
  if (length(auc) != 1L || length(prop_no_selected) != 1L || !is.finite(auc) ||
      !is.finite(prop_no_selected) || prop_no_selected < 0 || prop_no_selected > 1)
    .lt_marine_abort("Phase13 report decision requires current finite AUC and evaluated-fold zero-selection fraction.",
      "reporting_payload")
  if (auc >= .48 && auc <= .52 && prop_no_selected >= .95) "strict_null_collapse" else
    if (auc < .60 || prop_no_selected >= .80) "near_collapse" else
      if (auc < .75) "strong_attenuation" else "retained_sparse_prediction"
}

.lt_marine_table_boolean_text <- function(x) {
  if (!is.logical(x) || anyNA(x))
    .lt_marine_abort("S6 defined-status display requires complete logical values from this run.", "reporting_payload")
  # Original final S6 TSV parser/maybeNumber preserves True/False as text.
  ifelse(x, "True", "False")
}

.lt_marine_table_s6 <- function(sensitivity, fig5b, worlds, turnover) {
  a <- .lt_marine_sensitivity_tables(sensitivity)
  fold <- lapply(seq_along(sensitivity), function(j) {
    result <- sensitivity[[j]]
    evaluated <- vapply(result$folds, function(f) length(f$identity$test_species) > 0L, TRUE)
    n_selected <- vapply(result$folds[evaluated], function(f) length(f$beta), 1L)
    row <- a$table[j, ]
    row$nested_AUC_raw <- a$unrounded$AUC_unrounded[j]
    row$median_selected_predictors_per_fold <- stats::median(n_selected)
    row$IQR_selected_predictors_per_fold <- .lt_s2_iqr_text(n_selected)
    row$decision_label <- .lt_marine_phase13_report_decision(row$nested_AUC_raw, mean(n_selected == 0))
    row
  })
  # Terminal label counts come from the same admitted persisted world ledger,
  # not the post-ASR branch counts in fig5b$screening.
  .lt_assert_unique_ids(paste(worlds$perm_id, worlds$species_id, sep = "::"), "Fig5B reporting world/taxon keys")
  perm <- fig5b$screening
  if (!identical(perm$perm_id, unique(worlds$perm_id)) || any(!worlds$perm_state %in% c(0, 1)))
    .lt_marine_abort("S6 permutation/world ledgers differ.", "reporting_keys")
  counts <- lapply(perm$perm_id, function(id) {
    w <- worlds[worlds$perm_id == id, , drop = FALSE]
    if (length(unique(w$seed)) != 1L || anyNA(w))
      .lt_marine_abort("Missing/ambiguous persisted permutation identity.", "reporting_keys")
    c(positive = sum(w$perm_state == 1), negative = sum(w$perm_state == 0), seed = unique(w$seed))
  })
  counts <- do.call(rbind, counts)
  if (any(perm$seed != counts[, "seed"]) || length(unique(counts[, "positive"])) != 1L ||
      length(unique(counts[, "negative"])) != 1L)
    .lt_marine_abort("S6 cannot silently reconcile inconsistent persisted terminal counts/seeds.", "reporting_keys")
  perm$matched_terminal_positive_count <- counts[, "positive"]
  perm$matched_terminal_negative_count <- counts[, "negative"]
  perm$slow_proportion_defined <- .lt_marine_table_boolean_text(perm$slow_proportion_defined)
  p <- fig5b$observed_vs_null
  if (nrow(p) != 1L || p$n_permutations != nrow(perm))
    .lt_marine_abort("S6 observed/null summary is not the complete current ensemble.", "reporting_dependency")
  p$observed_slow_percent <- .lt_marine_table_percent(p$observed_slow_proportion)
  p$matched_terminal_positive_count <- counts[1, "positive"]
  p$matched_terminal_negative_count <- counts[1, "negative"]
  main <- c("marine baseline vs whale-only", "binary aquatic baseline vs aquatic no Cetacea")
  cross <- "marine baseline vs drop-cetaceans slow genes"
  if (!setequal(turnover$observed$comparison, c(main, cross)) || anyDuplicated(turnover$observed$comparison))
    .lt_marine_abort("Missing/ambiguous main versus cross-layer turnover comparisons.", "reporting_keys")
  cross_data <- turnover$summary[turnover$summary$comparison == cross, , drop = FALSE]
  origin <- turnover$observed[turnover$observed$comparison == cross, , drop = FALSE]
  cross_data$set_A_name <- origin$set_A_name
  cross_data$set_B_name <- origin$set_B_name
  list(Fig5A_sensitivity = a$table, Fig5A_fold_summary = do.call(rbind, fold),
    Fig5B_perm_summary = p, Fig5B_perm_distribution = perm,
    Fig5C_turnover_metrics = turnover$observed[turnover$observed$comparison %in% main, ],
    Fig5C_null_summary = turnover$summary[turnover$summary$comparison %in% main, ],
    Fig5C_cross_layer_context = cross_data,
    provenance = list(clarification = .lt_marine_reporting_clarification(),
      Fig5B_count_scope = "persisted_terminal_worlds_not_ASR_branch_counts",
      Fig5B_qualifier = fig5b$provenance$qualifier,
      Phase13_qualifier = "historical_phase13_screening_mapping_mismatch_preserved",
      original_overview_scope = "Original no-analysis-rerun text describes historical table compilation, not the fresh run.",
      main_vs_cross_layer_preserved = TRUE, new_world_or_fit_generated = FALSE))
}
