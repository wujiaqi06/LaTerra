.lt_s2_numeric_tsv_roundtrip <- function(x) {
  text <- character()
  con <- textConnection("text", "w", local = TRUE)
  on.exit(close(con), add = TRUE)
  utils::write.table(data.frame(value = x), con, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
  utils::read.delim(text = paste(text, collapse = "\n"))$value
}

.lt_s2_iqr_text <- function(x, digits = 0L) {
  x <- x[is.finite(x)]
  if (!length(x)) return("")
  q <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f-%.", digits, "f"), q[[1]], q[[2]])
}

.lt_marine_nested_fold_table <- function(result, trait, global_features) {
  rows <- lapply(result$folds, function(f) {
    x <- f$identity; ntest <- length(x$test_species)
    overlap <- length(intersect(f$genes, global_features))
    skipped <- !ntest
    notes <- if (skipped) "No evaluable held-out species for this model after exclusion." else
      switch(f$status,
        null_intercept_only_training_one_class = "Training fold has only one response class.",
        null_intercept_only_no_FDR_features_in_training_ttest = "Training-only t-test/FDR selected zero candidate genes.",
        null_intercept_only_no_features_after_all_missing_drop = "No usable features remain after all-missing training drop.",
        null_intercept_only_no_features_after_zero_variance_drop = "No usable features remain after zero-variance training drop.",
        "Fold-specific training t-test/FDR features; fold-wise preprocessing; glmnet standardize=FALSE.")
    data.frame(trait = trait, model = paste0(trait, "_nested_ttest_baseline"),
      fold_id = x$fold_id, held_out_genus = x$genus,
      n_training_samples = length(x$train_species), n_test_samples = ntest,
      n_train_positive = sum(x$y_train == 1), n_train_negative = sum(x$y_train == 0),
      n_test_positive = sum(x$y_test == 1), n_test_negative = sum(x$y_test == 0),
      n_ttest_genes_tested = f$n_tested, n_fold_candidate_genes = length(f$genes),
      n_overlap_with_global_FDR_candidates = overlap,
      prop_nested_candidates_in_global_FDR = if (length(f$genes)) overlap / length(f$genes) else NA_real_,
      prop_global_FDR_recovered_in_nested_candidates = overlap / length(global_features),
      n_heldout_terminal_branches_removed_from_ttest = f$n_removed,
      n_branches_screened_in_training_ttest = f$n_screened,
      n_features_dropped_all_missing = f$n_dropped_all,
      n_features_dropped_zero_variance = f$n_dropped_zero,
      n_features_used_for_lasso = f$n_used, n_selected_predictors = length(f$beta),
      lambda_min = f$lambda_min, fold_model_status = f$status,
      corrected_preprocessing_boundary = if (skipped) "not_run_no_evaluable_test_species" else
        "training-only imputation/scaling/lambda/model fitting; terminal species only",
      fold_seed = x$seed, notes = notes, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

.lt_marine_nested_tables <- function(nested, phase11, models, traits, legacy_reference) {
  trait_ids <- c("marine_binary", "binary_aquatic_dependence")
  run_ids <- c("fix_marine_binary", "fix_aquatic_v2")
  if (length(nested) != 2L || length(phase11) != 2L)
    .lt_marine_abort("Both computed nested and Phase11 baseline results are required.", "dependency")
  oof <- folds <- auc <- selected <- phase11_evidence <- vector("list", 2L)
  for (j in 1:2) {
    id <- trait_ids[j]; model <- models[[run_ids[j]]]
    terminals <- model$terminal
    key <- match(terminals$species, traits$species)
    if (anyNA(key)) .lt_marine_abort("Missing explicit trait identity in nested output join.", "trait")
    metadata <- data.frame(species = terminals$species, genus = terminals$genus,
      marine_binary = traits$marine_binary[key],
      aquatic_v2 = traits$aquatic_v2_status[key],
      aquaticity_score = traits$aquaticity_score_sum_0_18[key], stringsAsFactors = FALSE)
    metadata$trait_category <- ifelse(metadata$marine_binary == 1, "marine",
      ifelse(metadata$aquatic_v2 == 1, "non-marine aquatic",
        ifelse(metadata$aquatic_v2 == 0.5, "semi-aquatic", "terrestrial")))
    fold <- .lt_marine_nested_fold_table(nested[[j]], id, model$global_screen$significant$gene)
    folds[[j]] <- fold
    oof_rows <- list(); selected_rows <- list()
    for (f in nested[[j]]$folds) {
      x <- f$identity
      if (length(x$test_species)) {
        m <- metadata[match(x$test_species, metadata$species), , drop = FALSE]
        oof_rows[[length(oof_rows) + 1L]] <- data.frame(species = m$species, genus = m$genus,
          model = paste0(id, "_nested_ttest_baseline"), trait_category = m$trait_category,
          marine_binary = m$marine_binary, binary_aquatic_endpoint = ifelse(m$aquatic_v2 == 0.5, NA, m$aquatic_v2),
          aquaticity_score = m$aquaticity_score, held_out_genus = x$genus,
          corrected_globalFDR_OOF_prediction = NA_real_, nested_ttest_OOF_prediction = f$prediction,
          prediction_delta = NA_real_, nested_OOF_available = TRUE, exclusion_reason = "",
          fold_n_candidate_features_FDR = length(f$genes), fold_n_features_after_preprocessing = f$n_used,
          fold_n_selected_predictors = length(f$beta), fold_model_status = f$status, stringsAsFactors = FALSE)
      }
      if (length(f$beta)) selected_rows[[length(selected_rows) + 1L]] <- data.frame(
        trait = id, model = paste0(id, "_nested_ttest_baseline"), fold_id = x$fold_id,
        held_out_genus = x$genus, gene = names(f$beta), coefficient = as.numeric(f$beta),
        coefficient_sign = ifelse(f$beta > 0, "positive", "negative"), lambda_min = f$lambda_min,
        stringsAsFactors = FALSE)
    }
    excluded <- metadata[terminals$response == 0.5, , drop = FALSE]
    if (nrow(excluded)) oof_rows[[length(oof_rows) + 1L]] <- data.frame(
      species = excluded$species, genus = excluded$genus, model = paste0(id, "_nested_ttest_baseline"),
      trait_category = excluded$trait_category, marine_binary = excluded$marine_binary,
      binary_aquatic_endpoint = NA_real_, aquaticity_score = excluded$aquaticity_score,
      held_out_genus = excluded$genus, corrected_globalFDR_OOF_prediction = NA_real_,
      nested_ttest_OOF_prediction = NA_real_, prediction_delta = NA_real_, nested_OOF_available = FALSE,
      exclusion_reason = "semi-aquatic aquatic_v2=0.5 excluded from binary aquatic-dependence LASSO/AUC",
      fold_n_candidate_features_FDR = NA_integer_, fold_n_features_after_preprocessing = NA_integer_,
      fold_n_selected_predictors = NA_integer_, fold_model_status = "excluded_no_binary_aquatic_endpoint",
      stringsAsFactors = FALSE)
    all_oof <- do.call(rbind, oof_rows)
    all_oof <- all_oof[order(all_oof$species), , drop = FALSE]
    ref <- do.call(rbind, lapply(phase11[[j]]$folds, function(f) {
      if (!length(f$prediction)) return(NULL)
      data.frame(species = f$identity$test_species, prediction = f$prediction,
        response = f$identity$y_test, stringsAsFactors = FALSE)
    }))
    .lt_assert_unique_ids(ref$species, "computed Phase11 species")
    # Phase11's AUC is computed before its probability table is serialized.
    phase_auc <- .lt_s2_numeric_tsv_roundtrip(.lt_s2_auc(ref$response, ref$prediction))
    # Preserve the source boundary: Phase11 writes then Phase12B reads these
    # numbers. Primary nested results are not rounded by this comparator join.
    ref$prediction <- .lt_s2_numeric_tsv_roundtrip(ref$prediction)
    ref <- ref[order(ref$species), , drop = FALSE]
    all_oof$corrected_globalFDR_OOF_prediction <- ref$prediction[match(all_oof$species, ref$species)]
    all_oof$prediction_delta <- all_oof$nested_ttest_OOF_prediction - all_oof$corrected_globalFDR_OOF_prediction
    eval <- all_oof[all_oof$nested_OOF_available, , drop = FALSE]
    y <- if (j == 1L) eval$marine_binary else eval$binary_aquatic_endpoint
    value <- .lt_s2_auc(y, eval$nested_ttest_OOF_prediction)
    # Narrow reporting-only exception: this field never reaches fit functions.
    old <- legacy_reference$gLOOCV_AUC[match(run_ids[j], legacy_reference$run_id)]
    if (length(old) != 1L || is.na(old)) .lt_marine_abort("Missing admitted historical AUC report reference.", "reference")
    contributing <- fold$n_test_samples > 0
    auc[[j]] <- data.frame(model = paste0(id, "_nested_ttest_baseline"),
      trait_axis = c("marine specialization", "binary aquatic dependence")[j],
      n_species_evaluated = nrow(eval), n_positive = sum(y == 1, na.rm = TRUE), n_negative = sum(y == 0, na.rm = TRUE),
      semi_aquatic_handling = if (j == 1L) "retained as non-marine/background 0 for marine binary model" else
        "semi-aquatic aquatic_v2=0.5 labels excluded from training, testing, imputation, scaling, and AUC",
      AUC_legacy_global_preprocessing = as.numeric(old), AUC_corrected_globalFDR_foldwise_preprocessing = phase_auc,
      AUC_nested_ttest = value, delta_vs_corrected_globalFDR = value - phase_auc,
      n_folds = sum(contributing),
      n_folds_zero_candidate_features = sum(contributing & fold$n_fold_candidate_genes == 0, na.rm = TRUE),
      n_folds_null_intercept_only = sum(contributing & grepl("^null_intercept_only", fold$fold_model_status)),
      median_candidate_features_per_fold = stats::median(fold$n_fold_candidate_genes[contributing], na.rm = TRUE),
      IQR_candidate_features_per_fold = .lt_s2_iqr_text(fold$n_fold_candidate_genes[contributing]),
      median_selected_predictors_per_fold = stats::median(fold$n_selected_predictors[contributing], na.rm = TRUE),
      IQR_selected_predictors_per_fold = .lt_s2_iqr_text(fold$n_selected_predictors[contributing]),
      evidence_level = "terminal-species genus-level LOOCV with supervised t-test/FDR candidate discovery repeated inside each outer training fold; global deterministic ASR branch states retained; fold-wise training-only imputation/scaling/lambda/model fitting",
      trait = id, comparison_model = paste0(id, "_corrected_preprocessing_gLOOCV"), stringsAsFactors = FALSE)
    oof[[j]] <- all_oof; selected[[j]] <- do.call(rbind, selected_rows)
    phase11_evidence[[j]] <- list(trait = id, serialized_computed_oof = ref, serialized_computed_auc = phase_auc)
  }
  list(auc = do.call(rbind, auc), oof = do.call(rbind, oof),
    folds = do.call(rbind, folds), selected = do.call(rbind, selected),
    phase11 = phase11_evidence,
    provenance = list(actual_nested = "computed", phase11 = "computed_auxiliary",
      legacy_AUC = "historical_reference_imported_not_recomputed",
      prediction_delta = "computed_nested_minus_serialized_computed_Phase11"))
}
