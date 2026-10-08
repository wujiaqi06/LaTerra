.lt_marine_sensitivity_ids <- function() {
  c("fix_marine_binary", "fix_drop_whale", "fix_whale_only", "fix_pinniped_only",
    "fix_aquatic_v2", "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia",
    "fix_aquatic_v2_noMarineEdge")
}

.lt_marine_sensitivity_tables <- function(results) {
  ids <- .lt_marine_sensitivity_ids()
  if (!identical(names(results), ids))
    .lt_marine_abort("Sensitivity results must retain all eight frozen model identities/order.", "model_order")
  # Literal submitted Figure5A display metadata. These labels are historical
  # reproduction fields, not a new result-dependent biological adjudication.
  labels <- c("Marine baseline", "Drop cetaceans", "Cetacea only", "Pinnipedia only",
    "Aquatic baseline", "No Cetacea", "No Pinnipedia", "No marine-edge taxa")
  interpretation <- c("robust sparse prediction", "substantially attenuated; retained",
    "strongest marine sparse prediction", "weaker marine sparse prediction",
    "moderate sparse prediction", "retained discrimination", "retained, stronger discrimination",
    "retained discrimination")
  out <- vector("list", length(ids)); detail <- vector("list", length(ids))
  for (j in seq_along(ids)) {
    folds <- results[[j]]$folds
    contributing <- vapply(folds, function(x) length(x$identity$test_species) > 0L, logical(1))
    p <- unlist(lapply(folds, `[[`, "prediction"), use.names = FALSE)
    y <- unlist(lapply(folds, function(x) x$identity$y_test), use.names = FALSE)
    selected <- vapply(folds[contributing], function(x) length(x$beta), integer(1))
    med <- stats::median(selected)
    iqr <- .lt_s2_iqr_text(selected)
    auc <- .lt_s2_auc(y, p)
    out[[j]] <- data.frame(trait_axis = if (j <= 4L) "Marine specialization" else "Binary aquatic dependence",
      run = labels[j], trait_contrast = paste(sum(y == 1), "/", sum(y == 0)),
      # Python round(..., 3) in the historical table only; unrounded evidence
      # stays in the computed fold/OOF object and auxiliary table below.
      nested_AUC = round(.lt_s2_numeric_tsv_roundtrip(auc), 3),
      # Preserve the source's UTF-8 bytes even with LC_CTYPE=C; marked UTF-8
      # strings are otherwise escaped as literal <U+2013> by write.table.
      median_selected_predictors_per_fold_IQR = paste0(as.integer(med), " (",
        gsub("-", "\u2013", iqr, fixed = TRUE, useBytes = TRUE), ")"),
      model_status_interpretation = interpretation[j], run_id = ids[j],
      n_species_evaluated = length(y), n_positive = sum(y == 1), n_negative = sum(y == 0),
      n_folds = sum(contributing), n_folds_no_evaluable_test_species = sum(!contributing),
      evidence_level = "nested supervised t-test feature-selection gLOOCV; GBI/ASR frozen; fold-wise training-terminal-only imputation/scaling/lambda/model fitting",
      stringsAsFactors = FALSE)
    detail[[j]] <- data.frame(run_id = ids[j], AUC_unrounded = auc,
      median_selected_predictors = med, IQR_selected_predictors = iqr,
      n_null_intercept_only = sum(vapply(folds[contributing], function(x)
        startsWith(x$status, "null_intercept_only"), logical(1))),
      n_no_selected = sum(selected == 0L), stringsAsFactors = FALSE)
  }
  list(table = do.call(rbind, out), unrounded = do.call(rbind, detail),
    provenance = list(rounded_table_digits = 3L, raw_results_preserved = TRUE,
      display_interpretations = "historical_submitted_labels_not_new_classification",
      stage = "nested_phase13", no_baseline_fit_substitution = TRUE))
}
