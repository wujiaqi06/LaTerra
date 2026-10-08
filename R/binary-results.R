# Saved-result views only: never call baseline, ASR, screen or fitting operators.
.lt_binary_files <- function() c("supplied_inputs.rds", "configuration.rds", "S2_baseline_GBI.rds",
  "branch_ledger.tsv", "branch_states.tsv", "gene_effect.tsv", "branch_effect.tsv",
  "screen_tested.tsv", "screen_significant.tsv", "screen_gene_ledger.tsv", "screen_summary.tsv",
  "validation.rds", "OOF_predictions.tsv", "validation_summary.tsv", "fold_ledger.tsv",
  "seed_ledger.tsv", "input_hashes.tsv", "sessionInfo.txt", "RUN_STATUS.yaml", "run.rds")

.lt_binary_saved_table <- function(file) {
  tab <- utils::read.delim(file, colClasses = "character", check.names = FALSE,
    quote = "", comment.char = "", na.strings = character())
  keys <- c("gene", "gene_id", "branch", "branch_id", "species", "genus", "group", "fold_id", "run_id",
    "trait", "trait_id", "trait_name", "canonical_split_key", "branch_type", "terminal_taxon",
    "taxon_id", "scope", "role", "path", "sha256")
  for (n in setdiff(names(tab), keys)) tab[[n]] <- utils::type.convert(tab[[n]], as.is = TRUE)
  if (all(c("branch_type", "terminal_taxon") %in% names(tab))) {
    if (any(!tab$branch_type %in% c("terminal", "internal")))
      stop("Saved branch ledger requires terminal/internal branch_type values.", call. = FALSE)
    internal <- tab$branch_type == "internal"
    if (any(internal & !tab$terminal_taxon %in% c("", "NA")))
      stop("Saved branch ledger has an unexpected terminal_taxon on an internal branch; no identity is silently erased.", call. = FALSE)
    # Only the writer's internal missing sentinel is decoded. A terminal named
    # literally NA remains a scientific key, as do 001, TRUE and quoted keys.
    tab$terminal_taxon[internal] <- NA_character_
  }
  tab
}

.lt_binary_integrity <- function(path, required, verify) {
  absent <- required[!file.exists(file.path(path, required))]
  if (length(absent)) stop("Completed-run integrity error: missing ", paste(absent, collapse = ", "), call. = FALSE)
  manifest <- file.path(path, "SHA256SUMS")
  if (!file.exists(manifest)) stop("Completed-run integrity error: missing SHA256SUMS.", call. = FALSE)
  if (!verify) return(invisible(NULL))
  lines <- readLines(manifest, warn = FALSE)
  if (!length(lines) || any(!grepl("^[0-9a-f]{64}  [^/\\\\]+$", lines)))
    stop("Completed-run integrity error: malformed manifest.", call. = FALSE)
  names <- substring(lines, 67L); hashes <- substring(lines, 1L, 64L)
  if (anyDuplicated(names) || !all(required %in% names) || any(names %in% c(".", "..", "SHA256SUMS")) ||
      !setequal(names, setdiff(list.files(path, all.files = TRUE, no.. = TRUE), "SHA256SUMS")))
    stop("Completed-run integrity error: manifest inventory differs from run files.", call. = FALSE)
  actual <- vapply(file.path(path, names), .lt_marine_sha, character(1))
  if (!identical(unname(actual), unname(hashes))) stop("Completed-run integrity error: output hash mismatch.", call. = FALSE)
  invisible(NULL)
}

lt_read_binary_run <- function(run, verify = TRUE) {
  if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) stop("verify must be TRUE or FALSE.", call. = FALSE)
  path <- if (inherits(run, "lt_run")) run$results$output_dir else run
  .lt_assert_scalar_character(path, "run directory")
  path <- normalizePath(path, mustWork = TRUE)
  status_file <- file.path(path, "RUN_STATUS.yaml")
  if (!file.exists(status_file)) stop("Missing RUN_STATUS.yaml; not a saved binary run.", call. = FALSE)
  status <- yaml::read_yaml(status_file)
  completed <- identical(status$status, "COMPUTED_USER_S2_NOT_CERTIFIED")
  if (!completed && !identical(status$status, "FAILED")) stop("Not a supported user S2 run status.", call. = FALSE)
  if (completed) .lt_binary_integrity(path, .lt_binary_files(), verify)
  load <- function(name) {
    file <- file.path(path, name)
    if (!file.exists(file)) return(NULL)
    if (endsWith(name, ".rds")) readRDS(file) else .lt_binary_saved_table(file)
  }
  supplied <- load("supplied_inputs.rds"); baseline <- load("S2_baseline_GBI.rds")
  validation <- load("validation.rds")
  tables <- stats::setNames(lapply(c("branch_ledger", "branch_states", "gene_effect", "branch_effect",
    "screen_tested", "screen_significant", "screen_gene_ledger", "screen_summary", "OOF_predictions",
    "validation_summary", "fold_ledger", "seed_ledger", "input_hashes"), function(n) load(paste0(n, ".tsv"))),
    c("branch_ledger", "branch_states", "gene_effect", "branch_effect", "screen_tested", "screen_significant",
      "screen_gene_ledger", "screen_summary", "OOF_predictions", "validation_summary", "fold_ledger", "seed_ledger", "input_hashes"))
  diagnostics <- list()
  if (!is.null(supplied)) {
    qc <- lt_qc_matrix(supplied$matrix)
    diagnostics$input_cells <- qc$cell_summary$state_counts
    diagnostics$input_genes <- qc$gene_summary
    diagnostics$input_branches <- qc$branch_summary
    diagnostics$input_matrix <- qc$matrix_summary$overall
    diagnostics$input_distribution <- qc$matrix_summary$value_distribution
    diagnostics$input_domains <- qc$matrix_summary$metric_domains
    diagnostics$input_point_masses <- qc$matrix_summary$point_masses
    diagnostics$input_exact_extrema <- qc$matrix_summary$extreme_exact_repetition
    diagnostics$trait_support <- as.data.frame(table(state = supplied$trait$values), stringsAsFactors = FALSE)
  }
  if (!is.null(baseline)) {
    d <- baseline$gene_effect
    d$trim_cutoff_q0975_type7 <- unname(baseline$trim_cutoff[d$gene])
    d$n_branch_denominator <- ncol(baseline$trimmed_values)
    d$n_gbi_available <- rowSums(is.finite(baseline$gbi))
    d$n_denominator_unavailable <- rowSums(is.finite(baseline$trimmed_values) & !is.finite(baseline$gbi))
    d$removed_fraction_available_before <- ifelse(d$n_present_before_trim > 0, d$n_trimmed / d$n_present_before_trim, NA_real_)
    d$cutoff_status <- ifelse(is.finite(d$trim_cutoff_q0975_type7), "computed", "not_estimable")
    diagnostics$trim_genes <- d
  }
  if (!is.null(tables$branch_states)) diagnostics$branch_state_support <- as.data.frame(
    table(state = tables$branch_states$state), stringsAsFactors = FALSE)
  if (!is.null(tables$screen_gene_ledger)) {
    diagnostics$screen_support <- tables$screen_gene_ledger
    diagnostics$screen_support$disposition <- ifelse(!diagnostics$screen_support$tested, "not_estimable",
      ifelse(diagnostics$screen_support$FDR_selected, "selected", "not_selected"))
  }
  if (!is.null(validation)) {
    diagnostics$folds <- do.call(rbind, lapply(validation$folds, function(f) data.frame(
      fold_id = f$identity$fold_id, genus = f$identity$genus, seed = f$identity$seed, status = f$status,
      n_train_0 = sum(f$identity$y_train == 0), n_train_1 = sum(f$identity$y_train == 1),
      n_test_0 = sum(f$identity$y_test == 0), n_test_1 = sum(f$identity$y_test == 1),
      n_predicted = sum(is.finite(f$prediction)), n_tested_genes = f$n_tested,
      n_FDR_genes = length(f$genes), n_removed_branches = f$n_removed, n_screened_branches = f$n_screened,
      n_dropped_all_missing = f$n_dropped_all, n_dropped_zero_variance = f$n_dropped_zero,
      n_design_features = f$n_used, n_nonzero_coefficients = length(f$beta), solver_status = f$solver_status,
      warnings = paste(unique(c(f$warnings, f$fit_warnings)), collapse = " | "), stringsAsFactors = FALSE)))
    diagnostics$feature_usage <- do.call(rbind, lapply(validation$folds, function(f) {
      if (!length(f$genes)) return(NULL)
      data.frame(fold_id = f$identity$fold_id, gene = f$genes,
        selected_in_screen = TRUE, nonzero_coefficient = f$genes %in% names(f$beta), stringsAsFactors = FALSE)
    }))
    if (is.null(diagnostics$feature_usage)) diagnostics$feature_usage <- data.frame(
      fold_id = character(), gene = character(), selected_in_screen = logical(), nonzero_coefficient = logical())
    p <- tables$OOF_predictions
    eligible <- sum(supplied$trait$values %in% c(0, 1))
    diagnostics$oof_coverage <- data.frame(n_eligible_terminals = eligible,
      n_prediction_rows = if (is.null(p)) 0L else nrow(p),
      n_unique_predicted_terminals = if (is.null(p)) 0L else length(unique(p$species[is.finite(p$probability)])),
      AUC_status = if (is.null(tables$validation_summary)) "not_run_or_incomplete" else
        if (is.finite(tables$validation_summary$AUC)) "computed" else "not_estimable")
  }
  modules <- c("input", "trim", "screen", "validation")
  failed_module <- switch(if (is.null(status$stage)) "" else status$stage, baseline = "trim", screen = "screen",
    grouped_validation = "validation", "")
  diagnostics$dependencies <- data.frame(module = modules,
    status = ifelse(c(!is.null(supplied), !is.null(baseline), !is.null(tables$screen_summary), !is.null(validation)),
      "saved_outputs_available", ifelse(modules == failed_module, "incomplete_failed_stage", "not_run")), stringsAsFactors = FALSE)
  structure(list(run_dir = path, status = status, completed = completed,
    integrity = if (completed && verify) "verified_file_hashes" else "not_verified",
    configuration = load("configuration.rds"), tables = tables, diagnostics = diagnostics), class = "lt_binary_result")
}

summary.lt_binary_result <- function(object, ...) {
  list(status = object$status, integrity = object$integrity,
    configuration = object$configuration, screen = object$tables$screen_summary,
    validation = object$tables$validation_summary, oof_coverage = object$diagnostics$oof_coverage,
    dependencies = object$diagnostics$dependencies, source = object$run_dir)
}

print.lt_binary_result <- function(x, ...) {
  cat("La Terra saved binary results\n", x$status$status, " | ", x$integrity, "\n", x$run_dir, "\n", sep = "")
  if (!is.null(x$tables$screen_summary)) print(x$tables$screen_summary)
  if (!is.null(x$tables$validation_summary)) print(x$tables$validation_summary)
  invisible(x)
}

lt_report_binary <- function(run, output_dir) {
  # Always reload and verify; do not trust a cached view of a modified run.
  source <- if (inherits(run, "lt_binary_result")) run$run_dir else run
  x <- lt_read_binary_run(source)
  .lt_assert_scalar_character(output_dir, "output_dir")
  if (file.exists(output_dir)) stop("Report output_dir exists; choose a new directory.", call. = FALSE)
  parent <- normalizePath(dirname(output_dir), mustWork = TRUE)
  destination <- file.path(parent, basename(output_dir))
  if (identical(parent, x$run_dir) || startsWith(parent, paste0(x$run_dir, .Platform$file.sep)))
    stop("Reports must be outside the scientific run directory.", call. = FALSE)
  if (!dir.create(destination)) stop("Cannot create report directory.", call. = FALSE)
  sections <- list(
    D1_input = c("input_cells", "input_matrix", "input_genes", "input_branches", "input_distribution", "input_domains", "input_point_masses", "input_exact_extrema"),
    D2_recipe = "trim_genes", D3_screen = c("trait_support", "branch_state_support", "screen_support"),
    D4_validation = c("folds", "feature_usage", "oof_coverage", "dependencies"))
  lines <- c("# La Terra saved binary result report", "", paste("Source:", x$run_dir),
    paste("Status:", x$status$status), paste("Integrity:", x$integrity),
    "", "Descriptive evidence only. No refit, automatic exclusion, quality score or biological/release certification.")
  if (!x$completed) lines <- c(lines, "", paste("PARTIAL failed run; stage:", x$status$stage),
    paste("Failure:", x$status$message), "Dependencies distinguish not_run from incomplete_failed_stage; absent partial products are not_run_or_incomplete, not optional completed-run missingness.")
  lines <- c(lines, "", "## Input identity locations", "",
    "input_hashes.tsv indexes the supplied_inputs.rds snapshot, not the original input files.",
    "For file-backed inputs, read supplied_inputs.rds and inspect:",
    "- matrix$metadata$splitaligner_import[c('input_file', 'input_file_sha256')] (exchange file)",
    "- trait$metadata$input_source (trait TSV)",
    "- attr(terminal_groups, 'input_source') (groups TSV)",
    "- attr(folds, 'input_source') (folds TSV)",
    "Objects supplied directly in memory need not have an original file identity; no file hash is invented.")
  notes <- c(
    "Raw coordinate/value/eligibility availability; all-cell and available-value denominators are distinct. Standard lt_qc_matrix is threshold-free. Exact repetitions do not establish upstream causes.",
    "Selected submitted_S2 footprint: inclusive >= gene q0.975/type7; arithmetic GE/BE (medians auxiliary); divide BE then GE. Removed fractions are observed, not guaranteed 2.5%. Finite trimmed cells with unavailable GBI report denominator-domain loss, not new exclusions.",
    "Explicit 0 reference, 1 focal, 0.5 excluded; historical castor descendant-node ASR. Welch t and beta are reference-minus-focal; strict-prefix FDR0.01 is preserved (not substituted BH). Screen ledger has finite group sizes and untested reasons. FDR nonselection is not an unavailable test.",
    "Saved held-out terminal probabilities for class 1 and AUC. Fold n_removed/n_screened are BRANCH counts, not gene counts. Feature flow: tested -> screen selected -> all-missing/zero-variance drops -> design -> nonzero coefficients. NA special-fold fields stay NA. Fold usage is descriptive, not independent replication; per-cell imputation and per-gene design-drop reasons were not persisted.")
  science_tables <- list(D1_input = character(), D2_recipe = c("gene_effect", "branch_effect"),
    D3_screen = c("branch_ledger", "branch_states", "screen_summary", "screen_gene_ledger", "screen_tested", "screen_significant"),
    D4_validation = c("OOF_predictions", "validation_summary", "fold_ledger", "seed_ledger", "input_hashes"))
  for (i in seq_along(sections)) {
    section <- names(sections)[i]
    lines <- c(lines, "", paste0("## ", section), "", notes[i], "")
    for (n in sections[[i]]) {
      tab <- x$diagnostics[[n]]
      if (is.data.frame(tab)) {
        name <- paste0(section, "_", n, ".tsv")
        .lt_marine_write_tsv(tab, file.path(destination, name))
        lines <- c(lines, paste0("- [", n, "](", name, ")"))
      } else lines <- c(lines, paste0("- ", n, ": not_run_or_incomplete"))
    }
    for (n in science_tables[[section]]) if (!is.null(x$tables[[n]])) {
      name <- paste0(n, ".tsv")
      .lt_marine_write_tsv(x$tables[[n]], file.path(destination, name))
      lines <- c(lines, paste0("- [", n, "](", name, ")"))
    }
  }
  saveRDS(summary(x), file.path(destination, "summary.rds"), version = 3)
  writeLines(lines, file.path(destination, "report.md"))
  files <- list.files(destination, full.names = TRUE)
  writeLines(paste(vapply(files, .lt_marine_sha, character(1)), basename(files), sep = "  "), file.path(destination, "SHA256SUMS"))
  invisible(normalizePath(destination, mustWork = TRUE))
}

lt_plot_binary <- function(run, view = c("oof", "roc", "trim")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 for editable statistical plots.", call. = FALSE)
  view <- match.arg(view)
  x <- lt_read_binary_run(if (inherits(run, "lt_binary_result")) run$run_dir else run)
  caption <- function(file, note) paste(strwrap(paste0("Source: ", basename(x$run_dir), "/", file,
    "; ", note, "; saved outputs only; no refit"), width = 80), collapse = "\n")
  with_source <- function(plot) {
    if (identical(x$completed, FALSE)) {
      stage <- x$status$stage
      if (!is.character(stage) || length(stage) != 1L || is.na(stage) || !nzchar(stage)) stage <- "not recorded"
      marker <- paste0("PARTIAL / FAILED (stage: ", stage, ")")
      plot <- plot + ggplot2::labs(title = paste(marker, plot$labels$title, sep = "\n"),
        caption = paste(marker, plot$labels$caption, sep = "\n"))
    }
    attr(plot, "lt_source_run_dir") <- x$run_dir
    plot
  }
  if (view == "trim") {
    d <- x$diagnostics$trim_genes
    if (is.null(d)) stop("Trim outputs not_run_or_incomplete in this partial run.", call. = FALSE)
    return(with_source(ggplot2::ggplot(d, ggplot2::aes(x = n_present_before_trim, y = n_trimmed)) + ggplot2::geom_point() +
      ggplot2::labs(x = "Available input branches per gene (count)", y = "Actually removed branches per gene (count)",
        title = "Inclusive gene q0.975/type7 recipe footprint",
        caption = caption("S2_baseline_GBI.rds", "denominator: available input branches in each gene"))))
  }
  d <- x$tables$OOF_predictions
  if (is.null(d)) stop("OOF outputs not_run_or_incomplete in this partial run.", call. = FALSE)
  d <- d[is.finite(d$probability) & d$response %in% c(0, 1), , drop = FALSE]
  if (!nrow(d)) stop("OOF figure not_estimable: no finite eligible held-out probabilities.", call. = FALSE)
  if (view == "oof") return(with_source(ggplot2::ggplot(d, ggplot2::aes(x = factor(response), y = probability)) +
    ggplot2::geom_point(ggplot2::aes(group = genus), shape = 16,
      position = ggplot2::position_dodge(width = 0.3)) +
    ggplot2::labs(x = "Explicit terminal class (0 reference / 1 focal)", y = "Held-out probability of class 1 (unitless)",
      title = "Saved grouped held-out predictions",
      subtitle = paste(nrow(d), "finite eligible prediction rows across", length(unique(d$genus)), "supplied groups"),
      caption = caption("OOF_predictions.tsv", "denominator: finite eligible prediction rows; groups determine horizontal offsets"))))
  np <- sum(d$response == 1); nn <- sum(d$response == 0)
  if (!np || !nn) stop("ROC not_estimable: both classes required.", call. = FALSE)
  d <- d[order(d$probability, decreasing = TRUE), , drop = FALSE]
  last <- !duplicated(d$probability, fromLast = TRUE)
  roc <- data.frame(FPR = c(0, cumsum(d$response == 0)[last] / nn),
    TPR = c(0, cumsum(d$response == 1)[last] / np))
  with_source(ggplot2::ggplot(roc, ggplot2::aes(x = FPR, y = TPR)) + ggplot2::geom_path() +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "False positive fraction (reference terminals)", y = "True positive fraction (focal terminals)",
      title = "Saved OOF ROC (larger probability predicts class 1)",
      subtitle = paste("AUC =", format(.lt_s2_auc(d$response, d$probability))),
      caption = caption("OOF_predictions.tsv", "pooled finite eligible held-out rows; not independent fold replication")))
}

utils::globalVariables(c("n_present_before_trim", "n_trimmed", "response", "probability", "genus", "FPR", "TPR"))
