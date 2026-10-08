.lt_marine_target_artifacts <- function(workbook_names) {
  # Computation/evidence locations, not expected answers or acceptance overrides.
  science <- list(
    c("M1/raw_matrix.rds", "M1/upstream_tokens.rds", "M1/branch_ledger.tsv"),
    c("M1/S2_baseline_GBI.rds", "M1/GBI_consumption_boundary.rds", "M1/GBI.generated_consumed.oldlabels.tsv"),
    "M1/marine_binary.branch_states.tsv", "M1/aquatic_v2.branch_states.tsv",
    c("M1/marine_binary.tested.tsv", "M1/marine_binary.significant.tsv", "M1/marine_binary.gene_ledger.tsv"),
    c("M1/aquatic_v2.tested.tsv", "M1/aquatic_v2.significant.tsv", "M1/aquatic_v2.gene_ledger.tsv"),
    c("nested_phase12B_fix_marine_binary.rds", "nested_auc.tsv", "nested_oof.tsv", "nested_folds.tsv", "nested_selected.tsv"),
    c("nested_phase12B_fix_aquatic_v2.rds", "nested_auc.tsv", "nested_oof.tsv", "nested_folds.tsv", "nested_selected.tsv"),
    c("Fig5A_sensitivity.tsv", "Fig5A_unrounded.tsv", paste0("nested_phase13_", .lt_marine_sensitivity_ids(), ".rds")),
    c("Fig5B.rds", "M2_input_hashes.tsv"), c("Fig5B_screening.tsv", "Fig5B_observed_vs_null.tsv"),
    c("full_data.rds", "full_coefficients.tsv", "full_summary.tsv"),
    c("full_data.rds", "full_coefficients.tsv", "full_summary.tsv"),
    "full_partition.tsv", c("display_annotation.tsv", "display_counts.tsv", "display_provenance.rds"),
    c("turnover_observed.tsv", "turnover_universes.tsv"),
    c("turnover.rds", "turnover_worlds.tsv", "turnover_summary.tsv"), "nc_projection.tsv",
    c("nc_fingerprints.tsv", "focal_projections.tsv", "all_fingerprints.tsv"),
    c("internal_projections.tsv", "internal_targets.tsv"))
  rows <- lapply(seq_along(science), function(i) data.frame(target_id = sprintf("TGT-%03d", i),
    source_root = "scientific", path = science[[i]], stringsAsFactors = FALSE))
  rows[[21L]] <- data.frame(target_id = "TGT-021", source_root = "reporting",
    path = c(unname(workbook_names), "logical_sheets.rds", "reporting_provenance.rds"))
  rows[[22L]] <- data.frame(target_id = "TGT-022", source_root = c(rep("scientific", 5L), rep("reporting", 2L)),
    path = c("run.rds", "input_hashes.tsv", "seed_ledger.tsv", "fold_ledger.tsv", "M2_authority.rds",
      "reporting_provenance.rds", "reporting_recipe.rds"))
  result <- do.call(rbind, rows)
  result$required <- TRUE
  result$computation_status <- "COMPUTED_NOT_REPLAY_VALIDATED"
  result$independent_acceptance <- "NOT_CLAIMED"
  result
}

.lt_marine_complete_reporting_run <- function(reader, reporting, run_dir, output_dir, workbook_names) {
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  manifest_path <- file.path(run_dir, "SHA256SUMS")
  .lt_marine_verify(manifest_path, reader$manifest_sha256, "unchanged_scientific_run_manifest")
  lines <- readLines(manifest_path)
  science <- data.frame(role = paste0("scientific/", substring(lines, 67L)),
    path = substring(lines, 67L), sha256 = substr(lines, 1L, 64L), stringsAsFactors = FALSE)
  # Revalidate the entire declared source output set at the terminal boundary,
  # including linked targets that the worksheet assembler itself did not read.
  for (i in seq_len(nrow(science))) {
    p <- normalizePath(file.path(run_dir, science$path[i]), mustWork = TRUE)
    if (!startsWith(p, paste0(run_dir, "/")))
      .lt_marine_abort("Scientific target evidence escapes its source root.", "reporting_dependency")
    .lt_marine_verify(p, science$sha256[i], paste0("source_target_", science$path[i]))
  }
  recipe <- list(profile = "marine_NC_submitted", source_configuration_sha256 = reader$run$provenance$configuration_sha256,
    scientific_run_manifest_sha256 = reader$manifest_sha256,
    reporting_authority = reporting$reporting_authority, clarification = reporting$clarification,
    precision_clarification = reporting$precision_clarification,
    no_new_scientific_calculation = TRUE, no_archived_result_workbook_input = TRUE,
    required_targets = sprintf("TGT-%03d", 1:22), independent_acceptance = "NOT_CLAIMED")
  saveRDS(recipe, file.path(output_dir, "reporting_recipe.rds"), version = 3)
  targets <- .lt_marine_target_artifacts(workbook_names)
  targets$sha256 <- vapply(seq_len(nrow(targets)), function(i) {
    base <- if (targets$source_root[i] == "scientific") run_dir else output_dir
    if (targets$source_root[i] == "scientific" && !targets$path[i] %in% science$path)
      .lt_marine_abort("Required target evidence is not declared by the scientific run.", "reporting_dependency")
    .lt_marine_sha(file.path(base, targets$path[i]))
  }, "")
  .lt_marine_write_tsv(targets, file.path(output_dir, "target_evidence_ledger.tsv"))
  target_status <- data.frame(target_id = sprintf("TGT-%03d", 1:22), required = TRUE,
    computation_status = "COMPUTED_NOT_REPLAY_VALIDATED", validation_status = "NOT_RUN_BY_REPORTING_ENTRY",
    independent_acceptance = "NOT_CLAIMED", evidence_ledger = "target_evidence_ledger.tsv")
  .lt_marine_write_tsv(target_status, file.path(output_dir, "target_status.tsv"))
  p <- reader$run$provenance
  p$input_sha256 <- rbind(p$input_sha256, data.frame(role = c("scientific_source_run", "scientific_source_manifest"),
    path = c(file.path(run_dir, "run.rds"), manifest_path),
    sha256 = c(.lt_marine_sha(file.path(run_dir, "run.rds")), reader$manifest_sha256)))
  p$configuration_sha256 <- .lt_marine_sha(file.path(output_dir, "reporting_recipe.rds"))
  p$software_version <- list(package = "LaTerra", version = reporting$reporting_package_version,
    installed_runtime_files = reporting$reporting_installed_runtime_files,
    scientific_source = reader$run$provenance$software_version, implementation_certification = "NOT_CLAIMED")
  p$environment <- list(scientific = reader$run$provenance$environment, reporting = reporting$reporting_environment)
  # The scientific layer's relative references remain scoped to the source root.
  # Do not relabel them as files inside the terminal-report directory.
  p$eligibility_masks <- list(source_root = "scientific", layers = reader$run$provenance$eligibility_masks)
  p$commands <- c(p$commands, reporting$command)
  authorities <- c(stats::setNames(reporting$reporting_authority$authority_sha256,
    reporting$reporting_authority$authority_id), stats::setNames(reporting$clarification$sha256, reporting$clarification$id),
    stats::setNames(reporting$precision_clarification$sha256, reporting$precision_clarification$id))
  p$specification_ledger <- rbind(p$specification_ledger,
    data.frame(specification_id = names(authorities), kind = "frozen_terminal_reporting_authority",
      version = "TARGET001", definition = "reporting_provenance.rds", sha256 = unname(authorities)))
  report_files <- list.files(output_dir, full.names = TRUE)
  scientific_outputs <- science
  scientific_outputs$path <- file.path(run_dir, scientific_outputs$path)
  p$output_hashes <- rbind(scientific_outputs,
    data.frame(role = paste0("reporting/", basename(report_files)), path = report_files,
      sha256 = vapply(report_files, .lt_marine_sha, "")))
  result <- lt_run(paste0(reader$run$run_id, "::reporting"), provenance = p, status = "completed",
    results = list(profile = "marine_NC_submitted", milestone = "M3_reporting",
      scientific_run_id = reader$run$run_id, scientific_source_run_sha256 = utils::tail(p$input_sha256$sha256, 2L)[1L],
      scientific_source_configuration_sha256 = reader$run$provenance$configuration_sha256,
      dependency_roots = list(scientific = run_dir, reporting = output_dir),
      target_status = target_status, target_evidence = targets,
      replay_certification = "NOT_CLAIMED", scientific_results_modified = FALSE,
      self_hash_rule = "Final SHA256SUMS binds this run.rds; its own output ledger excludes self."))
  saveRDS(result, file.path(output_dir, "run.rds"), version = 3)
  result
}
