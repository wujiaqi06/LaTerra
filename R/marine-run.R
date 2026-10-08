#' Inspect the study-specific frozen Marine/S2 profile
#'
#' This is the TARGET001 M1 recipe, not a generic QC or statistical default.
#' The returned list is for inspection; editing it does not change the explicit
#' frozen workflow implemented by [lt_run_marine()].
#' @return A named list containing scientific choices and input identities.
#' @export
lt_marine_profile <- function() {
  path <- system.file("profiles", "marine_S2_frozen_m1.yaml", package = "LaTerra")
  if (!nzchar(path)) .lt_marine_abort("Installed Marine profile is missing.", "installation")
  .lt_marine_verify(path, "548e532d2bd880b89b17941f4e8576e53b01153caa1d20607ea9805283eac9c8",
    "canonical_M1_profile_not_caller_supplied_authority")
  yaml::read_yaml(path, eval.expr = FALSE)
}

.lt_marine_write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

#' Run the minimal frozen Marine/S2 workflow
#'
#' Starting from the frozen raw matrix, independently computes S2 arithmetic
#' marginals and GBI, deterministic marine/aquatic branch annotations, and the
#' two baseline Welch screens. This M1 entry point does not implement M2 stages
#' or certify all 22 TARGET001 gates. It never reads frozen GBI or screen output
#' authorities; a separate replay validator compares outputs after execution.
#'
#' Input SHA checks are specific to this historical reproduction profile, not a
#' validity requirement for ordinary user-facing La Terra objects.
#'
#' @param data_root Directory containing the frozen Dryad layout, including
#'   `17_large_matrix_archives`, `00_traits_and_species`, and the approved
#'   complete-state ledger in `supplemental_authority`. See the installed
#'   profile for exact relative paths. No output-oracle inputs are required.
#' @param output_dir New output directory. Existing files/directories are never
#'   overwritten. Do not place it within `data_root`.
#' @param environment `"frozen"` requires the recorded R/ape/castor stack.
#'   `"compatibility"` is an explicit, warned non-certifying execution mode.
#' @return A completed `lt_run` with output paths and provenance. Completion is
#'   computational status only, not scientific or target-replay certification.
#' @export
lt_run_marine <- function(data_root, output_dir,
                          environment = c("frozen", "compatibility")) {
  environment <- match.arg(environment)
  .lt_assert_scalar_character(data_root, "data_root")
  .lt_assert_scalar_character(output_dir, "output_dir")
  if (!dir.exists(data_root)) .lt_marine_abort("`data_root` must be an existing directory.", "missing_input")
  data_root <- normalizePath(data_root, mustWork = TRUE)
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    .lt_marine_abort("`output_dir` already exists. Choose a new run directory; no overwrite is permitted.", "output_exists")
  }
  profile <- lt_marine_profile()
  env <- .lt_marine_environment(profile, environment)
  extraction <- tempfile("laterra-marine-m1-")
  if (!dir.create(extraction)) .lt_marine_abort("Cannot create private extraction directory.", "io")
  on.exit(unlink(extraction, recursive = TRUE), add = TRUE)
  message("[Marine M1] Verifying frozen inputs (no output-oracle reads).")
  inputs <- .lt_marine_inputs(data_root, profile, extraction)
  # Resolve the output parent before creating the run directory, including
  # symlink parents. Frozen source data is never used as an output location.
  parent <- dirname(output_dir)
  if (!dir.exists(parent)) .lt_marine_abort("Output parent does not exist. Create it or choose an existing parent.", "io")
  parent <- normalizePath(parent, mustWork = TRUE)
  output_dir <- file.path(parent, basename(output_dir))
  if (identical(parent, data_root) || startsWith(parent, paste0(data_root, "/"))) {
    .lt_marine_abort("Outputs must be outside the frozen data root.", "output_location")
  }
  if (!dir.create(output_dir)) .lt_marine_abort("Cannot reserve new output directory.", "io")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  run_id <- basename(output_dir)
  stage <- "input_parse"
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  records <- list()
  record <- function(id, detail) {
    records[[length(records) + 1L]] <<- data.frame(stage = id, status = "COMPUTED",
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), detail = detail)
    .lt_marine_write_tsv(do.call(rbind, records), file.path(output_dir, "stages.tsv"))
  }
  tryCatch({
    message("[Marine M1] Reading raw payload, historical partial and supplemental complete state provenance.")
    gbi_intake <- .lt_marine_gbi_intake(inputs$matrix, inputs$crosswalk)
    raw <- gbi_intake$raw
    classified <- .lt_marine_parse_matrix(inputs$classified)
    complete <- .lt_marine_parse_complete_state(inputs$complete_state)
    token_ledger <- .lt_marine_state_layers(raw, complete, classified)
    token_ledger$authority_addendum <- profile$authority_addendum
    old_keys <- .lt_marine_read_tsv(inputs$old_keys)
    tree <- ape::read.tree(inputs$tree)
    traits <- .lt_marine_read_tsv(inputs$traits)
    join <- .lt_marine_branch_join(tree, old_keys, raw$branch_ids)
    if (length(raw$gene_ids) != 17432L || length(raw$branch_ids) != 601L ||
        length(tree$tip.label) != 302L) {
      .lt_marine_abort("Input dimensions do not match the frozen S2 profile.", "axis")
    }
    saveRDS(token_ledger, file.path(output_dir, "upstream_tokens.rds"), version = 3)
    saveRDS(raw$values, file.path(output_dir, "raw_matrix.rds"), version = 3)
    .lt_marine_write_tsv(join, file.path(output_dir, "branch_ledger.tsv"))
    .lt_marine_write_tsv(data.frame(order = seq_along(raw$gene_ids), gene = raw$gene_ids),
                          file.path(output_dir, "gene_ledger.tsv"))
    .lt_marine_write_tsv(traits, file.path(output_dir, "terminal_traits.tsv"))
    rm(classified, complete, token_ledger)
    record(stage, "Raw values unchanged; partial provenance retains its own axis; approved A2U states preserve full axis. Old labels not converted.")

    stage <- "S2_baseline_GBI"
    message("[Marine M1] Computing frozen S2 trimming, arithmetic marginals and GBI.")
    baseline <- .lt_marine_historical_baseline(gbi_intake, output_dir)
    rm(gbi_intake)
    saveRDS(baseline, file.path(output_dir, "S2_baseline_GBI.rds"), version = 3)
    .lt_marine_write_tsv(baseline$gene_effect, file.path(output_dir, "gene_effect.tsv"))
    .lt_marine_write_tsv(baseline$branch_effect, file.path(output_dir, "branch_effect.tsv"))
    counts <- data.frame(metric = c("genes", "branches", "cells", "raw_finite", "raw_zero",
      "trimmed_cells", "trimmed_finite", "ratio_available"),
      value = c(nrow(raw$values), ncol(raw$values), length(raw$values),
        sum(is.finite(raw$values)), sum(raw$values == 0, na.rm = TRUE),
        sum(baseline$trim_mask), sum(is.finite(baseline$trimmed_values)),
        sum(is.finite(baseline$gbi))))
    .lt_marine_write_tsv(counts, file.path(output_dir, "matrix_counts.tsv"))
    record(stage, "Unchanged S2 in admitted native computational order; original writer, text-only old-axis restoration and read-back; all downstream consumers use generated/read-back GBI; no C3.")

    outputs <- list()
    for (trait_id in names(profile$traits)) {
      stage <- paste0("annotation_", trait_id)
      message("[Marine M1] Deterministic ASR: ", trait_id)
      bs <- .lt_marine_annotate(tree, traits, join, profile$traits[[trait_id]], trait_id)
      .lt_marine_write_tsv(bs, file.path(output_dir, paste0(trait_id, ".branch_states.tsv")))
      record(stage, "Frozen numeric trait+1 ASR; first tied column; descendant-node branch join.")
      stage <- paste0("screen_", trait_id)
      message("[Marine M1] Frozen Welch screen: ", trait_id)
      screen <- .lt_marine_screen(baseline$gbi, bs)
      prefix <- file.path(output_dir, trait_id)
      for (key in names(screen)) {
        .lt_marine_write_tsv(screen[[key]], paste0(prefix, ".", key, ".tsv"))
      }
      outputs[[trait_id]] <- screen$summary
      record(stage, paste0(nrow(screen$tested), " tested; ", nrow(screen$significant),
                           " selected by frozen strict-prefix rule."))
    }
    stage <- "provenance"
    profile_file <- system.file("profiles", "marine_S2_frozen_m1.yaml", package = "LaTerra")
    if (!file.copy(profile_file, file.path(output_dir, "frozen_profile.yaml"))) {
      .lt_marine_abort("Could not retain the executed profile.", "io")
    }
    prov <- lt_empty_provenance()
    prov$input_sha256 <- inputs$hash_ledger
    # Extraction locations are disposable; retain an actionable archive/member URI.
    member_rows <- startsWith(prov$input_sha256$path, paste0(extraction, "/"))
    prov$input_sha256$path[member_rows] <- paste0(inputs$branch_archive, "::",
      substring(prov$input_sha256$path[member_rows], nchar(extraction) + 2L))
    prov$configuration_sha256 <- .lt_marine_sha(profile_file)
    prov$software_version <- list(package = "LaTerra", version =
      as.character(utils::packageVersion("LaTerra")),
      certified_base_version = profile$certified_source_base_version,
      certified_base_tar_sha256 = profile$certified_source_base_sha256,
      implementation_certification = "NOT_CLAIMED")
    # Hash the code actually installed and loaded, independently of mutable
    # version strings. The review manifest separately pins the source tarball.
    software_files <- c("DESCRIPTION", "NAMESPACE", "R/LaTerra",
      "R/LaTerra.rdb", "R/LaTerra.rdx")
    software_paths <- vapply(software_files, function(p)
      system.file(p, package = "LaTerra", mustWork = TRUE), character(1L))
    software_hashes <- data.frame(path = software_files,
      sha256 = unname(vapply(software_paths, .lt_marine_sha, character(1L))))
    prov$software_version$installed_runtime_sha256 <- digest::digest(
      serialize(software_hashes, NULL, version = 3), algo = "sha256", serialize = FALSE)
    prov$software_version$installed_runtime_files <- software_hashes
    .lt_marine_write_tsv(software_hashes, file.path(output_dir, "software_hashes.tsv"))
    prov$environment <- env
    prov$ordered_gene_ledger <- data.frame(order = seq_along(raw$gene_ids), gene_id = raw$gene_ids)
    prov$ordered_branch_ledger <- data.frame(order = seq_along(raw$branch_ids), branch_id = raw$branch_ids)
    prov$taxon_domain_ledger <- data.frame(order = seq_along(tree$tip.label),
      taxon_id = tree$tip.label, eligible = TRUE, exclusion_reason = NA_character_)
    prov$eligibility_masks <- list(
      raw_availability = "raw_matrix.rds finite cells",
      original_literal_states = "upstream_tokens.rds",
      trim_mask = "S2_baseline_GBI.rds::trim_mask",
      ratio_availability = "S2_baseline_GBI.rds::gbi finite cells",
      trait_domains = "*.branch_states.tsv::screened",
      scientific_rule = "Frozen S2 domains only; no new generic Layer-2 assignment")
    prov$specification_ledger <- data.frame(specification_id = c("marine_S2_frozen_M1",
      profile$authority_addendum$id, profile$execution_path_addendum$id),
      kind = c("study_specific_reproduction", "supplemental_state_provenance_only",
        "historical_GBI_execution_order_and_generated_consumption_boundary"),
      version = "TARGET001_with_approved_addendum",
      definition = c("frozen_profile.yaml", "frozen_profile.yaml", "GBI_consumption_boundary.rds"),
      sha256 = c(prov$configuration_sha256, prov$configuration_sha256,
        .lt_marine_sha(file.path(output_dir, "GBI_consumption_boundary.rds"))))
    prov$commands <- c(paste0("lt_run_marine(data_root=", encodeString(data_root, quote = "\""),
      ", output_dir=", encodeString(output_dir, quote = "\""),
      ", environment=", encodeString(environment, quote = "\""), ")"),
      paste(commandArgs(), collapse = " "),
      "Deterministic ASR; caller RNG state preserved, including absent seed. castor/Rcpp may initialize RNG bookkeeping; no stochastic choices, new worlds or folds.")
    .lt_marine_write_tsv(prov$input_sha256, file.path(output_dir, "input_hashes.tsv"))
    writeLines(env$session_info, file.path(output_dir, "sessionInfo.txt"))
    record(stage, "Input/member, axes, state, specification, software, environment and output evidence retained.")
    yaml::write_yaml(list(status = "COMPUTED_M1_NOT_REPLAY_CERTIFIED", run_id = run_id,
      started = started, completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      target_scope = "TGT-001 through TGT-006; validation separate",
      full_V1_certification = FALSE, RNG_used = FALSE,
      environment_mode = environment), file.path(output_dir, "RUN_STATUS.yaml"))
    files <- list.files(output_dir, full.names = TRUE)
    prov$output_hashes <- data.frame(role = basename(files), path = basename(files),
      sha256 = vapply(files, .lt_marine_sha, character(1L)))
    run <- lt_run(run_id = run_id, provenance = prov, status = "completed",
      results = list(profile = "marine_S2_frozen", milestone = "M1", output_dir = output_dir,
        screens = outputs, replay_certification = "NOT_CLAIMED"))
    saveRDS(run, file.path(output_dir, "run.rds"), version = 3)
    files <- list.files(output_dir, full.names = TRUE)
    writeLines(paste(vapply(files, .lt_marine_sha, character(1L)), basename(files), sep = "  "),
      file.path(output_dir, "SHA256SUMS"))
    message("[Marine M1] Complete: ", output_dir, ". Replay validation is a separate step.")
    run
  }, error = function(e) {
    yaml::write_yaml(list(status = "FAILED", stage = stage,
      condition_class = class(e), message = conditionMessage(e),
      completed_outputs_are_not_certified = TRUE), file.path(output_dir, "RUN_STATUS.yaml"))
    stop(e)
  })
}
