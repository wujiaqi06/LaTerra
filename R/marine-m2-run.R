# Study-specific orchestration; numerical work lives in reusable operators.
.lt_marine_m2_models <- function(gbi, tree, branch_join, inputs) {
  marine_cols <- names(inputs$drop_marine)[-1L]
  aquatic_cols <- names(inputs$drop_aquatic)[-1L]
  ids <- c("fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear", "fix_drop_sea_otter",
    "fix_drop_pinniped", "fix_whale_only", "fix_pinniped_only", "fix_aquatic_v2", paste0("fix_", aquatic_cols[-1L]))
  if (length(ids) != 14L || length(marine_cols) != 7L || length(aquatic_cols) != 7L)
    .lt_marine_abort("Incomplete admitted drop-table schema.", "trait")
  models <- stats::setNames(vector("list", length(ids)), ids)
  for (j in seq_along(ids)) {
    d <- if (j <= 7L) inputs$drop_marine else inputs$drop_aquatic
    column <- if (j <= 7L) marine_cols[j] else aquatic_cols[j - 7L]
    state <- .lt_marine_annotate(tree, d, branch_join, column, sub("^fix_", "", ids[j]))
    models[[j]] <- list(run_id = ids[j], trait_column = column, state = state,
      global_screen = .lt_marine_screen(gbi, state),
      terminal = .lt_marine_terminal_frame(branch_join, d, column))
  }
  models
}

.lt_marine_m2_gloocv <- function(gbi, model, stage, index, cores, progress = NULL, data_root = NULL) {
  if (stage == "nested_phase13" && model$run_id %in% .lt_marine_phase13_contract()$run_id) {
    if (is.null(data_root)) .lt_marine_abort("Historical Phase13 replay requires its declared input bundle.", "phase13_authority")
    return(.lt_marine_phase13_gloocv(gbi, model, index, data_root, cores, progress))
  }
  lt_s2_gloocv(gbi, stats::setNames(model$state$state, model$state$branch_id),
    model$terminal, .lt_marine_model_folds(model$terminal, stage, index),
    recipe = if (stage == "phase11") "global_screen_foldwise" else stage,
    global_features = model$global_screen$significant$gene, cores = cores, progress = progress)
}

.lt_s2_fold_provenance <- function(result, scope) {
  terminal <- result$provenance$terminal_ledger
  fold <- lapply(result$folds, function(f) {
    id <- f$identity
    data.frame(fold_id = paste(scope, id$fold_id, sep = "::"), taxon_id = terminal$species,
      role = ifelse(terminal$response == 0.5, "excluded",
        ifelse(terminal$species %in% id$test_species, "test", "train")),
      group = terminal$genus, stringsAsFactors = FALSE)
  })
  seeds <- result$provenance$fold_ledger
  rng <- result$provenance$RNGkind
  list(folds = do.call(rbind, fold), seeds = data.frame(
    scope = paste(scope, seeds$fold_id, sep = "::"), seed = seeds$seed,
    rng_kind = rng[1], normal_kind = rng[2], sample_kind = rng[3],
    r_version = as.character(getRversion()), stringsAsFactors = FALSE))
}

.lt_marine_module_display <- function(annotation, full) {
  .lt_assert_unique_ids(annotation$gene, "final display gene keys")
  selected <- unique(unlist(lapply(full$models[c("fix_marine_binary", "fix_aquatic_v2")],
    function(x) names(x$beta)[x$beta != 0]), use.names = FALSE))
  if (!setequal(selected, annotation$gene))
    .lt_marine_abort("Computed baseline predictor union and admitted display annotation keys differ.", "annotation_join")
  flag <- annotation$counted_in_Figure4C_circle_size
  if (anyNA(flag) || !all(flag %in% c("True", "False")))
    .lt_marine_abort("Unknown admitted display-count flag.", "annotation_fields")
  counted <- flag == "True"
  list(annotation = annotation, counts = data.frame(metric = c("union_predictors", "counted", "not_counted", "display_modules"),
    value = c(length(selected), sum(counted), sum(!counted), length(unique(annotation$display_module[counted])))),
    provenance = list(role = "admitted_annotation_join_to_computed_predictor_union",
      source_numeric_annotation_columns_used_for_fitting = FALSE))
}

#' Run the complete frozen Marine/S2 M1-M2 workflow
#'
#' Executes the submitted study profile from raw inputs through screens,
#' historical permutation, grouped validation, full-data fits, turnover and
#' descriptive projections. No expected fitted-output tables are read.
#' Scientific/replay certification remains a separate review.
#'
#' @param data_root Frozen input bundle root, including M1 and the explicitly
#'   admitted M2 supplemental inputs. See [lt_marine_profile()].
#' @param output_dir New directory outside data_root; no overwrite or implicit
#'   reuse of an earlier fitted run is permitted.
#' @param environment `"frozen"` requires the reference environment including
#'   glmnet 4.1-10. `"compatibility"` is warned and non-certifying.
#' @param cores Positive integer worker count for ordered grouped-validation
#'   batches. Windows requires 1. Explicit historical seeds are unchanged.
#' @return A completed `lt_run`, output tables, fitted stage objects and full
#'   input/configuration/software/seed/fold/output provenance on disk.
#' @export
lt_run_marine_m2 <- function(data_root, output_dir,
                             environment = c("frozen", "compatibility"), cores = 1L) {
  environment <- match.arg(environment)
  .lt_assert_scalar_character(data_root, "data_root")
  .lt_assert_scalar_character(output_dir, "output_dir")
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) || cores < 1 || cores != trunc(cores))
    stop("cores must be a positive integer.", call. = FALSE)
  if (cores > 1L && .Platform$OS.type == "windows")
    stop("Forked execution is unavailable here; explicitly choose cores = 1.", call. = FALSE)
  # ADDENDUM005 admission and boundary checks resolve the historical-only gate.
  # Unknown authority is still forbidden; no editable inspection copy is used.
  if (!identical(.lt_marine_m2_authority()$execution_status, "ADMITTED_PHASE13_SOURCE_BOUND_REPLAY"))
    .lt_marine_abort("Marine M2 lacks recognized source-bound Phase13 replay authority.", "phase13_authority")
  if (!requireNamespace("glmnet", quietly = TRUE))
    .lt_marine_abort("Marine M2 requires glmnet; install it before running.", "dependency")
  compatible <- utils::compareVersion(as.character(utils::packageVersion("glmnet")), "4.1-10") == 0L
  if (!compatible && environment == "frozen")
    .lt_marine_abort("Frozen M2 requires glmnet 4.1-10; no silent version fallback.", "environment")
  if (!compatible) warning("Explicit compatibility execution with non-reference glmnet; not certified.", call. = FALSE)
  if (!identical(RNGkind(), c("Mersenne-Twister", "Inversion", "Rejection")))
    .lt_marine_abort("Set the explicit frozen RNGkind before M2 execution; no silent RNG-kind replacement.", "rng")
  if (!dir.exists(data_root)) .lt_marine_abort("data_root must be an existing input bundle.", "missing_input")
  data_root <- normalizePath(data_root, mustWork = TRUE)
  profile <- lt_marine_profile()
  .lt_marine_environment(profile, environment)
  inputs <- .lt_marine_m2_inputs(data_root)
  if (file.exists(output_dir)) .lt_marine_abort("output_dir already exists; choose a new run directory.", "output_exists")
  parent <- normalizePath(dirname(output_dir), mustWork = TRUE)
  if (identical(parent, data_root) || startsWith(parent, paste0(data_root, "/")))
    .lt_marine_abort("Outputs must be outside the frozen input root.", "output_location")
  output_dir <- file.path(parent, basename(output_dir))
  if (!dir.create(output_dir)) .lt_marine_abort("Cannot reserve output directory.", "io")
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  started <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  stage <- "M1"; stages <- list(); warning_ledger <- list(); seed_parts <- fold_parts <- list()
  write <- function(x, name) .lt_marine_write_tsv(x, file.path(output_dir, name))
  save <- function(x, name) saveRDS(x, file.path(output_dir, name), version = 3)
  execute <- function(id, f) {
    stage <<- id
    message("[Marine M2] ", id)
    result <- withCallingHandlers(f(), warning = function(w) {
      warning_ledger[[length(warning_ledger) + 1L]] <<- data.frame(stage = id, message = conditionMessage(w))
      # Record and still emit, including historical optimizer warnings.
    })
    stages[[length(stages) + 1L]] <<- data.frame(stage = id, status = "COMPUTED",
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    write(do.call(rbind, stages), "stages.tsv")
    result
  }
  tryCatch({
    write(inputs$input_ledger, "M2_input_hashes.tsv")
    save(inputs$authority, "M2_authority.rds")
    save(inputs$legacy_reference_provenance, "legacy_reference_role.rds")
    m1 <- execute("M1", function() lt_run_marine(data_root, file.path(output_dir, "M1"), environment))
    baseline <- readRDS(file.path(output_dir, "M1/S2_baseline_GBI.rds"))
    gbi <- .lt_marine_consumed_gbi(baseline, file.path(output_dir, "M1"))
    rm(baseline)
    branches <- .lt_marine_read_tsv(file.path(output_dir, "M1/branch_ledger.tsv"))
    traits <- .lt_marine_read_tsv(file.path(output_dir, "M1/terminal_traits.tsv"))
    scratch <- tempfile("laterra-m2-input-"); dir.create(scratch)
    on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
    original <- .lt_marine_inputs(data_root, profile, scratch)
    tree <- ape::read.tree(original$tree)
    models <- execute("fourteen_global_screens", function() .lt_marine_m2_models(gbi, tree, branches, inputs))
    save(models, "computed_model_inputs.rds")
    # Preserve complete historical-vs-physical disclosures before any model fit.
    for (id in names(inputs$phase13)) {
      x <- .lt_marine_phase13_bound_inputs(gbi, models[[id]], inputs$phase13[[id]]$ordinal, data_root)$context
      directory <- file.path(output_dir, "phase13_authority", id)
      dir.create(directory, recursive = TRUE)
      saveRDS(x, file.path(directory, "source_bound_context.rds"), version = 3)
      for (field in c("screening_ledger", "terminal_ledger", "source_ledger"))
        .lt_marine_write_tsv(x[[field]], file.path(directory, paste0(field, ".tsv")))
      writeLines(x$qualifier, file.path(directory, "QUALIFIER.txt"))
    }
    for (id in names(models)) {
      directory <- file.path(output_dir, "screens", id); dir.create(directory, recursive = TRUE)
      .lt_marine_write_tsv(models[[id]]$state, file.path(directory, "branch_states.tsv"))
      for (field in names(models[[id]]$global_screen))
        .lt_marine_write_tsv(models[[id]]$global_screen[[field]], file.path(directory, paste0(field, ".tsv")))
    }
    fig5b <- execute("Fig5B_historical_200_worlds", function() .lt_marine_fig5b(gbi, tree, branches,
      inputs$drop_marine, inputs$worlds, models$fix_drop_whale$global_screen,
      function(i, n, x) if (i %% 20L == 0L) message("  Fig5B ", i, "/", n)))
    save(fig5b, "Fig5B.rds"); write(fig5b$screening, "Fig5B_screening.tsv")
    write(fig5b$observed_vs_null, "Fig5B_observed_vs_null.tsv")
    grouped <- function(stage_id, ids) {
      result <- stats::setNames(vector("list", length(ids)), ids)
      for (j in seq_along(ids)) {
        scope <- paste(stage_id, ids[j], sep = "_")
        result[[j]] <- execute(scope, function() .lt_marine_m2_gloocv(gbi, models[[ids[j]]], stage_id, j, cores,
          function(i, n, x) if (i %% 24L == 0L || i == n) message("  ", scope, " ", i, "/", n),
          data_root = data_root))
        save(result[[j]], paste0(scope, ".rds"))
        p <- .lt_s2_fold_provenance(result[[j]], scope)
        seed_parts[[scope]] <<- p$seeds; fold_parts[[scope]] <<- p$folds
      }
      result
    }
    baseline_ids <- c("fix_marine_binary", "fix_aquatic_v2")
    nested <- grouped("nested_phase12B", baseline_ids)
    phase11 <- grouped("phase11", baseline_ids)
    nested_tables <- execute("nested_output_assembly", function()
      .lt_marine_nested_tables(nested, phase11, models, traits, inputs$legacy_auc))
    save(nested_tables, "nested_tables.rds")
    for (field in c("auc", "oof", "folds", "selected")) write(nested_tables[[field]], paste0("nested_", field, ".tsv"))
    rm(nested, phase11)
    sensitivity <- grouped("nested_phase13", .lt_marine_sensitivity_ids())
    sensitivity_tables <- .lt_marine_sensitivity_tables(sensitivity)
    write(sensitivity_tables$table, "Fig5A_sensitivity.tsv")
    write(sensitivity_tables$unrounded, "Fig5A_unrounded.tsv")
    save(sensitivity_tables$provenance, "Fig5A_provenance.rds")
    rm(sensitivity)
    full <- execute("full_data_fourteen_fits", function() .lt_marine_full_data(gbi, models,
      function(i, n, x) message("  full data ", i, "/", n)))
    save(full, "full_data.rds")
    for (field in c("coefficients", "summary", "partition")) write(full[[field]], paste0("full_", field, ".tsv"))
    annotation <- .lt_marine_read_tsv(inputs$final_module_path, TRUE)
    display <- execute("module_display_join", function() .lt_marine_module_display(annotation, full))
    write(display$annotation, "display_annotation.tsv"); write(display$counts, "display_counts.tsv")
    save(display$provenance, "display_provenance.rds")
    turnover <- execute("turnover_exact_streams", function() .lt_marine_turnover(models, full, inputs$turnover$annotation))
    save(turnover, "turnover.rds")
    for (field in c("observed", "universes", "worlds", "summary")) write(turnover[[field]], paste0("turnover_", field, ".tsv"))
    internal <- execute("terminal_only_internal_projection", function() .lt_marine_internal_projections(gbi, models))
    save(internal, "internal_projections.rds")
    write(internal$projections, "internal_projections.tsv"); write(internal$targets, "internal_targets.tsv")
    nc_oof <- .lt_marine_nc_oof(nested_tables$oof)
    fingerprints <- execute("fitted_fingerprints", function() .lt_marine_fingerprints(gbi, models, full, traits, annotation, nc_oof))
    save(fingerprints, "fingerprints.rds")
    write(nc_oof, "nc_oof.tsv"); write(fingerprints$projections, "focal_projections.tsv")
    write(fingerprints$fingerprints, "all_fingerprints.tsv")
    write(.lt_marine_nc_fingerprints(fingerprints), "nc_fingerprints.tsv")
    write(.lt_marine_nc_projection_profiles(fingerprints$projections, internal$projections), "nc_projection.tsv")
    stage <- "provenance"
    prov <- m1$provenance
    prov$specification_ledger$definition <- paste0("M1/", prov$specification_ledger$definition)
    prov$eligibility_masks <- list(M1 = list(directory = "M1", layers = m1$provenance$eligibility_masks))
    prov$input_sha256 <- rbind(prov$input_sha256, inputs$input_ledger[c("role", "path", "sha256")])
    configuration <- list(M1 = profile, M2 = inputs$authority, environment = environment, cores = cores,
      scope = "frozen_Marine_S2_reproduction", historical_qualifiers = inputs$authority$qualifiers)
    save(configuration, "configuration.rds")
    prov$configuration_sha256 <- .lt_marine_sha(file.path(output_dir, "configuration.rds"))
    prov$environment$glmnet <- as.character(utils::packageVersion("glmnet"))
    prov$environment$RNGkind <- RNGkind()
    prov$seed_ledger <- do.call(rbind, seed_parts)
    extra <- c(stats::setNames(vapply(full$models, `[[`, integer(1), "seed"), paste0("full_data::", names(full$models))),
      stats::setNames(vapply(turnover$provenance, `[[`, integer(1), "seed"), paste0("turnover::", names(turnover$provenance))))
    # Internal projection uses the independent frozen seed formula retained
    # with its fitted model provenance, not any baseline fit substitute.
    extra <- c(extra, stats::setNames(vapply(internal$models, `[[`, integer(1), "seed"),
      c("internal_projection::marine", "internal_projection::aquatic")))
    rng <- RNGkind()
    prov$seed_ledger <- rbind(prov$seed_ledger, data.frame(scope = names(extra), seed = as.integer(extra),
      rng_kind = rng[1], normal_kind = rng[2], sample_kind = rng[3], r_version = as.character(getRversion())))
    prov$fold_ledger <- do.call(rbind, fold_parts)
    prov$eligibility_masks$M2 <- list(traits = "computed_model_inputs.rds::terminal/state",
      fold_domains = "fold_ledger.tsv", world_domains = "Fig5B.rds", preprocessing = "per-stage fitted RDS objects")
    prov$commands <- c(prov$commands, paste0("lt_run_marine_m2(data_root=", encodeString(data_root, quote = '"'),
      ", output_dir=", encodeString(output_dir, quote = '"'), ", environment=", encodeString(environment, quote = '"'),
      ", cores=", cores, ")"), "Fig5B uses persisted worlds; turnover uses the exact historical stream. No new generator.")
    spec <- inputs$authority$addenda
    prov$specification_ledger <- rbind(prov$specification_ledger, data.frame(specification_id = names(spec),
      kind = "supplemental_frozen_M2_authority", version = "TARGET001",
      definition = "M2_authority.rds", sha256 = unname(spec)))
    write(prov$input_sha256, "input_hashes.tsv"); write(prov$seed_ledger, "seed_ledger.tsv")
    write(prov$fold_ledger, "fold_ledger.tsv")
    write(prov$software_version$installed_runtime_files, "software_hashes.tsv")
    writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
    warnings <- if (length(warning_ledger)) do.call(rbind, warning_ledger) else
      data.frame(stage = character(), message = character())
    write(warnings, "warnings.tsv")
    writeLines(c(inputs$authority$qualifiers, "full_data_fit_and_projection_are_not_validation",
      "historical_Task001_Task003_labels_are_inert_submitted_source_labels",
      "No scientific/replay/M2/release acceptance is implied by computational completion."),
      file.path(output_dir, "QUALIFIERS.txt"))
    yaml::write_yaml(list(status = "COMPUTED_M2_NOT_REPLAY_CERTIFIED", started = started,
      completed = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), environment = environment,
      target_scope = "TGT-001 through TGT-020; validation separate", release_certification = FALSE),
      file.path(output_dir, "RUN_STATUS.yaml"))
    files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    relative <- substring(files, nchar(output_dir) + 2L)
    prov$output_hashes <- data.frame(role = relative, path = relative, sha256 = vapply(files, .lt_marine_sha, character(1)))
    run <- lt_run(basename(output_dir), provenance = prov, status = "completed",
      results = list(profile = "marine_NC_submitted", milestone = "M2", output_dir = output_dir,
        replay_certification = "NOT_CLAIMED", stage_count = length(stages)))
    save(run, "run.rds")
    files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    writeLines(paste(vapply(files, .lt_marine_sha, character(1)), substring(files, nchar(output_dir) + 2L), sep = "  "),
      file.path(output_dir, "SHA256SUMS"))
    message("[Marine M2] Computation complete; independent replay review remains required: ", output_dir)
    run
  }, error = function(e) {
    yaml::write_yaml(list(status = "FAILED", stage = stage, message = conditionMessage(e),
      condition_class = class(e), partial_outputs_preserved = TRUE, certification = "NOT_CLAIMED"),
      file.path(output_dir, "RUN_STATUS.yaml"))
    stop(e)
  })
}
