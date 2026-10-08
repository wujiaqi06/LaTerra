# Resolve only the runtime belonging to this loaded namespace. A different
# installed copy must not mask an unsupported source/load_all session.
.lt_s2_installed_runtime <- function() {
  root <- getNamespaceInfo(environment(.lt_s2_installed_runtime), "path")
  runtime <- file.path(root, "R", "LaTerra.rdb")
  if (!file.exists(runtime) || dir.exists(runtime) || !file.exists(file.path(root, "Meta", "package.rds")))
    stop("The S2 run requires an installed LaTerra source package with its R/LaTerra.rdb runtime. pkgload::load_all() is not supported for audited runs. Install the source package into your chosen library and load it in a fresh R session before running; no output directory was created.", call. = FALSE)
  runtime
}

.lt_s2_user_coordinates <- function(x, tree, coordinates) {
  validate_lt_matrix(x)
  if (!identical(x$payload$type, "branch_length") || any(x$values < 0, na.rm = TRUE))
    stop("The explicit S2 baseline requires nonnegative branch-length payloads, not a precomputed rate.", call. = FALSE)
  if (!inherits(tree, "phylo") || !is.matrix(tree$edge) || ncol(tree$edge) != 2L)
    stop("Supply a phylo tree for annotation on existing coordinates.", call. = FALSE)
  .lt_assert_unique_ids(tree$tip.label, "tree tip keys")
  n_tip <- length(tree$tip.label)
  n_node <- tree$Nnode
  if (!is.numeric(n_node) || length(n_node) != 1L || !is.finite(n_node) || n_node < 1L || n_node != trunc(n_node) ||
      !is.numeric(tree$edge) || anyNA(tree$edge) || any(tree$edge != trunc(tree$edge)) ||
      any(tree$edge < 1L) || any(tree$edge > n_tip + n_node) ||
      nrow(tree$edge) != n_tip + n_node - 1L || anyDuplicated(tree$edge[, 2]) || any(tree$edge[, 1] <= n_tip))
    stop("Malformed annotation-tree node/edge structure.", call. = FALSE)
  roots <- setdiff(unique(tree$edge[, 1]), tree$edge[, 2])
  if (length(roots) != 1L || !setequal(as.vector(tree$edge), seq_len(n_tip + n_node)))
    stop("Annotation tree must have one root and a complete connected node domain.", call. = FALSE)
  seen <- integer(); frontier <- roots
  while (length(frontier)) {
    if (any(frontier %in% seen) || anyDuplicated(frontier)) stop("Cycle in annotation tree.", call. = FALSE)
    seen <- c(seen, frontier)
    frontier <- tree$edge[tree$edge[, 1] %in% frontier, 2]
  }
  if (length(seen) != n_tip + n_node) stop("Disconnected or cyclic annotation tree.", call. = FALSE)
  fields <- c("branch_id", "canonical_split_key", "side_A_taxa", "side_B_taxa", "branch_type", "terminal_taxon")
  if (!is.data.frame(coordinates) || !all(fields %in% names(coordinates)))
    stop(paste("Coordinate ledger requires:", paste(fields, collapse = ", ")), call. = FALSE)
  .lt_assert_unique_ids(coordinates$branch_id, "coordinate branch keys")
  .lt_assert_unique_ids(coordinates$canonical_split_key, "coordinate split keys")
  if (!setequal(coordinates$branch_id, x$branch_ids))
    stop("Matrix and supplied coordinate ledger branch keys disagree.", call. = FALSE)
  coordinates <- coordinates[match(x$branch_ids, coordinates$branch_id), fields, drop = FALSE]
  if (anyNA(coordinates$side_A_taxa) || anyNA(coordinates$side_B_taxa))
    stop("Coordinate split sides must be explicit.", call. = FALSE)
  for (i in seq_len(nrow(coordinates))) {
    a <- strsplit(coordinates$side_A_taxa[i], ";", fixed = TRUE)[[1]]
    b <- strsplit(coordinates$side_B_taxa[i], ";", fixed = TRUE)[[1]]
    if (!length(a) || !length(b) || anyDuplicated(c(a, b)) || !setequal(c(a, b), tree$tip.label))
      stop("Every supplied split must partition the exact tree-tip domain without overlap.", call. = FALSE)
  }
  source <- coordinates
  names(source)[names(source) == "branch_id"] <- "branch_label"
  names(source)[names(source) == "terminal_taxon"] <- "terminal_taxon_if_terminal"
  join <- .lt_marine_branch_join(tree, source, x$branch_ids)
  terminal <- join$offspring <= length(tree$tip.label)
  if (anyNA(join$branch_type) || any(join$branch_type != ifelse(terminal, "terminal", "internal")) ||
      anyNA(join$terminal_taxon[terminal]) ||
      !identical(join$terminal_taxon[terminal], tree$tip.label[join$offspring[terminal]]))
    stop("Supplied terminal/internal identities disagree with the annotation tree.", call. = FALSE)
  join
}

#' Run an explicitly chosen S2 recipe on user-owned data and traits
#'
#' Structural, keyed input validation only: no Marine snapshot SHA, species,
#' gene-count or branch-count requirement. This does not create coordinates or
#' choose a scientific method for the user. The historical trimming/ASR/screen
#' recipe must be explicitly selected; it is not recommended as universal QC.
#'
#' @param x Coordinate-valid branch-length `lt_matrix`.
#' @param tree Supplied `phylo` annotation tree compatible with the existing
#'   coordinate ledger. No alignment, tree inference or coordinate construction.
#' @param coordinates Data frame with branch_id, canonical_split_key,
#'   side_A_taxa, side_B_taxa (semicolon-delimited), branch_type and terminal_taxon.
#' @param trait Explicit binary `lt_trait` with numeric 0/1 and optional excluded
#'   0.5. Continuous traits and automatic recoding are unsupported here.
#' @param terminal_groups Data frame with unique species and explicit genus
#'   grouping keys. Grouping is not inferred from species spelling.
#' @param folds Complete ordered genus/fold_id/seed ledger.
#' @param recipe Required literal `"submitted_S2"`: inclusive q0.975/type7
#'   trimming, arithmetic marginals, GBI, historical ASR and strict-prefix screen.
#' @param validation_recipe Required explicit S2 grouped-model recipe.
#' @param output_dir New output directory; no overwrite.
#' @param cores Explicit positive integer grouped-validation worker count.
#' @return A computed, uncertified `lt_run` and auditable outputs on disk.
#' @export
lt_run_s2 <- function(x, tree, coordinates, trait, terminal_groups, folds,
                       recipe, validation_recipe, output_dir, cores = 1L) {
  if (missing(recipe) || !identical(recipe, "submitted_S2"))
    stop("Explicitly choose recipe = 'submitted_S2'; historical trimming is not a generic QC default.", call. = FALSE)
  if (missing(validation_recipe)) stop("Choose the grouped validation recipe explicitly.", call. = FALSE)
  validation_recipe <- match.arg(validation_recipe, c("nested_phase12B", "nested_phase13", "global_screen_foldwise"))
  runtime_file <- .lt_s2_installed_runtime()
  join <- .lt_s2_user_coordinates(x, tree, coordinates)
  validate_lt_trait(trait)
  if (!identical(trait$type, "binary") || !is.numeric(trait$values) || anyNA(trait$values) ||
      !all(trait$values %in% c(0, 0.5, 1)) || !all(c(0, 1) %in% trait$values))
    stop("Supply an explicit binary 0/1 trait (optional 0.5 excluded), with both endpoint classes; no automatic recoding.", call. = FALSE)
  if (!setequal(trait$taxon_ids, tree$tip.label)) stop("Trait and annotation tree taxon keys disagree.", call. = FALSE)
  coding <- c(positive = 1, negative = 0, intermediate = 0.5, excluded = 0.5)
  for (key in intersect(names(coding), names(trait$coding))) {
    declared <- trait$coding[[key]]
    if (!is.numeric(declared) || length(declared) != 1L || is.na(declared) || declared != coding[[key]])
      stop("Explicit trait coding conflicts with the selected 0/0.5/1 S2 recipe; recode deliberately before running.", call. = FALSE)
  }
  if (!is.data.frame(terminal_groups) || !all(c("species", "genus") %in% names(terminal_groups)))
    stop("terminal_groups requires species and explicit genus keys.", call. = FALSE)
  .lt_assert_unique_ids(terminal_groups$species, "terminal group species")
  if (!setequal(terminal_groups$species, tree$tip.label)) stop("Grouping and trait/tree taxon keys disagree.", call. = FALSE)
  for (package in c("ape", "castor", "glmnet")) if (!requireNamespace(package, quietly = TRUE))
    stop(paste("Install the required package", package), call. = FALSE)
  trait_table <- data.frame(species = trait$taxon_ids, supplied_state = trait$values)
  # Same explicitly selected historical ASR; no vocabulary inference from ID.
  states <- .lt_marine_annotate(tree, trait_table, join, "supplied_state", "explicit_user_binary")
  states$trait <- trait$trait_id; states$run_id <- trait$trait_id
  terminal <- join[join$branch_type == "terminal", , drop = FALSE]
  terminals <- data.frame(species = terminal$terminal_taxon, branch = terminal$branch_id,
    genus = terminal_groups$genus[match(terminal$terminal_taxon, terminal_groups$species)],
    response = trait$values[match(terminal$terminal_taxon, trait$taxon_ids)], stringsAsFactors = FALSE)
  # Validate folds/keys before output creation or baseline/model work.
  .lt_s2_model_inputs(x$values, stats::setNames(states$state, states$branch_id), terminals, folds)
  .lt_assert_scalar_character(output_dir, "output_dir")
  if (file.exists(output_dir)) stop("output_dir exists; choose a new run directory.", call. = FALSE)
  parent <- normalizePath(dirname(output_dir), mustWork = TRUE)
  output_dir <- file.path(parent, basename(output_dir))
  if (!dir.create(output_dir)) stop("Cannot create the new output directory.", call. = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  stage <- "baseline"
  write <- function(value, name) .lt_marine_write_tsv(value, file.path(output_dir, name))
  save <- function(value, name) saveRDS(value, file.path(output_dir, name), version = 3)
  tryCatch({
    supplied <- list(matrix = x, tree = tree, coordinates = coordinates, trait = trait,
      terminal_groups = terminal_groups, folds = folds)
    save(supplied, "supplied_inputs.rds")
    configuration <- list(recipe = recipe, validation_recipe = validation_recipe, cores = cores,
      trim = "q0.975_type7_inclusive", marginals = "arithmetic", operation_order = "divide_BE_then_GE",
      ASR = "castor_numeric_plus_one_first_tie_descendant_node", screen = "Welch_background_minus_focal_strict_prefix_0.01",
      input_role = "user_declared_not_Marine_snapshot")
    save(configuration, "configuration.rds")
    baseline <- .lt_marine_baseline(list(values = x$values, gene_ids = x$gene_ids, branch_ids = x$branch_ids))
    save(baseline, "S2_baseline_GBI.rds")
    write(join, "branch_ledger.tsv"); write(states, "branch_states.tsv")
    write(baseline$gene_effect, "gene_effect.tsv"); write(baseline$branch_effect, "branch_effect.tsv")
    stage <- "screen"
    screen <- .lt_marine_screen(baseline$gbi, states)
    for (field in names(screen)) write(screen[[field]], paste0("screen_", field, ".tsv"))
    stage <- "grouped_validation"
    validation <- lt_s2_gloocv(baseline$gbi, stats::setNames(states$state, states$branch_id), terminals,
      folds, validation_recipe, global_features = screen$significant$gene, cores = cores)
    save(validation, "validation.rds")
    predictions <- do.call(rbind, lapply(validation$folds, function(f) {
      if (!length(f$prediction)) return(NULL)
      data.frame(species = f$identity$test_species, fold_id = f$identity$fold_id,
        genus = f$identity$genus, response = f$identity$y_test, probability = f$prediction, status = f$status)
    }))
    write(predictions, "OOF_predictions.tsv")
    write(data.frame(trait_id = trait$trait_id, n_evaluated = nrow(predictions),
      AUC = .lt_s2_auc(predictions$response, predictions$probability),
      evidence = "grouped_held_out_validation"), "validation_summary.tsv")
    stage <- "provenance"
    provenance <- lt_empty_provenance()
    provenance$input_sha256 <- data.frame(role = "actual_user_supplied_inputs", path = "supplied_inputs.rds",
      sha256 = .lt_marine_sha(file.path(output_dir, "supplied_inputs.rds")))
    provenance$configuration_sha256 <- .lt_marine_sha(file.path(output_dir, "configuration.rds"))
    provenance$software_version <- list(package = "LaTerra", version = as.character(utils::packageVersion("LaTerra")),
      installed_runtime_sha256 = .lt_marine_sha(runtime_file))
    provenance$environment <- list(R = as.character(getRversion()), platform = R.version$platform,
      packages = lapply(c("ape", "castor", "glmnet"), function(p) list(package = p, version = as.character(utils::packageVersion(p)))))
    provenance$ordered_gene_ledger <- data.frame(order = seq_along(x$gene_ids), gene_id = x$gene_ids)
    provenance$ordered_branch_ledger <- data.frame(order = seq_along(x$branch_ids), branch_id = x$branch_ids)
    provenance$taxon_domain_ledger <- data.frame(order = seq_len(nrow(terminals)), taxon_id = terminals$species,
      eligible = terminals$response != 0.5, exclusion_reason = ifelse(terminals$response == 0.5, "explicit_0.5_excluded", NA_character_))
    p <- .lt_s2_fold_provenance(validation, trait$trait_id)
    provenance$fold_ledger <- p$folds; provenance$seed_ledger <- p$seeds
    provenance$eligibility_masks <- list(input = "supplied_inputs.rds::matrix coordinate/value-reason layers unchanged",
      trim = "S2_baseline_GBI.rds::trim_mask", ratio_available = "finite S2_baseline_GBI.rds::gbi",
      annotation = "branch_states.tsv", training_test_excluded = "fold_ledger.tsv")
    provenance$specification_ledger <- data.frame(specification_id = recipe, kind = "explicit_user_selected_historical_recipe",
      version = "S2", definition = "configuration.rds", sha256 = provenance$configuration_sha256)
    provenance$commands <- c(paste(deparse(match.call()), collapse = " "), paste(commandArgs(), collapse = " "))
    write(provenance$fold_ledger, "fold_ledger.tsv"); write(provenance$seed_ledger, "seed_ledger.tsv")
    write(provenance$input_sha256, "input_hashes.tsv")
    writeLines(utils::capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
    yaml::write_yaml(list(status = "COMPUTED_USER_S2_NOT_CERTIFIED", trait = trait$trait_id,
      scientific_default_selected_for_user = FALSE, null_generated = FALSE), file.path(output_dir, "RUN_STATUS.yaml"))
    files <- list.files(output_dir, full.names = TRUE)
    provenance$output_hashes <- data.frame(role = basename(files), path = basename(files),
      sha256 = vapply(files, .lt_marine_sha, character(1)))
    run <- lt_run(basename(output_dir), trait = trait, provenance = provenance, status = "completed",
      results = list(workflow = "user_declared_S2", output_dir = output_dir, screen_count = nrow(screen$significant),
        certification = "NOT_CLAIMED"))
    save(run, "run.rds")
    files <- list.files(output_dir, full.names = TRUE)
    writeLines(paste(vapply(files, .lt_marine_sha, character(1)), basename(files), sep = "  "), file.path(output_dir, "SHA256SUMS"))
    run
  }, error = function(e) {
    yaml::write_yaml(list(status = "FAILED", stage = stage, message = conditionMessage(e),
      partial_outputs_preserved = TRUE), file.path(output_dir, "RUN_STATUS.yaml"))
    stop(e)
  })
}
