# ADDENDUM005: historical screening labels are NOT upstream coordinate truth.
# Canonical sources and complete vectors, never caller-editable attestations.
.lt_marine_phase13_sources <- function() {
  stems <- c("endpointfix_stage2_deterministic_manifest.tsv", "mammal302.anno.nwk",
    "mammal302.anno.BL_support.nwk", "mammal.branch.txt")
  hashes <- c("14fe9c3ace27f88eefbcc7ff61d80aea5e9024b888fcac2f9fbc49c21edbf4fd",
    "eb07d58b0ef08194523570639ed23b6a6f3019e91a3950ca9797a8d59f0978cf",
    "b032a99796f146816e93d3efd585f21fe908e9fabe57366e8c6ac24f410f287e",
    "40caeb6ad333214ed6cc4945b1cefdc01b5332aedbf143434ec239caabfa5065")
  stats::setNames(lapply(seq_along(stems), function(i) list(
    path = file.path("supplemental_authority/phase13", stems[i]), sha256 = hashes[i],
    role = "historical_Phase13_stage2_screening_source_only")),
    c("phase13_manifest", "phase13_main", "phase13_support", "phase13_branches"))
}

.lt_marine_phase13_contract <- function() {
  data.frame(run_id = c("fix_drop_whale", "fix_whale_only", "fix_pinniped_only",
    "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia", "fix_aquatic_v2_noMarineEdge"),
    ordinal = c(2L, 3L, 4L, 6L, 7L, 8L),
    trait_column = c("drop_whale", "whale_only", "pinniped_only", "aquatic_v2_noCetacea",
      "aquatic_v2_noPinnipedia", "aquatic_v2_noMarineEdge"),
    trait_input = c(rep("drop_marine", 3L), rep("drop_aquatic", 3L)),
    literal_sha256 = c("9cf3e6b99b37b6528265614beaf9e524aaf6d3e94305e2f4ac76ebfb1e10d3a2",
      "ac9e715fb76e2b1523d04d616e0a92bae24726910b94a0152126cb431f047b38",
      "67d1d6a37824e6f007d1ba31b14fee8223c9579c1028144722758468eeccbd0b",
      "697f5503644d8c5877474a76a8af8beddc74fba07a98813e26e612ffde2ee863",
      "58a33ff65f88eff05fc6f31962ebc13103febed54c2a42f92eddf3f97bce261e",
      "60d87ac050cd1dd6d0b9ea47e4a22d41c0f35d86a1048904e77b9b3f816a903a"),
    physical_sha256 = c("10bf1ffe9045e80cdb877bf9ee0bd44d2c8a56df8ee8f1c6382d7644d47c0c31",
      "1a2f3d6126726a6bd0d5a3c03edddc7460d83c0a58caa4bced01a36b67e4d0a3",
      "75feb41d01ae51e43e1e0330c21acd322a8aa528d588d2e2600eafc2990df553",
      "2f458bfdb007d09c29f93554f646fbbd7e622001fbee5d81eaafdeab144df623",
      "f509b4c63e5d7d23326832911ed76a610f12072ba6578dce8a32ab299c85af40",
      "601b8050db22252b8916f1657c68da2b087ccaaeeabee6914c4655d2980a3a66"),
    branch_disagreements = c(72L, 18L, 20L, 111L, 88L, 81L),
    terminal_disagreements = c(37L, 6L, 10L, 66L, 53L, 48L), stringsAsFactors = FALSE)
}

.lt_marine_phase13_state_hash <- function(x) {
  digest::digest(paste(paste(names(x), x, sep = "="), collapse = "\n"),
    algo = "sha256", serialize = FALSE)
}

.lt_marine_phase13_materialize <- function(data_root, run_id) {
  .lt_assert_scalar_character(run_id, "historical run_id")
  contract <- .lt_marine_phase13_contract()
  route <- contract[contract$run_id == run_id, , drop = FALSE]
  if (nrow(route) != 1L) .lt_marine_abort("Run is outside the six admitted Phase13 stage2 routes.", "phase13_authority")
  spec <- .lt_marine_m2_authority()$inputs
  spec <- spec[c(names(.lt_marine_phase13_sources()), route$trait_input)]
  paths <- lapply(spec, function(x) {
    p <- file.path(data_root, x$path)
    .lt_marine_verify(p, x$sha256, "canonical Phase13 source")
    normalizePath(p, mustWork = TRUE)
  })
  manifest <- .lt_marine_read_tsv(paths$phase13_manifest)
  mr <- manifest[manifest$run_id == run_id, , drop = FALSE]
  if (nrow(mr) != 1L || !identical(mr$trait_column, route$trait_column))
    .lt_marine_abort("Historical manifest route differs from canonical admission.", "phase13_authority")
  # The manifest's absolute paths are inert provenance, never file-open targets.
  main <- ape::read.tree(paths$phase13_main)
  support <- ape::read.tree(paths$phase13_support)
  parser_hash <- function(t) digest::digest(paste(c(
    paste(seq_along(t$tip.label), t$tip.label, sep = "="),
    paste(seq_len(nrow(t$edge)), t$edge[, 1L], t$edge[, 2L], sep = ":")), collapse = "\n"),
    algo = "sha256", serialize = FALSE)
  ph <- c(main = parser_hash(main), support = parser_hash(support))
  if (!all(ph == "e8bf7d0302939a181e037bc247ba49c84ba9e822c8ddb9befccff86d6a8f5684"))
    .lt_marine_abort("Historical tree parser node/edge identity changed.", "phase13_authority")
  traits <- .lt_marine_read_tsv(paths[[route$trait_input]])
  branch <- .lt_marine_read_tsv(paths$phase13_branches)
  y <- traits[[route$trait_column]][match(support$tip.label, traits$species)]
  anc <- .lt_marine_asr(support, y)
  node_states <- c(y, max.col(anc$ancestral_likelihoods, ties.method = "first") - 1)
  edge_key <- function(t) paste(t$edge[, 1L], t$edge[, 2L], sep = "-")
  states <- node_states[support$edge[match(edge_key(main), edge_key(support)), 2L]]
  names(states) <- paste0("B", seq_len(nrow(main$edge)))
  states <- states[branch$Branch]
  if (!identical(names(states), paste0("B", 1:601)) || anyNA(states) ||
      !identical(.lt_marine_phase13_state_hash(states), route$literal_sha256))
    .lt_marine_abort("Complete literal Phase13 screening vector differs from canonical authority.", "phase13_authority")
  terminals <- branch[branch$Species != "internal", , drop = FALSE]
  terminal <- data.frame(species = terminals$Species, branch = terminals$Branch,
    genus = sub("_.*$", "", terminals$Species),
    response = as.numeric(traits[[route$trait_column]][match(terminals$Species, traits$species)]),
    stringsAsFactors = FALSE)
  if (nrow(terminal) != 302L || anyNA(terminal))
    .lt_marine_abort("Incomplete historical terminal-response domain.", "phase13_authority")
  list(run_id = run_id, ordinal = route$ordinal, trait_column = route$trait_column,
    state = states, terminal = terminal, parser_sha256 = ph,
    source_ledger = do.call(rbind, lapply(names(spec), function(id) data.frame(
      source_id = id, bundle_relative_path = spec[[id]]$path,
      actual_path = paths[[id]], sha256 = spec[[id]]$sha256))),
    historical_manifest_row = mr,
    authority_id = "TARGET001-AUTHORITY-ADDENDUM005-PHASE13-SCREENING-MAP001",
    authority_sha256 = "c5d9bf5ad9b17720413005caca83fa9e7f2028aa6d1f180191d7a9d2a7d34236",
    qualifier = "historical_phase13_screening_mapping_mismatch_preserved")
}

.lt_marine_phase13_validate <- function(x, data_root) {
  if (!is.list(x) || is.null(x$run_id))
    .lt_marine_abort("Missing source-bound historical context.", "phase13_authority")
  expected <- .lt_marine_phase13_materialize(data_root, x$run_id)
  # No caller-updated hash, authority string or coordinated vector edit can
  # replace the exact source-derived object. This restriction is historical-only.
  if (!identical(x, expected))
    .lt_marine_abort("Historical Phase13 context differs from canonical source materialization.", "phase13_authority")
  invisible(x)
}

.lt_marine_phase13_bound_inputs <- function(gbi, model, index, data_root) {
  x <- .lt_marine_phase13_materialize(data_root, model$run_id)
  route <- .lt_marine_phase13_contract()
  route <- route[route$run_id == model$run_id, , drop = FALSE]
  coherent <- stats::setNames(model$state$state, model$state$branch_id)
  if (!is.numeric(index) || length(index) != 1L || is.na(index) || index != x$ordinal ||
      !identical(model$trait_column, x$trait_column) ||
      !identical(model$terminal, x$terminal) ||
      !identical(names(coherent), names(x$state)) ||
      !identical(.lt_marine_phase13_state_hash(coherent), route$physical_sha256))
    .lt_marine_abort("Phase13 run ordinal, terminal response or physical state dependency changed.", "phase13_authority")
  folds <- .lt_marine_model_folds(x$terminal, "nested_phase13", x$ordinal)
  .lt_s2_model_inputs(gbi, coherent, x$terminal, folds)
  x$screening_ledger <- data.frame(branch_id = names(x$state),
    literal_screen_state = unname(x$state), physical_state = unname(coherent),
    mismatch = unname(x$state != coherent))
  x$terminal_ledger <- transform(x$terminal,
    literal_screen_state = unname(x$state[x$terminal$branch]),
    mismatch = unname(x$state[x$terminal$branch] != x$terminal$response))
  if (sum(x$screening_ledger$mismatch) != route$branch_disagreements ||
      sum(x$terminal_ledger$mismatch) != route$terminal_disagreements)
    .lt_marine_abort("Historical discrepancy ledger changed.", "phase13_authority")
  list(context = x, folds = folds)
}

.lt_marine_phase13_gloocv <- function(gbi, model, index, data_root, cores, progress = NULL) {
  checked <- .lt_marine_phase13_bound_inputs(gbi, model, index, data_root)
  x <- checked$context
  result <- .lt_s2_gloocv_execute(gbi, x$state, x$terminal, checked$folds,
    "nested_phase13", model$global_screen$significant$gene, progress, cores,
    paste(deparse(match.call()), collapse = " "))
  result$provenance$input_role <- "canonical_historical_Phase13_stage2_replay_only"
  result$provenance$historical_screening <- x
  result
}
