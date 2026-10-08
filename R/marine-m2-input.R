# Trusted compiled profile identities, not caller-supplied audit specifications.
# lt_marine_profile() returns inspection copies; no public run accepts those
# copies as replacement authority. Other datasets use structural operators.
.lt_marine_m2_authority <- function() {
  list(profile_id = "marine_NC_submitted", milestone = "M2_development",
    execution_status = "ADMITTED_PHASE13_SOURCE_BOUND_REPLAY",
    addenda = c(
      TARGET001_AUTHORITY_ADDENDUM002_FIG5B_LEGACY001 = "aa0cdd0522876d8ca889c28e2bf06013cd3e61d0d0ae379523221a1c9515a165",
      TARGET001_AUTHORITY_ADDENDUM003_TURNOVER_MODULE001 = "04b55295b7f34b98864e25a5f138b440d45688007f1bcfd0683ff37f30df755d",
      TARGET001_AUTHORITY_ADDENDUM004_M2_DEPENDENCIES001 = "8cd9931e24a8df84f94d4f8ce587ea264755ef7dfba9e0e67b2b0b731c01623a",
      TARGET001_AUTHORITY_ADDENDUM005_PHASE13_SCREENING_MAP001 = "c5d9bf5ad9b17720413005caca83fa9e7f2028aa6d1f180191d7a9d2a7d34236",
      TARGET001_AUTHORITY_ADDENDUM006_GBI_ORDER_SERIALIZATION001 = "cebaafc6fdfb019d36d9fcb0960f171a765f65cb6a767f06c9a156a5d9ec99c5"),
    GBI_dependency = "fresh_ADDENDUM006_generated_consumed_bits_no_old_fit_or_feature_cache_reuse",
    inputs = c(list(
      drop_marine = list(path = "supplemental_authority/trait_input.marine_binary_drop_sets.combined.tsv",
        sha256 = "bc440fc8008bc82e18f466b907a44e2fbe1caf30c85929d8b4d8c815d6011489", role = "trait_input"),
      drop_aquatic = list(path = "supplemental_authority/trait_input.aquatic_v2_drop_sets.combined.tsv",
        sha256 = "e5463d8ee3f563f97edc160612a45ba89ef65bb49cb0afdccff4e1dffe60fca2", role = "trait_input"),
      turnover = list(path = "supplemental_authority/marine_NC_submitted_turnover_module_lookup_131x4.tsv",
        sha256 = "5a05051d86b9e23c1ea0b3c9d45dd573ce4dc6463ad9933979cf66cd8228d243", role = "annotation_only"),
      legacy_auc = list(path = "supplemental_authority/marine_NC_submitted_legacy_AUC_reference_2x2.tsv",
        sha256 = "9a07bcd4721650f203159341d56b46f9125396262c775768c7119e7fce7b4f7e", role = "historical_reference_imported_not_recomputed"),
      fig5b_worlds = list(path = "10_permutation_controls/Fig5B_positive_count_matched_permutation/endpointfix_permutation_label_sets_long.tsv",
        sha256 = "6858293118ba95dc73cb8a5239a8c232657a197a379a7246daadcec033a0e003", role = "exact_persisted_terminal_worlds"),
      final_modules = list(path = "14_supplementary_tables/TableS5_tsv_exports/Supplementary_Table_S5_Figure4C_predictor_annotation.tsv",
        sha256 = "85b7d121a54cc7e2e86ae2c522823a69492974ff2766584dfe540a2115065ca1", role = "final_display_annotation_not_feature_selection")),
      .lt_marine_phase13_sources()),
    turnover_source_sha256 = "3552b6f3876f4b3cd27c06e6a2f0ca0c70b95f97be332b64ae93b03eaed3a0d1",
    legacy_auc_source_sha256 = "9be31eaa34be030e87ac322895a3bcaaf333bb92c6c9b503e7e32e496deccdb2",
    reference_environment = list(R = "4.4.2", glmnet = "4.1-10"),
    qualifiers = c("historical_direction_mismatch_preserved",
      "historical_reference_imported_not_recomputed",
      "historical_phase13_screening_mapping_mismatch_preserved"))
}

.lt_marine_m2_inputs <- function(data_root) {
  authority <- .lt_marine_m2_authority()
  paths <- lapply(authority$inputs, function(x) {
    path <- file.path(data_root, x$path)
    .lt_marine_verify(path, x$sha256, x$role)
    normalizePath(path, mustWork = TRUE)
  })
  annotation <- .lt_marine_turnover_annotation_input(paths$turnover, list(
    delivery_sha256 = authority$inputs$turnover$sha256,
    full_source_sha256 = authority$turnover_source_sha256,
    authority_id = "TARGET001-AUTHORITY-ADDENDUM003-TURNOVER-MODULE001",
    authority_sha256 = unname(authority$addenda[[2L]])))
  read <- function(p) utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE)
  reference <- .lt_read_keyed_character_fields(paths$legacy_auc,
    c("run_id", "gLOOCV_AUC"), "run_id")
  phase13 <- lapply(.lt_marine_phase13_contract()$run_id,
    function(id) .lt_marine_phase13_materialize(data_root, id))
  names(phase13) <- .lt_marine_phase13_contract()$run_id
  list(drop_marine = read(paths$drop_marine), drop_aquatic = read(paths$drop_aquatic),
    phase13 = phase13,
    turnover = annotation, legacy_auc = reference,
    # Keep the reference out of every model/screen argument; reporting only.
    legacy_reference_provenance = list(
      role = "historical_reference_imported_not_recomputed",
      source_sha256 = authority$legacy_auc_source_sha256,
      view_sha256 = authority$inputs$legacy_auc$sha256,
      rows = reference$run_id, fields = names(reference),
      sole_downstream_field = "AUC_legacy_global_preprocessing"),
    worlds = read(paths$fig5b_worlds),
    final_module_path = paths$final_modules,
    input_ledger = do.call(rbind, lapply(names(paths), function(id) data.frame(
      role = authority$inputs[[id]]$role, input_id = id, path = paths[[id]],
      sha256 = authority$inputs[[id]]$sha256))), authority = authority)
}

.lt_marine_terminal_frame <- function(branch_join, traits, column) {
  rows <- branch_join[branch_join$branch_type == "terminal", , drop = FALSE]
  if (!column %in% names(traits) || !setequal(rows$terminal_taxon, traits$species))
    .lt_marine_abort("Terminal mapping and explicit trait keys disagree.", "trait")
  idx <- match(rows$terminal_taxon, traits$species)
  data.frame(species = rows$terminal_taxon, branch = rows$branch_id,
    genus = vapply(strsplit(rows$terminal_taxon, "_", fixed = TRUE), `[`, character(1), 1L),
    response = as.numeric(traits[[column]][idx]), stringsAsFactors = FALSE)
}

.lt_marine_model_folds <- function(terminals, stage, model_index) {
  stage <- match.arg(stage, c("phase11", "nested_phase12B", "nested_phase13"))
  genera <- sort(unique(terminals$genus))
  offset <- if (stage == "phase11") 0L else model_index * 100000L
  base <- if (stage == "nested_phase13") 20260524L else 20260523L
  data.frame(fold_id = seq_along(genera), genus = genera,
    seed = base + seq_along(genera) + offset, stringsAsFactors = FALSE)
}
