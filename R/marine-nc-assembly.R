.lt_marine_nc_oof <- function(nested_oof) {
  allowed <- c("marine_binary_nested_ttest_baseline", "binary_aquatic_dependence_nested_ttest_baseline")
  if (!all(nested_oof$model %in% allowed)) .lt_marine_abort("Unexpected nested model identity in NC export.", "model_order")
  .lt_assert_unique_ids(paste(nested_oof$species, nested_oof$model, sep = "||"), "NC species/model keys")
  p <- .lt_s2_numeric_tsv_roundtrip(nested_oof$nested_ttest_OOF_prediction)
  clipped <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  # Literal NC source: log(p/(1-p)), before either independent .10g export.
  logit <- log(clipped / (1 - clipped))
  format10 <- function(x) ifelse(is.finite(x), sprintf("%.10g", x), NA_character_)
  data.frame(species = nested_oof$species,
    model = ifelse(nested_oof$model == "marine_binary_nested_ttest_baseline", "marine", "binary_aquatic_dependence"),
    probability = as.numeric(format10(p)), logit = as.numeric(format10(logit)),
    probability_text = format10(p), logit_text = format10(logit), stringsAsFactors = FALSE)
}

.lt_marine_nc_fingerprints <- function(fingerprint_result) {
  # Explicit source-stage write/read boundaries, using R's original CSV writer.
  fp <- fingerprint_result$fingerprints
  proj <- fingerprint_result$projections
  for (field in names(fp)[vapply(fp, is.numeric, logical(1))])
    fp[[field]] <- .lt_s2_numeric_tsv_roundtrip(fp[[field]])
  for (field in c("fitted_probability", "fitted_logit"))
    proj[[field]] <- .lt_s2_numeric_tsv_roundtrip(proj[[field]])
  key <- paste(proj$species, proj$model, sep = "||")
  .lt_assert_unique_ids(key, "fitted projection species/model keys")
  match_id <- match(paste(fp$species, fp$model, sep = "||"), key)
  if (anyNA(match_id)) .lt_marine_abort("Missing computed projection for NC fingerprint metadata join.", "dependency")
  fp$fitted_probability <- proj$fitted_probability[match_id]
  fp$fitted_logit <- proj$fitted_logit[match_id]
  species <- c("Orcinus_orca", "Zalophus_californianus", "Leptonychotes_weddellii",
    "Odobenus_rosmarus_divergens", "Enhydra_lutris_kenyoni", "Ursus_maritimus")
  labels <- c("killer whale", "California sea lion", "Weddell seal", "walrus", "sea otter", "polar bear")
  roles <- c("core cetacean", "core pinniped", "heterogeneous pinniped", "heterogeneous pinniped",
    "marine-edge decoupled", "marine-edge decoupled")
  id <- match(fp$species, species)
  fp$display_label[!is.na(id)] <- labels[id[!is.na(id)]]
  fp$portrait_role <- ifelse(is.na(id), "", roles[id])
  # pandas' historical source-table export spells boolean payloads this way.
  fp$imputed_flag <- ifelse(fp$imputed_flag, "True", "False")
  fp[!is.na(id), , drop = FALSE]
}

.lt_marine_nc_projection_profiles <- function(focal, internal) {
  species <- c("Orcinus_orca", "Zalophus_californianus", "Leptonychotes_weddellii",
    "Odobenus_rosmarus_divergens", "Dugong_dugon", "Ursus_maritimus", "Enhydra_lutris_kenyoni",
    "Platanista_minor", "Inia_geoffrensis", "Lipotes_vexillifer", "Hippopotamus_amphibius",
    "Aonyx_cinereus", "Pteronura_brasiliensis")
  terminal_group <- c("Core cetacean", "Core pinniped", rep("Heterogeneous pinniped", 2),
    "Sirenian / edge", rep("Marine edge", 2), rep("River dolphin bridge", 3), rep("Non-marine aquatic controls", 3))
  terminal_label <- c("Killer whale", "California sea lion", "Weddell seal", "Walrus", "Dugong",
    "Polar bear", "Sea otter", "Platanista", "Inia", "Baiji", "Hippopotamus", "Small-clawed otter", "Giant otter")
  targets <- c("common_ancestor_Cetacea_plus_Hippopotamus", "crown_Cetacea", "Odontoceti", "Pinnipedia",
    "Phocidae", "Otarioidea", "Otariidae", "Sirenia", "Dugong_plus_Hydrodamalis")
  internal_group <- c(rep("Cetacean ancestry", 3), rep("Pinniped ancestry", 4), rep("Sirenian decoupling context", 2))
  internal_label <- c("Cetacea + hippo ancestor", "Crown Cetacea", "Odontoceti", "Crown Pinnipedia",
    "Phocidae", "Otarioidea", "Otariidae", "Sirenia", "Dugong + Hydrodamalis")
  # Original source reads serialized model tables; numerical rounding here is
  # that existing stage boundary, not a newly chosen precision or tolerance.
  focal$fitted_probability <- .lt_s2_numeric_tsv_roundtrip(focal$fitted_probability)
  internal$predicted_probability <- .lt_s2_numeric_tsv_roundtrip(internal$predicted_probability)
  rows <- list()
  for (i in seq_along(species)) {
    d <- focal[focal$species == species[i], , drop = FALSE]
    if (nrow(d) != 2L || !setequal(d$model, c("marine", "binary_aquatic_dependence")))
      .lt_marine_abort("Missing/ambiguous paired terminal projection; no silent missing row.", "dependency")
    rows[[i]] <- data.frame(row_type = "terminal", group = terminal_group[i], label = terminal_label[i],
      source_id = species[i], marine_probability = d$fitted_probability[d$model == "marine"],
      aquatic_probability = d$fitted_probability[d$model == "binary_aquatic_dependence"],
      trait_category = d$trait_category[1], aquaticity_score = d$aquaticity_score[1],
      # Preserve the submitted table's historical task label, not this run ID.
      source_layer = "Task003 corrected full-data terminal fitted projection",
      branch_id = NA_character_, target_type = NA_character_, target_note = NA_character_, stringsAsFactors = FALSE)
  }
  for (i in seq_along(targets)) {
    d <- internal[internal$target_name == targets[i], , drop = FALSE]
    if (nrow(d) != 2L || !setequal(d$trait, c("marine_binary", "binary_aquatic_dependence")))
      .lt_marine_abort("Missing/ambiguous paired internal projection; no silent missing row.", "dependency")
    rows[[length(rows) + 1L]] <- data.frame(row_type = "internal", group = internal_group[i], label = internal_label[i],
      source_id = targets[i], marine_probability = d$predicted_probability[d$trait == "marine_binary"],
      aquatic_probability = d$predicted_probability[d$trait == "binary_aquatic_dependence"],
      trait_category = NA_character_, aquaticity_score = NA_real_,
      source_layer = "Task001 terminal-only internal branch projection; descriptive, not OOF validation",
      branch_id = d$branch_id[1], target_type = d$target_type[1], target_note = d$target_note[1], stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  d$axis_delta_aquatic_minus_marine <- d$aquatic_probability - d$marine_probability
  m <- d$marine_probability; a <- d$aquatic_probability
  d$projection_summary <- ifelse(m >= 0.75 & a >= 0.75, "high on both fitted axes",
    ifelse(m < 0.25 & a >= 0.75, "aquatic-axis high, marine-axis low",
      ifelse(m >= 0.75 & a < 0.25, "marine-axis high, aquatic-axis low",
        ifelse(m < 0.25 & a < 0.25, "low on both fitted axes", "intermediate or divergent"))))
  d
}
