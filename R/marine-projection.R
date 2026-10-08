.lt_marine_internal_targets <- function(branch_ids) {
  targets <- data.frame(
    target_name = c("common_ancestor_Cetacea_plus_Hippopotamus", "crown_Cetacea", "Mysticeti", "Odontoceti",
      "Platanista_terminal_branch", "Inia_plus_Lipotes", "Pinnipedia", "Phocidae", "Otarioidea", "Otariidae",
      "Otariinae_sampled", "Sirenia", "Dugong_plus_Hydrodamalis", "Trichechus_terminal_branch",
      "Ursus_arctos_plus_Ursus_maritimus", "Lutrinae_sampled", "Enhydra_terminal_branch"),
    native_branch_label = c("B593", "B592", "B568", "B591", "B281", "B573", "B501", "B500", "B493",
      "B492", "B491", "B313", "B312", "B13", "B470", "B483", "B197"),
    branch_id = c("B593", "B592", "B544", "B591", "B548", "B558", "B410", "B409", "B394", "B393",
      "B392", "B26", "B25", "B22", "B349", "B380", "B374"),
    target_type = c(rep("internal_clade", 4), "terminal_lineage", rep("internal_clade", 8),
      "terminal_lineage", rep("internal_clade", 2), "terminal_lineage"),
    target_note = c(
      "Exact internal clade: Hippopotamus_amphibius plus sampled Cetacea; native B593 maps to old B593.",
      "Exact internal clade: sampled crown Cetacea; native B592 maps to old B592.",
      "Exact internal clade: sampled Mysticeti; native B568 maps to old B544.",
      "Exact internal clade: sampled Odontoceti; native B591 maps to old B591.",
      "Single sampled Platanista branch; native B281 maps to old B548.",
      "Exact internal clade: Inia_geoffrensis plus Lipotes_vexillifer; native B573 maps to old B558.",
      "Exact internal clade: sampled Pinnipedia; native B501 maps to old B410.",
      "Exact internal clade: sampled Phocidae; native B500 maps to old B409.",
      "Exact internal clade: Odobenus plus sampled otariids; native B493 maps to old B394.",
      "Exact internal clade: sampled Otariidae; native B492 maps to old B393.",
      "Exact internal clade: Eumetopias plus Zalophus; native B491 maps to old B392.",
      "Exact internal clade: sampled Sirenia; native B313 maps to old B26.",
      "Exact internal clade: Dugong_dugon plus Hydrodamalis_gigas; native B312 maps to old B25.",
      "Single sampled Trichechus branch; native B13 maps to old B22.",
      "Exact internal clade: sampled brown bear plus polar bear; native B470 maps to old B349.",
      "Exact internal clade: sampled Pteronura/Lontra/Enhydra/Aonyx/Lutra; native B483 maps to old B380.",
      "Single sampled sea otter branch; native B197 maps to old B374."), stringsAsFactors = FALSE)
  targets$branch_label_found <- targets$branch_id %in% branch_ids
  if (!all(targets$branch_label_found))
    .lt_marine_abort("A frozen descriptive projection target is absent; no branch substitution.", "projection_target")
  targets
}

.lt_marine_internal_projections <- function(gbi, models) {
  targets <- .lt_marine_internal_targets(colnames(gbi))
  result <- fits <- vector("list", 2L)
  for (j in 1:2) {
    trait <- c("marine_binary", "binary_aquatic_dependence")[j]
    model <- models[[c("fix_marine_binary", "fix_aquatic_v2")[j]]]
    terminal <- model$terminal
    g0 <- .lt_s2_group_stats(gbi, terminal$branch[terminal$response == 0])
    g1 <- .lt_s2_group_stats(gbi, terminal$branch[terminal$response == 1])
    tt <- .lt_s2_welch_stats(g0, g1, "nested_multiply")
    tt <- tt[is.finite(tt$pvalue), , drop = FALSE]
    tt <- tt[order(tt$pvalue, decreasing = FALSE), , drop = FALSE]
    genes <- tt$gene[.lt_marine_historical_fdr(tt$pvalue)]
    fit <- .lt_s2_fit_full(gbi, terminal, genes, 20260701L + 900000L + j)
    if (is.null(fit$fit)) .lt_marine_abort("No usable fitted terminal-only projection model.", "projection_model")
    projected <- .lt_s2_project_fit(gbi, fit, targets$branch_id)
    fit$candidate_genes <- genes; fit$terminal_screen <- tt
    fits[[j]] <- fit
    result[[j]] <- cbind(targets, trait = trait,
      model = if (j == 1L) "marine_terminal_only_full_data_projection" else "binary_aquatic_dependence_terminal_only_full_data_projection",
      context = "terminal_only_internal_projection", predicted_probability = projected$probability,
      logit = projected$logit, n_candidate_features = length(genes),
      n_selected_predictors = sum(fit$beta != 0), lambda_min = fit$lambda_min,
      stringsAsFactors = FALSE)
  }
  list(projections = do.call(rbind, result), targets = targets, models = fits,
    provenance = list(evidence = "descriptive_not_OOF_or_ancestral_habitat_inference",
      feature_screen_domain = "eligible_terminals_only", baseline_71_101_models_used = FALSE))
}

.lt_marine_focal_clades <- function() {
  c(Orcinus_orca = "cetacean", Zalophus_californianus = "pinniped", Leptonychotes_weddellii = "pinniped",
    Odobenus_rosmarus_divergens = "pinniped", Ursus_maritimus = "marine edge", Enhydra_lutris_kenyoni = "marine edge",
    Dugong_dugon = "sirenian", Trichechus_manatus_latirostris = "sirenian", Hydrodamalis_gigas = "sirenian",
    Hippopotamus_amphibius = "non-marine aquatic control", Lutra_lutra = "non-marine aquatic control",
    Lontra_canadensis = "non-marine aquatic control", Pteronura_brasiliensis = "non-marine aquatic control",
    Aonyx_cinereus = "non-marine aquatic control", Platanista_minor = "freshwater cetacean",
    Inia_geoffrensis = "freshwater cetacean", Lipotes_vexillifer = "freshwater cetacean")
}

.lt_marine_fingerprints <- function(gbi, models, full_data, traits, annotations, nc_oof) {
  # A serialized fit does not load its S3 methods in a new R session.
  if (!requireNamespace("glmnet", quietly = TRUE))
    .lt_marine_abort("Install glmnet to use a stored fitted-model dependency.", "dependency")
  clades <- .lt_marine_focal_clades()
  species <- names(clades)
  .lt_assert_unique_ids(annotations$gene, "final Table S5 annotation keys")
  if (!all(c("lasso_group", "display_module") %in% names(annotations)))
    .lt_marine_abort("Frozen final Table S5 display columns are absent.", "annotation_fields")
  fingerprints <- projections <- list()
  for (j in 1:2) {
    run_id <- c("fix_marine_binary", "fix_aquatic_v2")[j]
    model_name <- c("marine", "binary_aquatic_dependence")[j]
    fit <- full_data$models[[run_id]]
    genes <- names(fit$beta)[fit$beta != 0]
    coefficients <- fit$beta[genes]
    intercept <- unname(as.matrix(stats::coef(fit$fit))[, 1]["(Intercept)"])
    ordered <- genes[order(coefficients, decreasing = FALSE)]
    for (s in species) {
      terminal <- models[[run_id]]$terminal
      branch <- terminal$branch[match(s, terminal$species)]
      if (is.na(branch)) .lt_marine_abort("Missing declared focal species.", "projection_target")
      raw <- as.numeric(gbi[genes, branch]); names(raw) <- genes
      imputed <- is.na(raw); filled <- raw
      filled[imputed] <- fit$design$impute_means[genes][imputed]
      scaled <- (filled - fit$design$scale_means[genes]) / fit$design$scale_sds[genes]
      contribution <- scaled * coefficients
      logit <- intercept + sum(contribution)
      probability <- 1 / (1 + exp(-logit))
      ti <- match(s, traits$species)
      category <- if (traits$marine_binary[ti] == 1) "marine" else
        if (traits$aquatic_v2_status[ti] == 1) "non-marine aquatic" else
          if (traits$aquatic_v2_status[ti] == 0.5) "semi-aquatic" else "terrestrial"
      parts <- strsplit(s, "_", fixed = TRUE)[[1]]
      label <- paste(parts[1:2], collapse = " ")
      ann <- annotations[match(genes, annotations$gene), , drop = FALSE]
      predictor_class <- ifelse(is.na(ann$lasso_group) | ann$lasso_group == "", "not_classified", ann$lasso_group)
      module <- ifelse(!is.na(ann$display_module) & ann$display_module != "", ann$display_module, "Table S5-only / unassigned")
      row <- data.frame(species = s, display_label = label, clade_group = unname(clades[s]),
        trait_category = category, aquaticity_score = traits$aquaticity_score_sum_0_18[ti],
        model = model_name, gene = genes, gene_order = match(genes, ordered),
        predictor_class = predictor_class, module = module, coefficient = as.numeric(coefficients),
        scaled_GBI = as.numeric(scaled), contribution = as.numeric(contribution),
        abs_contribution = abs(as.numeric(contribution)),
        contribution_sign = ifelse(contribution > 0, "positive", ifelse(contribution < 0, "negative", "zero")),
        GBI_raw_if_available = raw, imputed_flag = imputed, stringsAsFactors = FALSE)
      row <- row[order(row$gene_order), , drop = FALSE]
      fingerprints[[length(fingerprints) + 1L]] <- row
      o <- nc_oof[nc_oof$model == model_name & nc_oof$species == s, , drop = FALSE]
      if (nrow(o) != 1L) .lt_marine_abort("Missing/ambiguous computed nested NC export row.", "dependency")
      projections[[length(projections) + 1L]] <- data.frame(species = s, display_label = label,
        clade_group = unname(clades[s]), trait_category = category,
        aquaticity_score = traits$aquaticity_score_sum_0_18[ti], model = model_name,
        n_predictors = nrow(row), n_imputed_predictors = sum(row$imputed_flag), intercept = intercept,
        contribution_sum = sum(row$contribution), fitted_logit = logit, fitted_probability = probability,
        OOF_probability_if_available = o$probability, OOF_logit_if_available = o$logit,
        fitted_minus_OOF = probability - o$probability,
        evidence_label = "descriptive full-data fitted projection after corrected preprocessing; not held-out validation",
        stringsAsFactors = FALSE)
    }
  }
  fp <- do.call(rbind, fingerprints)
  proj <- do.call(rbind, projections)
  proj <- proj[order(match(proj$species, species), proj$model), , drop = FALSE]
  list(fingerprints = fp, projections = proj,
    provenance = list(model_source = "same_run_full_data_baseline_fits",
      terminal_only_model_used = FALSE, evidence = "descriptive_not_validation"))
}
