.lt_s2_fit_full <- function(gbi, terminals, genes, seed) {
  eligible <- terminals$response != 0.5
  train <- terminals$branch[eligible]
  y <- as.numeric(terminals$response[eligible])
  # The full-data source checks the single-class case before variance filtering.
  raw <- t(gbi[genes, train, drop = FALSE])
  means <- colMeans(raw, na.rm = TRUE)
  null <- function(status, all, zero = NA_integer_) list(status = status,
    beta = stats::setNames(rep(0, length(genes)), genes), lambda_min = NA_real_,
    fit = NULL, design = NULL, n_dropped_all = all, n_dropped_zero = zero,
    y = y, train = train, seed = seed, warnings = character())
  if (!any(is.finite(means)) || length(unique(y)) < 2)
    return(null("Null fit: no usable features or one training class.", sum(!is.finite(means))))
  design <- .lt_s2_prepare_design(gbi, train, character(), genes)
  if (!design$ok) return(null("Null fit: no usable nonzero-variance features.",
    design$n_dropped_all, design$n_dropped_zero))
  model <- .lt_s2_fit_design(design, y, seed)
  list(status = "fit", beta = model$beta, lambda_min = model$lambda_min,
    fit = model$fit, design = design,
    n_dropped_all = design$n_dropped_all, n_dropped_zero = design$n_dropped_zero,
    y = y, train = train, seed = seed, warnings = model$warnings,
    fit_warnings = model$fit_warnings, solver_status = model$solver_status)
}

.lt_s2_project_fit <- function(gbi, model, target_branches) {
  if (is.null(model$fit) || is.null(model$design))
    stop("Projection requires a fitted source model; no substitute fit is used.", call. = FALSE)
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Install glmnet to use a stored fitted-model dependency.", call. = FALSE)
  design <- model$design
  x <- t(gbi[design$feature_names, target_branches, drop = FALSE])
  for (k in seq_along(design$impute_means))
    if (anyNA(x[, k])) x[is.na(x[, k]), k] <- design$impute_means[[k]]
  x <- sweep(sweep(x, 2, design$scale_means, "-"), 2, design$scale_sds, "/")
  p <- as.numeric(stats::predict(model$fit, newx = as.matrix(x), type = "response")[, 1])
  list(probability = p, logit = stats::qlogis(pmin(pmax(p, 1e-12), 1 - 1e-12)),
       scaled_values = x)
}

.lt_marine_full_data <- function(gbi, models, progress = NULL) {
  run_ids <- c("fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear",
    "fix_drop_sea_otter", "fix_drop_pinniped", "fix_whale_only", "fix_pinniped_only",
    "fix_aquatic_v2", "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia",
    "fix_aquatic_v2_noMarineEdge", "fix_aquatic_v2_noNonMarineCarnivores",
    "fix_aquatic_v2_noNonMarineRodents", "fix_aquatic_v2_noHippoLechwe")
  if (!identical(names(models), run_ids))
    .lt_marine_abort("Full-data models must retain their frozen 14-run order.", "model_order")
  axes <- c(rep("marine", 7L), rep("binary aquatic-dependence", 7L))
  drop <- c("baseline", "Cetacea", "Ursus_maritimus", "Enhydra_lutris_kenyoni", "Pinnipedia",
    "Cetacea only", "Pinnipedia only", "baseline", "Cetacea", "Pinnipedia", "MarineEdge",
    "NonMarineCarnivores", "NonMarineRodents", "HippoLechwe")
  fits <- coefficients <- summaries <- vector("list", 14L)
  names(fits) <- run_ids
  for (j in seq_along(models)) {
    model <- models[[j]]; genes <- model$global_screen$significant$gene
    seed <- 20260523L + 900000L + length(genes) + j
    fit <- .lt_s2_fit_full(gbi, model$terminal, genes, seed)
    fits[[j]] <- fit; beta <- fit$beta
    handling <- if (j > 7L)
      "semi-aquatic aquatic_v2=0.5 labels excluded from training, testing, imputation, scaling, and AUC" else
      if (any(model$terminal$response == 0.5))
        "0.5 labels for dropped/held-out marine taxa excluded from training, testing, imputation, scaling, and AUC" else
        "all terminal species included; all non-positive terminal species coded 0"
    coefficients[[j]] <- data.frame(run_id = run_ids[j], trait_axis = axes[j],
      trait_column = model$trait_column, drop_rule = drop[j], gene = names(beta),
      coefficient = as.numeric(beta), selected = beta != 0,
      coefficient_sign = ifelse(beta > 0, "positive", ifelse(beta < 0, "negative", "zero")),
      lambda_used = fit$lambda_min, stringsAsFactors = FALSE)
    max_i <- if (!length(beta) || max(abs(beta)) == 0) NA_integer_ else which.max(abs(beta))
    notes <- if (fit$status != "fit") fit$status else paste0("Dropped all-missing features: ",
      fit$n_dropped_all, "; dropped zero-variance features: ", fit$n_dropped_zero,
      ". Full-data fit is not validation.")
    summaries[[j]] <- data.frame(run_id = run_ids[j], trait_axis = axes[j], drop_rule = drop[j],
      n_global_FDR_candidate_features = length(genes), n_terminal_samples_used = length(fit$y),
      n_positive = sum(fit$y == 1), n_negative = sum(fit$y == 0), semi_aquatic_handling = handling,
      n_selected_predictors = sum(beta != 0), n_positive_coefficients = sum(beta > 0),
      n_negative_coefficients = sum(beta < 0),
      max_abs_coef_gene = ifelse(is.na(max_i), "", names(beta)[max_i]),
      max_abs_coef = ifelse(is.na(max_i), 0, abs(beta[max_i])), lambda_used = fit$lambda_min,
      model_family = "binomial glmnet",
      preprocessing_policy = "eligible terminal samples only; terminal mean imputation and scaling computed over eligible terminal samples; glmnet standardize=FALSE",
      notes = notes, stringsAsFactors = FALSE)
    if (!is.null(progress)) progress(j, 14L, summaries[[j]])
  }
  marine <- names(fits[[1]]$beta)[fits[[1]]$beta != 0]
  aquatic <- names(fits[[8]]$beta)[fits[[8]]$beta != 0]
  sets <- list(setdiff(marine, aquatic), setdiff(aquatic, marine), intersect(marine, aquatic))
  partition <- data.frame(partition = c("marine_only", "aquatic_only", "shared"),
    n_genes = lengths(sets), genes = vapply(sets, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE)
  list(models = fits, coefficients = do.call(rbind, coefficients),
    summary = do.call(rbind, summaries), partition = partition,
    provenance = list(evidence_class = "full_data_fit_not_validation",
      seed_rule = "20260523 + 900000 + computed_global_feature_count + ordered_run_index",
      run_ids = run_ids, archived_coefficients_used_as_input = FALSE))
}
