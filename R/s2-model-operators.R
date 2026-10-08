# Pure data operators shared by explicit S2 recipes. Marine identities and
# historical file paths do not belong here. Arithmetic variants are deliberate:
# the two historical screen sources used different floating-point expressions.
.lt_s2_group_stats <- function(values, branches) {
  if (!length(branches)) {
    z <- stats::setNames(rep(0, nrow(values)), rownames(values))
    return(list(n = z, sum = z, sumsq = z))
  }
  x <- values[, branches, drop = FALSE]
  list(n = rowSums(!is.na(x)), sum = rowSums(x, na.rm = TRUE),
       sumsq = rowSums(x * x, na.rm = TRUE))
}

.lt_s2_subtract_stats <- function(base, heldout) {
  list(n = base$n - heldout$n, sum = base$sum - heldout$sum,
       sumsq = base$sumsq - heldout$sumsq)
}

.lt_s2_welch_stats <- function(g0, g1, arithmetic) {
  arithmetic <- match.arg(arithmetic, c("nested_multiply", "phase13_power"))
  n0 <- g0$n; n1 <- g1$n
  ok <- n0 > 1 & n1 > 1
  m0 <- m1 <- v0 <- v1 <- tv <- pv <- rep(NA_real_, length(n0))
  m0[ok] <- g0$sum[ok] / n0[ok]
  m1[ok] <- g1$sum[ok] / n1[ok]
  if (arithmetic == "nested_multiply") {
    s0 <- g0$sum[ok] * g0$sum[ok]
    s1 <- g1$sum[ok] * g1$sum[ok]
  } else {
    s0 <- g0$sum[ok]^2
    s1 <- g1$sum[ok]^2
  }
  v0[ok] <- pmax((g0$sumsq[ok] - s0 / n0[ok]) / (n0[ok] - 1), 0)
  v1[ok] <- pmax((g1$sumsq[ok] - s1 / n1[ok]) / (n1[ok] - 1), 0)
  se2 <- v0 / n0 + v1 / n1
  keep <- ok & is.finite(se2) & se2 > 0
  tv[keep] <- (m0[keep] - m1[keep]) / sqrt(se2[keep])
  numerator <- if (arithmetic == "nested_multiply") se2[keep] * se2[keep] else se2[keep]^2
  df <- numerator / (((v0[keep] / n0[keep])^2 / (n0[keep] - 1)) +
                       ((v1[keep] / n1[keep])^2 / (n1[keep] - 1)))
  pv[keep] <- 2 * stats::pt(-abs(tv[keep]), df = df)
  data.frame(gene = names(n0), tvalue = tv, pvalue = pv,
    mean1 = m0, mean2 = m1, marine = n1, non_marine = n0,
    stringsAsFactors = FALSE, check.names = FALSE)
}

.lt_s2_screen_cache <- function(values, states) {
  list(state = states, base0 = .lt_s2_group_stats(values, names(states)[states == 0]),
       base1 = .lt_s2_group_stats(values, names(states)[states == 1]))
}

.lt_s2_fold_features <- function(values, cache, held_branches, arithmetic) {
  h0 <- held_branches[cache$state[held_branches] == 0]
  h1 <- held_branches[cache$state[held_branches] == 1]
  tt <- .lt_s2_welch_stats(
    .lt_s2_subtract_stats(cache$base0, .lt_s2_group_stats(values, h0)),
    .lt_s2_subtract_stats(cache$base1, .lt_s2_group_stats(values, h1)), arithmetic)
  tt <- tt[is.finite(tt$pvalue), , drop = FALSE]
  tt <- tt[order(tt$pvalue), , drop = FALSE]
  selected <- .lt_marine_historical_fdr(tt$pvalue)
  list(genes = tt$gene[selected], tested = tt, n_removed = length(h0) + length(h1),
       n_screened = sum(cache$state %in% c(0, 1)) - length(h0) - length(h1))
}

.lt_s2_prepare_design <- function(values, train, test, genes) {
  x <- t(values[genes, train, drop = FALSE])
  xt <- t(values[genes, test, drop = FALSE])
  impute <- colMeans(x, na.rm = TRUE)
  keep <- is.finite(impute)
  drop_all <- sum(!keep)
  x <- x[, keep, drop = FALSE]; xt <- xt[, keep, drop = FALSE]
  impute <- impute[keep]
  if (!length(impute)) return(list(ok = FALSE,
    status = "no_features_after_all_missing_drop", n_dropped_all = drop_all, n_dropped_zero = 0L))
  for (k in seq_along(impute)) {
    if (anyNA(x[, k])) x[is.na(x[, k]), k] <- impute[[k]]
    if (anyNA(xt[, k])) xt[is.na(xt[, k]), k] <- impute[[k]]
  }
  center <- colMeans(x)
  scale <- apply(x, 2, stats::sd)
  keep <- is.finite(scale) & scale > 0
  drop_zero <- sum(!keep)
  if (!any(keep)) return(list(ok = FALSE,
    status = "no_features_after_zero_variance_drop", n_dropped_all = drop_all,
    n_dropped_zero = drop_zero))
  x <- x[, keep, drop = FALSE]; xt <- xt[, keep, drop = FALSE]
  center <- center[keep]; scale <- scale[keep]; impute <- impute[keep]
  list(ok = TRUE,
    x_train = sweep(sweep(x, 2, center, "-"), 2, scale, "/"),
    x_test = sweep(sweep(xt, 2, center, "-"), 2, scale, "/"),
    feature_names = colnames(x), impute_means = impute,
    scale_means = center, scale_sds = scale,
    n_dropped_all = drop_all, n_dropped_zero = drop_zero)
}

.lt_s2_with_preserved_rng <- function(code) {
  kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(kind))
    if (had_seed) assign(".Random.seed", seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  force(code)
}

.lt_s2_fit_design <- function(design, y, seed) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("This explicitly selected S2 model needs the optional 'glmnet' package.", call. = FALSE)
  }
  .lt_s2_with_preserved_rng({
    warnings <- character()
    set.seed(seed)
    cv <- withCallingHandlers(glmnet::cv.glmnet(
      x = as.matrix(design$x_train), y = y, family = "binomial",
      foldid = seq_along(y), standardize = FALSE, type.measure = "deviance"),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      })
    fit_warnings <- character()
    fit <- withCallingHandlers(glmnet::glmnet(x = as.matrix(design$x_train), y = y,
      family = "binomial", lambda = cv$lambda.min, standardize = FALSE),
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        # Retain the original emitted warning as well as recording it.
      })
    beta <- as.matrix(fit$beta)[, 1]
    list(fit = fit, beta = beta, lambda_min = cv$lambda.min,
         warnings = unique(warnings), fit_warnings = unique(fit_warnings),
         solver_status = fit$jerr)
  })
}

.lt_s2_auc <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]; p <- p[keep]
  np <- sum(y == 1); nn <- sum(y == 0)
  if (!np || !nn) return(NA_real_)
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1]) - np * (np + 1) / 2) / (np * nn)
}

.lt_s2_model_inputs <- function(gbi, branch_states, terminals, folds) {
  if (!is.matrix(gbi) || !is.numeric(gbi) || any(is.infinite(gbi)) ||
      any(gbi < 0, na.rm = TRUE)) stop("Supply a numeric nonnegative GBI matrix with NA availability.", call. = FALSE)
  .lt_assert_unique_ids(rownames(gbi), "GBI gene keys")
  .lt_assert_unique_ids(colnames(gbi), "GBI branch keys")
  if (!nrow(gbi) || !ncol(gbi)) stop("GBI axes cannot be empty.", call. = FALSE)
  if (!is.numeric(branch_states) || !is.null(dim(branch_states)) ||
      anyNA(branch_states) || !all(branch_states %in% c(0, 0.5, 1)))
    stop("Declare keyed branch states as 0, 0.5 (excluded), or 1; continuous traits are unsupported.", call. = FALSE)
  .lt_assert_unique_ids(names(branch_states), "branch state keys")
  if (!setequal(names(branch_states), colnames(gbi)))
    stop("Branch states must exactly cover the GBI scientific keys.", call. = FALSE)
  fields <- c("species", "branch", "genus", "response")
  if (!is.data.frame(terminals) || !all(fields %in% names(terminals)))
    stop("Terminals need species, branch, genus and response fields.", call. = FALSE)
  .lt_assert_unique_ids(terminals$species, "terminal species keys")
  .lt_assert_unique_ids(terminals$branch, "terminal branch keys")
  if (!nrow(terminals) || !all(terminals$branch %in% colnames(gbi)) ||
      !is.character(terminals$genus) || anyNA(terminals$genus) || any(!nzchar(terminals$genus)))
    stop("Every terminal needs an existing branch and an explicit nonempty genus key.", call. = FALSE)
  if (!is.numeric(terminals$response) || anyNA(terminals$response) ||
      !all(terminals$response %in% c(0, 0.5, 1)))
    stop("Terminal responses must be explicit 0/0.5/1; no inferred discretization.", call. = FALSE)
  if (!identical(as.numeric(branch_states[terminals$branch]), as.numeric(terminals$response)))
    stop("Terminal response and corresponding branch state disagree.", call. = FALSE)
  if (!is.data.frame(folds) || !all(c("fold_id", "genus", "seed") %in% names(folds)))
    stop("Supply the full ordered fold_id/genus/seed ledger.", call. = FALSE)
  .lt_assert_unique_ids(folds$genus, "fold genus keys")
  .lt_assert_unique_ids(as.character(folds$fold_id), "fold identifiers")
  if (!setequal(folds$genus, terminals$genus) || !is.numeric(folds$seed) ||
      anyNA(folds$seed) || any(!is.finite(folds$seed)) ||
      any(folds$seed != trunc(folds$seed)) || any(abs(folds$seed) > .Machine$integer.max))
    stop("Fold ledger must cover all genera, including excluded-only genera, with valid integer seeds.", call. = FALSE)
  if (!any(terminals$response %in% c(0, 1))) stop("No eligible terminal responses.", call. = FALSE)
  branch_states[colnames(gbi)]
}

.lt_s2_fit_fold <- function(gbi, states, terminals, fold, recipe, global_features, cache) {
  eligible <- terminals$response != 0.5
  test <- which(eligible & terminals$genus == fold$genus)
  train <- which(eligible & terminals$genus != fold$genus)
  base <- list(fold_id = fold$fold_id, genus = fold$genus, seed = fold$seed,
    train_species = terminals$species[train], test_species = terminals$species[test],
    train_branches = terminals$branch[train], test_branches = terminals$branch[test],
    y_train = terminals$response[train], y_test = terminals$response[test])
  result <- list(identity = base, genes = character(), n_tested = NA_integer_,
    n_removed = 0L, n_screened = NA_integer_, prediction = numeric(),
    n_dropped_all = NA_integer_, n_dropped_zero = NA_integer_,
    n_used = 0L, beta = numeric(), lambda_min = NA_real_, warnings = character(),
    fit_warnings = character(), solver_status = NA_integer_)
  if (!length(test)) {
    result$status <- "no_evaluable_test_species_after_exclusion"
    return(result)
  }
  if (!length(train)) stop("A held-out genus leaves no eligible training terminals.", call. = FALSE)
  if (recipe == "global_screen_foldwise") result$genes <- global_features else {
    held <- terminals$branch[terminals$genus == fold$genus]
    fr <- .lt_s2_fold_features(gbi, cache, held,
      if (recipe == "nested_phase12B") "nested_multiply" else "phase13_power")
    result$genes <- fr$genes; result$n_tested <- nrow(fr$tested)
    result$n_removed <- fr$n_removed; result$n_screened <- fr$n_screened
  }
  null <- function(status, all = NA_integer_, zero = NA_integer_) {
    result$status <- status
    result$prediction <- rep(mean(base$y_train), length(test))
    result$n_dropped_all <- all; result$n_dropped_zero <- zero
    result
  }
  if (length(unique(base$y_train)) < 2)
    return(null("null_intercept_only_training_one_class"))
  if (!length(result$genes) && recipe != "global_screen_foldwise")
    return(null("null_intercept_only_no_FDR_features_in_training_ttest", 0L, 0L))
  design <- .lt_s2_prepare_design(gbi, base$train_branches, base$test_branches, result$genes)
  if (!design$ok) return(null(paste0("null_intercept_only_", design$status),
    design$n_dropped_all, design$n_dropped_zero))
  model <- .lt_s2_fit_design(design, base$y_train, fold$seed)
  result$prediction <- as.numeric(stats::predict(model$fit,
    newx = as.matrix(design$x_test), type = "response")[, 1])
  result$beta <- model$beta[model$beta != 0]
  result$lambda_min <- model$lambda_min; result$warnings <- model$warnings
  result$fit_warnings <- model$fit_warnings; result$solver_status <- model$solver_status
  result$n_dropped_all <- design$n_dropped_all
  result$n_dropped_zero <- design$n_dropped_zero
  result$n_used <- ncol(design$x_train)
  result$status <- if (length(result$beta)) "fit_nonzero_predictors" else "fit_no_selected_predictors"
  result
}

#' Execute an explicitly selected S2 grouped validation recipe
#'
#' Reuses the frozen S2 operator contracts with other keyed data and an explicitly
#' supplied binary/excluded trait and full fold ledger. This is not a Marine
#' snapshot validator, automatic trait inference, null generator or new default.
#' It does not compute a baseline or ASR. Internal branches participate in the
#' selected nested screen but never in terminal model preprocessing or fitting.
#'
#' @param gbi Nonnegative numeric gene-by-branch GBI matrix; NA is unavailable.
#' @param branch_states Named numeric branch states 0, 0.5 (excluded), or 1.
#' @param terminals Data frame with unique character species and branch keys,
#'   character genus, and numeric response. Row order defines sample order.
#' @param folds Complete ordered data frame with unique fold_id, genus and
#'   integer-valued seed, including excluded-only genera. No folds are generated.
#' @param recipe Required explicit choice: `"nested_phase12B"`,
#'   `"nested_phase13"`, or `"global_screen_foldwise"`. The first two retain
#'   their historical arithmetic and strict-prefix 0.01 screen; the last needs
#'   explicitly supplied global features and does not redo supervised selection.
#' @param global_features Ordered gene keys required for the global-screen
#'   recipe; optional reference keys for nested recipe overlap reporting.
#' @param progress Optional callback receiving completed fold count, total and
#'   current result. It cannot replace the computed result.
#' @param cores Positive integer worker count. More than one uses ordered fork
#'   batches on supported platforms; unsupported parallelism fails explicitly.
#' @return List of per-fold results and a separate input/specification provenance
#'   snapshot. Completion is not scientific or frozen-replay certification.
#' @export
lt_s2_gloocv <- function(gbi, branch_states, terminals, folds, recipe,
                        global_features = NULL, progress = NULL, cores = 1L) {
  if (missing(recipe)) stop("Choose an explicit S2 validation recipe; there is no automatic default.", call. = FALSE)
  recipe <- match.arg(recipe, c("nested_phase12B", "nested_phase13", "global_screen_foldwise"))
  states <- .lt_s2_model_inputs(gbi, branch_states, terminals, folds)
  .lt_s2_gloocv_execute(gbi, states, terminals, folds, recipe, global_features,
    progress, cores, paste(deparse(match.call()), collapse = " "))
}

# Private numerical execution, entered only after the ordinary keyed contract
# or the independently source-bound historical Marine adapter has validated.
# This is not an exported bypass or a recipe-selectable validation exception.
.lt_s2_gloocv_execute <- function(gbi, branch_states, terminals, folds, recipe,
                                global_features, progress, cores, command) {
  states <- branch_states[colnames(gbi)]
  if (recipe == "global_screen_foldwise" && is.null(global_features))
    stop("Supply the ordered global-screen feature set explicitly.", call. = FALSE)
  if (!is.null(global_features)) {
    if (!is.character(global_features) || anyNA(global_features) ||
        anyDuplicated(global_features) || !all(global_features %in% rownames(gbi)))
      stop("Global features must be unique existing gene keys in declared order.", call. = FALSE)
  }
  if (!is.null(progress) && !is.function(progress)) stop("progress must be a function or NULL.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) ||
      !is.finite(cores) || cores < 1 || cores != trunc(cores))
    stop("cores must be a positive integer.", call. = FALSE)
  if (cores > 1L && .Platform$OS.type == "windows")
    stop("Forked execution is unavailable here; explicitly choose cores = 1.", call. = FALSE)
  cache <- if (recipe == "global_screen_foldwise") NULL else .lt_s2_screen_cache(gbi, states)
  results <- vector("list", nrow(folds))
  for (chunk in split(seq_len(nrow(folds)), ceiling(seq_len(nrow(folds)) / cores))) {
    compute <- function(i) .lt_s2_fit_fold(gbi, states, terminals,
      folds[i, , drop = FALSE], recipe, global_features, cache)
    batch <- if (cores == 1L) lapply(chunk, compute) else
      parallel::mclapply(chunk, compute, mc.cores = min(cores, length(chunk)), mc.set.seed = FALSE)
    if (any(vapply(batch, inherits, logical(1), "try-error")))
      stop("A parallel S2 fold failed; no skipped-fold or null-model fallback was substituted.", call. = FALSE)
    for (k in seq_along(chunk)) {
      i <- chunk[[k]]; results[[i]] <- batch[[k]]
      if (!is.null(progress)) progress(i, nrow(folds), results[[i]])
    }
  }
  fingerprint <- function(x) digest::digest(serialize(x, NULL, version = 3),
    algo = "sha256", serialize = FALSE)
  list(folds = results, provenance = list(
    recipe = recipe, input_role = "user_declared_not_frozen_Marine_certification",
    input_sha256 = c(gbi = fingerprint(gbi), branch_states = fingerprint(branch_states),
      terminals = fingerprint(terminals), folds = fingerprint(folds),
      global_features = fingerprint(global_features)),
    ordered_gene_keys = rownames(gbi), ordered_branch_keys = colnames(gbi),
    terminal_ledger = terminals, fold_ledger = folds, RNGkind = RNGkind(),
    requested_cores = cores,
    glmnet_version = if (requireNamespace("glmnet", quietly = TRUE))
      as.character(utils::packageVersion("glmnet")) else NA_character_,
    command = command,
    inner_fold_rule = "seq_along(eligible training terminal responses)",
    caller_RNG_preserved = TRUE, silent_fallback = FALSE))
}
