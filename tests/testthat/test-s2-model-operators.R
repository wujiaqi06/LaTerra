s2_other_data <- function() {
  x <- outer(seq_len(4), seq_len(9), function(i, j) i + (j %% 4) / 3)
  dimnames(x) <- list(paste0("custom-g", 4:1), paste0("coord-", 9:1))
  x[1, 2] <- NA_real_
  states <- stats::setNames(c(0, 1, 0, 1, 0.5, 0.5, 0, 1, 0), colnames(x))
  terminals <- data.frame(species = paste0("taxon ", 1:6),
    branch = colnames(x)[1:6], genus = rep(c("group-Z", "group-A", "excluded"), each = 2),
    response = unname(states[1:6]))
  folds <- data.frame(fold_id = c("f9", "f2", "fX"),
    genus = c("group-A", "excluded", "group-Z"), seed = c(17L, 23L, 19L))
  list(gbi = x, branch_states = states, terminals = terminals, folds = folds)
}

test_that("other data and explicit traits run without Marine hashes or labels", {
  d <- s2_other_data()
  run <- do.call(lt_s2_gloocv, c(d, list(recipe = "global_screen_foldwise", global_features = character())))
  expect_length(run$folds, 3)
  expect_identical(vapply(run$folds, function(z) z$identity$fold_id, character(1)), d$folds$fold_id)
  expect_identical(run$folds[[2]]$status, "no_evaluable_test_species_after_exclusion")
  expect_identical(run$folds[[2]]$identity$seed, 23L)
  expect_identical(run$folds[[1]]$status, "null_intercept_only_no_features_after_all_missing_drop")
  expect_identical(run$folds[[1]]$prediction, c(0.5, 0.5))
  expect_identical(run$provenance$ordered_gene_keys, rownames(d$gbi))
  edited <- d; edited$gbi[2, 2] <- 99
  edited$terminals$display_label <- paste("edited", seq_len(6))
  expect_no_error(do.call(lt_s2_gloocv, c(edited,
    list(recipe = "global_screen_foldwise", global_features = character()))))
  reversed <- d; reversed$branch_states <- rev(d$branch_states)
  replay <- do.call(lt_s2_gloocv, c(reversed,
    list(recipe = "global_screen_foldwise", global_features = character())))
  expect_identical(run$folds, replay$folds)
})

test_that("explicit scientific keys and full fold coverage are mandatory", {
  d <- s2_other_data()
  call <- function(x) do.call(lt_s2_gloocv, c(x, list(recipe = "nested_phase12B")))
  expect_error(do.call(lt_s2_gloocv, d), "explicit")
  bad <- d; bad$terminals$response[1] <- 0.2
  expect_error(call(bad), "0/0.5/1")
  bad <- d; bad$terminals$response[1] <- 1
  expect_error(call(bad), "disagree")
  bad <- d; bad$folds <- bad$folds[-2, ]
  expect_error(call(bad), "all genera")
  bad <- d; bad$branch_states <- bad$branch_states[-1]
  expect_error(call(bad), "exactly cover")
  bad <- d; bad$gbi[1, 1] <- Inf
  expect_error(call(bad), "numeric nonnegative")
  bad <- d; rownames(bad$gbi)[2] <- rownames(bad$gbi)[1]
  expect_error(call(bad), "unique|duplicate")
  expect_error(do.call(lt_s2_gloocv, c(d, list(recipe = "global_screen_foldwise"))), "explicitly")
})

test_that("nested screening subtracts held-out terminals and retains internal branches", {
  d <- s2_other_data(); cache <- LaTerra:::.lt_s2_screen_cache(d$gbi, d$branch_states)
  fr <- LaTerra:::.lt_s2_fold_features(d$gbi, cache, d$terminals$branch[1:2], "nested_multiply")
  expect_equal(fr$n_removed, 2)
  expect_equal(fr$n_screened, 5)
  expect_identical(fr$tested$pvalue, sort(fr$tested$pvalue))
  result <- do.call(lt_s2_gloocv, c(d, list(recipe = "nested_phase12B")))
  expect_identical(result$folds[[1]]$identity$train_branches, d$terminals$branch[1:2])
  expect_equal(result$folds[[1]]$n_removed, 2)
  expect_equal(result$folds[[1]]$n_screened, 5)
  expect_identical(result$folds[[2]]$status, "no_evaluable_test_species_after_exclusion")
})

test_that("preprocessing is training-only and preserves declared feature order", {
  x <- rbind(g3 = c(1, NA, 3, 1e6), g1 = c(NA, NA, NA, 12),
    g2 = c(4, 4, 4, 500), g4 = c(2, 3, 4, NA))
  colnames(x) <- c("a", "b", "c", "held")
  prep <- LaTerra:::.lt_s2_prepare_design
  d <- prep(x, c("a", "b", "c"), "held", c("g4", "g2", "g1", "g3"))
  expect_true(d$ok)
  expect_identical(d$feature_names, c("g4", "g3"))
  expect_identical(unname(d$impute_means), c(3, 2))
  expect_identical(unname(d$scale_sds), c(1, 1))
  expect_equal(d$n_dropped_all, 1); expect_equal(d$n_dropped_zero, 1)
  x[, "held"] <- c(0, 15, 88, 50)
  e <- prep(x, c("a", "b", "c"), "held", c("g4", "g2", "g1", "g3"))
  expect_identical(d$x_train, e$x_train)
  expect_identical(d$impute_means, e$impute_means)
  expect_identical(d$scale_sds, e$scale_sds)
  expect_false(identical(d$x_test, e$x_test))
})

test_that("S2 model uses explicit seed and leaves caller RNG intact on success and error", {
  skip_if_not_installed("glmnet")
  x <- rbind(a = (1:24) / 10, b = ((1:24) %% 7) / 5 + 0.1)
  colnames(x) <- paste0("b", 1:24)
  y <- rep(c(0, 1), each = 12)
  design <- LaTerra:::.lt_s2_prepare_design(x, colnames(x), colnames(x)[1:2], rownames(x))
  set.seed(771); seed <- .Random.seed
  fit <- LaTerra:::.lt_s2_fit_design(design, y, 4201L)
  expect_identical(.Random.seed, seed)
  expect_length(fit$beta, 2)
  expect_true(is.finite(fit$lambda_min))
  expect_identical(fit$beta, LaTerra:::.lt_s2_fit_design(design, y, 4201L)$beta)
  expect_error(LaTerra:::.lt_s2_with_preserved_rng({set.seed(12); stop("test error")}), "test error")
  expect_identical(.Random.seed, seed)
  LaTerra:::.lt_s2_with_preserved_rng({
    rm(".Random.seed", envir = .GlobalEnv)
    LaTerra:::.lt_s2_with_preserved_rng(set.seed(3))
    expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  })
  expect_identical(.Random.seed, seed)
})

test_that("public other-data workflow actually fits nontrivial models and predicts", {
  skip_if_not_installed("glmnet")
  y <- rep(0:1, 12)
  x <- rbind(custom_rate_signal = y + 0.2 + (1:24)/100,
    custom_rate_other = abs(sin(1:24)) + 0.1)
  colnames(x) <- paste0("own-coordinate-", 1:24)
  x[2, 3] <- NA_real_
  states <- stats::setNames(as.numeric(y), colnames(x))
  terminal <- data.frame(species = paste0("own species ", 1:24), branch = colnames(x),
    genus = rep(c("set-Z", "set-K", "set-A", "set-W"), each = 6), response = as.numeric(y))
  folds <- data.frame(fold_id = c("user-4", "user-1", "user-7", "user-2"),
    genus = c("set-W", "set-Z", "set-A", "set-K"), seed = 401:404)
  run <- lt_s2_gloocv(x, states, terminal, folds,
    recipe = "global_screen_foldwise", global_features = rownames(x))
  expect_true(all(vapply(run$folds, function(f) length(f$beta) > 0, logical(1))))
  p <- unlist(lapply(run$folds, `[[`, "prediction"))
  expect_length(p, 24L)
  expect_true(all(is.finite(p) & p >= 0 & p <= 1))
  expect_gt(length(unique(p)), 2L)
  expect_identical(run$provenance$terminal_ledger, terminal)
  expect_identical(run$provenance$fold_ledger, folds)
  expect_identical(vapply(run$folds, function(f) f$identity$fold_id, character(1)), folds$fold_id)
  expect_identical(run$provenance$input_role, "user_declared_not_frozen_Marine_certification")
})
