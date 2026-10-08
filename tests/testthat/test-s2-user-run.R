s2_user_example <- function() {
  species <- paste0(rep(c("CustomA", "CustomB", "CustomC"), each = 8), "_", rep(1:8, 3))
  tree <- ape::read.tree(text = paste0("(", paste(vapply(split(species, rep(1:3, each = 8)),
    function(s) paste0("(", paste(s, collapse = ","), ")"), ""), collapse = ","), ");"))
  children <- split(tree$edge[, 2], tree$edge[, 1])
  descend <- function(node) if (node <= length(species)) tree$tip.label[node] else
    unlist(lapply(children[[as.character(node)]], descend), use.names = FALSE)
  sides <- lapply(tree$edge[, 2], descend)
  terminal <- tree$edge[, 2] <= length(species)
  coordinates <- data.frame(branch_id = paste0("suppliedEdge", seq_along(sides)),
    canonical_split_key = paste0("suppliedSplit", seq_along(sides)),
    side_A_taxa = vapply(sides, function(s) paste(sort(s), collapse = ";"), ""),
    side_B_taxa = vapply(sides, function(s) paste(sort(setdiff(tree$tip.label, s)), collapse = ";"), ""),
    branch_type = ifelse(terminal, "terminal", "internal"),
    terminal_taxon = ifelse(terminal, tree$tip.label[tree$edge[, 2]], NA_character_))
  y <- rep(rep(0:1, each = 4), 3)
  positive <- stats::setNames(y, species)[coordinates$terminal_taxon] == 1
  positive[is.na(positive)] <- FALSE
  values <- outer(1:32, seq_along(sides), function(g, b) 2 + g/100 + ((g*7+b*3) %% 11)/30)
  values[1:5, positive] <- values[1:5, positive] * 8
  dimnames(values) <- list(paste0("userGene", 1:32), coordinates$branch_id)
  values[32, 2] <- NA_real_
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  state[32, 2] <- "NA_struct"
  x <- lt_matrix(values, state, fixture_payload())
  trait <- lt_trait("different_binary_trait", species, as.numeric(y), type = "binary")
  groups <- data.frame(species = species, genus = rep(c("declaredGroup1", "declaredGroup2", "declaredGroup3"), each = 8))
  folds <- data.frame(fold_id = c("first", "second", "third"), genus = unique(groups$genus), seed = c(91L, 92L, 93L))
  list(x = x, tree = tree, coordinates = coordinates, trait = trait,
    terminal_groups = groups, folds = folds)
}

test_that("different raw matrix and trait complete the public S2 workflow without Marine hashes", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor"); skip_if_not_installed("glmnet")
  example <- s2_user_example()
  output <- tempfile("s2-user-")
  on.exit(unlink(output, recursive = TRUE))
  before <- serialize(example, NULL, version = 3)
  result <- do.call(lt_run_s2, c(example, list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", output_dir = output)))
  expect_s3_class(result, "lt_run")
  expect_identical(result$status, "completed")
  expect_identical(serialize(example, NULL, version = 3), before)
  expect_true(is.na(example$x$coordinate_provenance$reference_tree_sha256))
  expect_true(all(file.exists(file.path(output, c("S2_baseline_GBI.rds", "branch_states.tsv",
    "screen_tested.tsv", "OOF_predictions.tsv", "fold_ledger.tsv", "input_hashes.tsv", "SHA256SUMS")))))
  fitted <- readRDS(file.path(output, "validation.rds"))
  expect_true(any(vapply(fitted$folds, function(f) length(f$beta) > 0L, logical(1))))
  expect_identical(fitted$provenance$fold_ledger, example$folds)
  expect_true(all(vapply(fitted$folds, function(f) all(is.finite(f$prediction)), logical(1))))
  expect_equal(nrow(read.delim(file.path(output, "OOF_predictions.tsv"))), 24L)
  expect_error(do.call(lt_run_s2, c(example, list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", output_dir = output))), "exists")
})

test_that("user workflow validates keyed science rather than stale display fingerprints", {
  skip_if_not_installed("ape"); skip_if_not_installed("castor"); skip_if_not_installed("glmnet")
  e <- s2_user_example()
  reordered <- e$coordinates[nrow(e$coordinates):1, ]
  expect_identical(LaTerra:::.lt_s2_user_coordinates(e$x, e$tree, e$coordinates),
    LaTerra:::.lt_s2_user_coordinates(e$x, e$tree, reordered))
  renamed <- e$x; renamed$gene_display_labels <- paste("my label", 1:32)
  renamed$metadata <- list(user_note = "edited without freezing")
  expect_identical(LaTerra:::.lt_s2_user_coordinates(renamed, e$tree, e$coordinates),
    LaTerra:::.lt_s2_user_coordinates(e$x, e$tree, e$coordinates))
  e$coordinates$side_B_taxa[1] <- e$coordinates$side_A_taxa[1]
  expect_error(LaTerra:::.lt_s2_user_coordinates(e$x, e$tree, e$coordinates), "partition")
  e <- s2_user_example(); e$trait$type <- "continuous"
  expect_error(do.call(lt_run_s2, c(e, list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", output_dir = tempfile()))), "binary")
  expect_error(lt_run_s2(recipe = "automatic"), "Explicitly")
  e <- s2_user_example(); e$trait$coding <- list(positive = 0)
  expect_error(do.call(lt_run_s2, c(e, list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", output_dir = tempfile()))), "coding conflicts")
  e <- s2_user_example(); e$tree$edge[2, 2] <- e$tree$edge[1, 2]
  expect_error(LaTerra:::.lt_s2_user_coordinates(e$x, e$tree, e$coordinates), "Malformed")
})
