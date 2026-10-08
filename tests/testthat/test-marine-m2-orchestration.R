test_that("M2 public controls reject invalid workers without creating an output", {
  out <- tempfile()
  for (workers in list(NA_real_, Inf, 0, 1.2, "2", c(1, 2)))
    expect_error(lt_run_marine_m2("not-an-input", out, cores = workers), "positive integer")
  expect_false(file.exists(out))
  expect_false("profile" %in% names(formals(lt_run_marine_m2)))
  expect_false("fitted_inputs" %in% names(formals(lt_run_marine_m2)))
  expect_error(lt_run_marine_m2("not-an-input", out), class = "lt_error_marine_missing_input")
  expect_false(file.exists(out))
})

test_that("full fold provenance retains excluded and training species explicitly", {
  terminal <- data.frame(species = c("Own_a", "Other_b", "Excluded_c"), genus = c("Own", "Other", "Excluded"),
    response = c(1, 0, 0.5))
  x <- list(folds = list(list(identity = list(fold_id = 7L, test_species = "Own_a"))),
    provenance = list(terminal_ledger = terminal, fold_ledger = data.frame(fold_id = 7L, seed = 42L),
      RNGkind = c("Mersenne-Twister", "Inversion", "Rejection")))
  p <- LaTerra:::.lt_s2_fold_provenance(x, "supplied_example")
  expect_identical(p$folds$role, c("test", "train", "excluded"))
  expect_identical(p$folds$taxon_id, terminal$species)
  expect_identical(p$seeds$seed, 42L)
  expect_identical(p$seeds$scope, "supplied_example::7")
})

test_that("module display joins computed sets without using annotation coefficients", {
  full <- list(models = list(fix_marine_binary = list(beta = c(g1 = 1, g2 = 0)),
    fix_aquatic_v2 = list(beta = c(g3 = -1, g1 = 2))))
  annotation <- data.frame(gene = c("g3", "g1"), display_module = c("module", "unassigned"),
    counted_in_Figure4C_circle_size = c("True", "False"), marine_coef = c(999, -999))
  x <- LaTerra:::.lt_marine_module_display(annotation, full)
  expect_identical(x$annotation, annotation)
  expect_identical(x$counts$value, c(2L, 1L, 1L, 1L))
  expect_false(x$provenance$source_numeric_annotation_columns_used_for_fitting)
  annotation$gene[1] <- "unselected"
  expect_error(LaTerra:::.lt_marine_module_display(annotation, full), "union")
})

test_that("sensitivity display cannot substitute a missing or reordered model", {
  expect_error(LaTerra:::.lt_marine_sensitivity_tables(list()), "eight")
  expect_length(LaTerra:::.lt_marine_sensitivity_ids(), 8L)
  expect_identical(LaTerra:::.lt_marine_sensitivity_ids()[2], "fix_drop_whale")
})
