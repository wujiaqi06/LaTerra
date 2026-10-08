test_that("Marine frozen authority is not replaced by editable inspection copies", {
  original <- lt_marine_profile()
  edited <- original; edited$inputs$traits$sha256 <- paste(rep("0", 64), collapse = "")
  expect_identical(lt_marine_profile(), original)
  expect_false(identical(lt_marine_profile(), edited))
  a <- LaTerra:::.lt_marine_m2_authority()
  replacement <- a
  replacement$inputs$turnover$sha256 <- paste(rep("f", 64), collapse = "")
  replacement$turnover_source_sha256 <- replacement$inputs$turnover$sha256
  expect_identical(LaTerra:::.lt_marine_m2_authority(), a)
  expect_false(identical(LaTerra:::.lt_marine_m2_authority(), replacement))
  expect_false("specification" %in% names(formals(LaTerra:::.lt_marine_m2_inputs)))
  expect_false("profile" %in% names(formals(lt_run_marine)))
})

test_that("historical fold stages retain explicit seed offsets and excluded genera", {
  terminal <- data.frame(genus = c("z", "a", "excluded"), response = c(1, 0, 0.5))
  f <- LaTerra:::.lt_marine_model_folds
  expect_identical(f(terminal, "phase11", 2L)$genus, c("a", "excluded", "z"))
  expect_identical(f(terminal, "phase11", 2L)$seed, 20260523L + 1:3)
  expect_identical(f(terminal, "nested_phase12B", 2L)$seed, 20460523L + 1:3)
  expect_identical(f(terminal, "nested_phase13", 8L)$seed, 21060524L + 1:3)
})

test_that("historical intermediate numeric serialization is an explicit separate operation", {
  x <- c(0.12345678901234567, NA_real_, 1e-30, 1)
  parsed <- LaTerra:::.lt_s2_numeric_tsv_roundtrip(x)
  expect_identical(is.na(parsed), is.na(x))
  expect_equal(parsed[1], 0.123456789012346, tolerance = 0)
  expect_identical(x[1], 0.12345678901234567)
})

test_that("module metrics retain supplied annotations and explicit missing-gene behavior", {
  lookup <- c(a = "Module X", b = "Module X", c = "Module Y", d = "REVIEW_UNANNOTATED")
  m <- LaTerra:::.lt_s2_set_metrics(c("a", "b", "a", "unknown"), c("b", "c", "d"), lookup)
  expect_equal(m$set_A_size_total, 3)
  expect_equal(m$set_A_size_module_annotated, 2)
  expect_equal(m$set_B_size_module_annotated, 2)
  expect_equal(m$gene_overlap_Jaccard, 1/5)
  expect_equal(m$module_presence_Jaccard, 1/2)
  expect_equal(m$module_count_cosine_similarity, 1/sqrt(2))
  expect_identical(LaTerra:::.lt_s2_module_of("unknown", lookup), "REVIEW_UNANNOTATED")
})

test_that("full-data projection uses frozen fitted training transforms", {
  skip_if_not_installed("glmnet")
  x <- rbind(g2 = 1:24/10, g1 = ((1:24) %% 7)/5 + 0.1)
  colnames(x) <- paste0("custom", 1:24)
  terminal <- data.frame(branch = colnames(x), response = rep(0:1, each = 12))
  fit <- LaTerra:::.lt_s2_fit_full(x, terminal, rownames(x), 4201L)
  p <- LaTerra:::.lt_s2_project_fit(x, fit, colnames(x)[1:2])
  expect_equal(p$probability, as.numeric(stats::predict(fit$fit,
    newx = fit$design$x_train[1:2, , drop = FALSE], type = "response")[, 1]))
  expect_identical(p$logit, stats::qlogis(pmin(pmax(p$probability, 1e-12), 1-1e-12)))
  expect_error(LaTerra:::.lt_s2_project_fit(x, list(fit = NULL), colnames(x)[1]), "requires")
})
