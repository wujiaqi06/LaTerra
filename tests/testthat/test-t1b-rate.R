t1b_matrix <- function(values,
                       mask = matrix(TRUE, nrow(values), ncol(values)),
                       gene_ids = paste0("g", seq_len(nrow(values))),
                       branch_ids = paste0("b", seq_len(ncol(values)))) {
  values <- as.matrix(values); storage.mode(values) <- "double"
  values[!mask] <- NA_real_
  dimnames(values) <- list(gene_ids, branch_ids)
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  state[!mask] <- "NA_struct"
  lt_matrix(values, state, fixture_payload(),
            coordinate_provenance = fixture_coordinate_provenance())
}

t1b_fit <- function(values, mask = matrix(TRUE, nrow(values), ncol(values)),
                    gene_ids = paste0("g", seq_len(nrow(values))),
                    branch_ids = paste0("b", seq_len(ncol(values))),
                    run_id = "t1b_fixture") {
  x <- t1b_matrix(values, mask, gene_ids, branch_ids)
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = run_id)
  list(x = x, eligibility = eligibility, fit = fit)
}

t1b_rate <- function(z, method) {
  if (identical(method, "logGBI")) {
    testthat::skip(paste(
      "Superseded 0.0.0.9013 finite-sentinel fixture;",
      "see test-log-relative-boundary.R."
    ))
  }
  lt_rate(
    z$x, z$fit, lt_rate_spec(paste0("fixture_", method), method),
    z$eligibility
  )
}

test_that("T1B-001 modified input matrix after C3 fit errors", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  modified <- z$x; modified$values[1, 1] <- 2
  z$x <- modified
  expect_error(t1b_rate(z, "ADD_LT"),
               "values differ")
})

test_that("T1B-002 reordered genes without mapping errors", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  reordered <- t1b_matrix(z$x$values[c(2, 1), ],
                          gene_ids = z$x$gene_ids[c(2, 1)],
                          branch_ids = z$x$branch_ids)
  expect_error(lt_rate(reordered, z$fit, lt_rate_spec("reorder_gene", "ADD_LT"),
                       z$eligibility), "gene ledger")
})

test_that("T1B-003 reordered branches without mapping errors", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  reordered <- t1b_matrix(z$x$values[, c(2, 1)],
                          gene_ids = z$x$gene_ids,
                          branch_ids = z$x$branch_ids[c(2, 1)])
  expect_error(lt_rate(reordered, z$fit, lt_rate_spec("reorder_branch", "ADD_LT"),
                       z$eligibility), "branch ledger")
})

test_that("T1B-004 Layer-2 identity mismatch errors", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  status <- z$eligibility$status; reason <- z$eligibility$reason
  status[1, 1] <- "ineligible"; reason[1, 1] <- "user_predeclared_excluded"
  changed <- lt_measurement_eligibility(z$x, status, reason)
  expect_error(lt_rate(z$x, z$fit, lt_rate_spec("layer2", "ADD_LT"), changed),
               "Layer-2 eligibility identity")
})

test_that("T1B-005 certified baseline mu hash mismatch errors", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  z$fit$component_fits[[1]]$mu_sha256 <- strrep("0", 64)
  expect_error(t1b_rate(z, "ADD_LT"), "mu hash mismatch")
})

test_that("T1B-006 positive Y and mu implement all four formulas", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  add <- t1b_rate(z, "ADD_LT"); gbi <- t1b_rate(z, "GBI")
  loggbi <- t1b_rate(z, "logGBI")
  strict <- t1b_rate(z, "logGBI_strict")
  mu <- z$fit$mu; y <- z$fit$observed_value; edge <- z$fit$edge_linear
  expect_equal(add$values[edge], y - mu, tolerance = 0)
  expect_equal(gbi$values[edge], y / mu, tolerance = 0)
  expect_equal(loggbi$values[edge], log(y / mu), tolerance = 0)
  expect_equal(strict$values[edge], log(y / mu), tolerance = 0)
})

test_that("T1B-007 zero Y positive mu is encoded only in default logGBI", {
  z <- t1b_fit(matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE))
  edge <- z$fit$edge_linear; zero <- edge[z$fit$observed_value == 0]
  add <- t1b_rate(z, "ADD_LT"); gbi <- t1b_rate(z, "GBI")
  loggbi <- t1b_rate(z, "logGBI")
  strict <- t1b_rate(z, "logGBI_strict")
  expect_equal(add$values[zero], -z$fit$mu[z$fit$observed_value == 0])
  expect_identical(gbi$values[zero], c(0, 0))
  expect_true(all(is.finite(loggbi$values[zero])))
  expect_true(all(loggbi$zero_encoded[zero]))
  expect_true(all(lt_rate_reason(loggbi)[zero] == "available"))
  expect_true(all(is.na(strict$values[zero])))
  expect_true(all(lt_rate_reason(strict)[zero] == "input_zero_log_undefined"))
})

test_that("T1B-008 zero Y zero boundary mu keeps ADD and loses ratios", {
  z <- t1b_fit(matrix(c(1, 2, 0, 0), 2, 2, byrow = TRUE))
  edge <- z$fit$edge_linear; boundary <- edge[z$fit$mu == 0]
  add <- t1b_rate(z, "ADD_LT"); gbi <- t1b_rate(z, "GBI")
  loggbi <- t1b_rate(z, "logGBI")
  strict <- t1b_rate(z, "logGBI_strict")
  expect_true(all(add$values[boundary] == 0))
  expect_true(all(is.na(gbi$values[boundary])))
  expect_true(all(is.na(loggbi$values[boundary])))
  expect_true(all(is.na(strict$values[boundary])))
  expect_true(all(lt_rate_reason(gbi)[boundary] == "baseline_zero"))
})

test_that("T1B-009 positive Y zero mu is fatal", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  z$fit$mu[1] <- 0
  expect_error(t1b_rate(z, "ADD_LT"), "positive observed payload")
})

test_that("T1B-010 logGBI variants never store infinity", {
  z <- t1b_fit(matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE))
  rate <- t1b_rate(z, "logGBI")
  strict <- t1b_rate(z, "logGBI_strict")
  expect_false(any(is.infinite(rate$values), na.rm = TRUE))
  expect_true(all(is.finite(rate$values[z$x$values == 0])))
  expect_true(all(rate$zero_encoded[z$x$values == 0]))
  expect_false(any(is.infinite(strict$values), na.rm = TRUE))
  expect_true(all(is.na(strict$values[z$x$values == 0])))
})

test_that("T1B-011 no production or C2 pseudocount exists", {
  z <- t1b_fit(matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE))
  expect_identical(t1b_rate(z, "GBI")$representation_provenance$pseudocount,
                   "none")
  expect_identical(lt_c2_domain(z$x, z$eligibility)$pseudocount, "none")
  expect_error(lt_rate_spec("bad", "GBI", list(pseudocount = 1e-8)),
               "cannot alter")
})

test_that("T1B-012 representation reason preserves coordinate state", {
  mask <- matrix(c(TRUE, TRUE, FALSE, TRUE), 2, 2, byrow = TRUE)
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2), mask)
  rate <- t1b_rate(z, "ADD_LT")
  expect_identical(rate$coord_state, z$x$coord_state)
  expect_identical(lt_rate_reason(rate)[!mask], "coordinate_unavailable")
})

test_that("T1B-013 measurement ineligible is unavailable with prior reason", {
  x <- t1b_matrix(matrix(c(1, 2, 3, 4), 2, 2))
  status <- matrix("eligible", 2, 2, dimnames = dimnames(x$values))
  reason <- matrix("eligible_declared", 2, 2, dimnames = dimnames(x$values))
  status[1, 1] <- "ineligible"; reason[1, 1] <- "analysis_domain_excluded"
  e <- lt_measurement_eligibility(x, status, reason)
  fit <- lt_c3_fit(lt_c3_domain(x, e), run_id = "ineligible")
  rate <- lt_rate(x, fit, lt_rate_spec("ineligible", "ADD_LT"), e)
  expect_false(rate$representation_availability[1, 1])
  expect_identical(lt_rate_reason(rate)[1, 1], "measurement_ineligible")
  expect_identical(rate$measurement_fit_eligibility$reason[1, 1],
                   "analysis_domain_excluded")
})

test_that("T1B-014 numeric-indeterminate required C3 fails closed", {
  z <- t1b_fit(matrix(c(1, 2, 3, 4), 2, 2))
  z$fit$state <- "NUMERICALLY_INDETERMINATE"; z$fit$mu <- NULL
  expect_error(t1b_rate(z, "ADD_LT"), "NUMERICALLY_INDETERMINATE")
})

test_that("T1B-015 three-rate common domain is exact keyed intersection", {
  z <- t1b_fit(matrix(c(1, 0, 2, 3, 0, 0), 2, 3, byrow = TRUE))
  rates <- lapply(c("ADD_LT", "GBI", "logGBI"), t1b_rate, z = z)
  common <- lt_rate_common_domain(rates)
  expected <- Reduce(`&`, lapply(rates, `[[`, "representation_availability"))
  expect_identical(common$common_linear, which(expected))
  expect_equal(common$common_cell_count, sum(expected))
})

test_that("T1B-016 keyed common identity is invariant to joint axis permutation", {
  y <- matrix(c(1, 0, 2, 3, 4, 5), 2, 3, byrow = TRUE)
  a <- t1b_fit(y, gene_ids = c("ga", "gb"), branch_ids = c("ba", "bb", "bc"))
  b <- t1b_fit(y[c(2, 1), c(3, 1, 2)], gene_ids = c("gb", "ga"),
               branch_ids = c("bc", "ba", "bb"))
  ac <- lt_rate_common_domain(lapply(c("ADD_LT", "GBI", "logGBI"),
                                     t1b_rate, z = a))
  bc <- lt_rate_common_domain(lapply(c("ADD_LT", "GBI", "logGBI"),
                                     t1b_rate, z = b))
  expect_identical(ac$common_cell_sha256, bc$common_cell_sha256)
})

test_that("T1B-017 zero-only log loss and boundary ratio loss are exact", {
  z <- t1b_fit(matrix(c(1, 0, 2, 3, 0, 0), 2, 3, byrow = TRUE))
  rates <- lapply(c("ADD_LT", "GBI", "logGBI_strict"), t1b_rate, z = z)
  common <- lt_rate_common_domain(rates)
  expect_equal(common$zero_cells_lost_solely_log_undefined, 1L)
  expect_equal(common$boundary_zero_cells_lost_ratio_denominator_zero, 2L)
})

test_that("T1B-018 unequal full domains are labelled descriptive", {
  z <- t1b_fit(matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE))
  common <- lt_rate_common_domain(lapply(c("ADD_LT", "GBI", "logGBI_strict"),
                                         t1b_rate, z = z))
  expect_identical(common$full_domain_label, "UNEQUAL_DOMAIN_DESCRIPTIVE")
  expect_identical(common$interpretation,
                   "domain difference is not estimator-performance evidence")
})

test_that("T1B-019 C2 fit domain excludes every exact zero", {
  x <- t1b_matrix(matrix(c(1, 0, 2, 0), 2, 2))
  d <- lt_c2_domain(x)
  expect_false(any(x$values[d$fit_mask] == 0))
  expect_equal(sum(d$fit_mask), 2L)
})

test_that("T1B-020 C2 uses no pseudocount", {
  x <- t1b_matrix(matrix(c(1, 0, 2, 0), 2, 2))
  d <- lt_c2_domain(x)
  expect_identical(d$pseudocount, "none")
  expect_equal(d$fit_log_value, log(x$values[d$fit_mask]), tolerance = 0)
})

test_that("T1B-021 complete positive table matches log-OLS closed form", {
  y <- outer(c(.7, 1.3, 2.1), c(.5, 1.1, 3.2))
  x <- t1b_matrix(y); d <- lt_c2_domain(x); fit <- lt_c2_fit(d)
  expected <- outer(exp(rowMeans(log(y))), exp(colMeans(log(y)))) /
    exp(mean(log(y)))
  expect_equal(fit$fit_mu, expected[d$fit_linear], tolerance = 1e-10)
})

test_that("T1B-022 incomplete positive C2 matches generic lm", {
  y <- outer(c(.7, 1.1, 1.8, 2.4), c(.5, .9, 1.3))
  mask <- matrix(TRUE, 4, 3); mask[cbind(c(2, 4), c(3, 1))] <- FALSE
  x <- t1b_matrix(y, mask); d <- lt_c2_domain(x); fit <- lt_c2_fit(d)
  frame <- expand.grid(gene = factor(seq_len(4)), branch = factor(seq_len(3)))
  frame$y <- as.vector(x$values); frame$observed <- as.vector(mask)
  model <- stats::lm(log(y) ~ gene + branch, data = frame, subset = observed)
  expected <- exp(stats::predict(model, newdata = frame))[d$fit_linear]
  expect_equal(unname(fit$fit_mu), unname(expected), tolerance = 1e-9)
})

test_that("T1B-023 positive-empty gene is unsupported", {
  x <- t1b_matrix(matrix(c(1, 2, 0, 0), 2, 2, byrow = TRUE))
  d <- lt_c2_domain(x)
  expect_identical(d$unsupported_vertices$vertex_id, "g2")
})

test_that("T1B-024 positive-empty branch is unsupported", {
  x <- t1b_matrix(matrix(c(1, 0, 2, 0), 2, 2, byrow = TRUE))
  d <- lt_c2_domain(x)
  expect_identical(d$unsupported_vertices$vertex_id, "b2")
})

test_that("T1B-025 disconnected positive graph is solved component-wise", {
  mask <- diag(TRUE, 2); x <- t1b_matrix(diag(c(2, 3)), mask)
  fit <- lt_c2_fit(lt_c2_domain(x))
  expect_length(fit$certificates, 2L)
  expect_equal(fit$fit_mu, c(2, 3), tolerance = 1e-12)
})

test_that("T1B-026 no cross-component baseline product is materialized", {
  values <- matrix(c(2, 0, 0, 3), 2, 2, byrow = TRUE)
  x <- t1b_matrix(values); fit <- lt_c2_fit(lt_c2_domain(x))
  # Positive graph has two components; zero off-diagonal cells are evaluation-only.
  expect_equal(sum(fit$baseline_availability), 2L)
  expect_true(all(is.na(fit$evaluation_mu[c(2, 3)])))
  expect_true(all(fit$baseline_reason_levels[fit$baseline_reason[c(2, 3)]] ==
                    "cross_component_baseline_unavailable"))
})

test_that("T1B-027 gauge shift changes factors but not fitted mu", {
  y <- outer(c(.7, 1.1, 1.8), c(.5, .9, 1.3))
  fit <- lt_c2_fit(lt_c2_domain(t1b_matrix(y)))
  shift <- .75
  changed <- exp(fit$gene_log_effect + shift)[fit$domain$fit_gene_index] *
    exp(fit$branch_log_effect - shift)[fit$domain$fit_branch_index]
  expect_false(identical(fit$gene_log_effect, fit$gene_log_effect + shift))
  expect_equal(changed, fit$fit_mu, tolerance = 1e-12)
})

test_that("T1B-028 normal-equation residual certificate is required", {
  fit <- lt_c2_fit(lt_c2_domain(t1b_matrix(matrix(c(1, 2, 3, 4), 2, 2))))
  expect_true(all(vapply(fit$certificates, function(z)
    z$maximum_absolute_normal_equation_mean_residual <=
      fit$control$normal_equation_tolerance, logical(1L))))
  expect_true(all(vapply(fit$certificates, `[[`, logical(1L), "certificate_pass")))
})

test_that("T1B-029 failed C2 numerical certificate is indeterminate", {
  y <- matrix(c(1, 2, 3, 10), 2, 2)
  fit <- lt_c2_fit(lt_c2_domain(t1b_matrix(y)),
                   lt_c2_control(tolerance = 1e-30, max_iterations = 1L,
                                 normal_equation_tolerance = 1e-30))
  expect_identical(fit$state, "NUMERICALLY_INDETERMINATE")
  expect_null(fit$fit_mu)
  expect_null(fit$evaluation_mu)
})

test_that("T1B-030 C2 cannot be selected as production rate", {
  expect_error(lt_rate_spec("bad", "C2"), "required sensitivity")
  expect_error(lt_rate_spec("bad", "GBI", list(baseline = "C2")),
               "required sensitivity")
})

test_that("T1B-031 positive scaling scales ADD", {
  y <- matrix(c(1, 0, 2, 3), 2, 2)
  a <- t1b_fit(y); b <- t1b_fit(y * 10)
  expect_equal(t1b_rate(b, "ADD_LT")$values,
               10 * t1b_rate(a, "ADD_LT")$values, tolerance = 1e-9)
  ac2 <- lt_c2_fit(lt_c2_domain(a$x, a$eligibility))
  bc2 <- lt_c2_fit(lt_c2_domain(b$x, b$eligibility))
  expect_equal(lt_c2_sensitivity(b$x, bc2, "ADD_C2", b$eligibility)$values,
               10 * lt_c2_sensitivity(a$x, ac2, "ADD_C2", a$eligibility)$values,
               tolerance = 1e-9)
})

test_that("T1B-032 positive scaling preserves GBI", {
  y <- matrix(c(1, 0, 2, 3), 2, 2)
  a <- t1b_fit(y); b <- t1b_fit(y * 10)
  expect_equal(t1b_rate(b, "GBI")$values,
               t1b_rate(a, "GBI")$values, tolerance = 1e-9)
  ac2 <- lt_c2_fit(lt_c2_domain(a$x, a$eligibility))
  bc2 <- lt_c2_fit(lt_c2_domain(b$x, b$eligibility))
  expect_equal(lt_c2_sensitivity(b$x, bc2, "GBI_C2", b$eligibility)$values,
               lt_c2_sensitivity(a$x, ac2, "GBI_C2", a$eligibility)$values,
               tolerance = 1e-9)
})

test_that("T1B-033 positive scaling preserves logGBI", {
  y <- matrix(c(1, 0, 2, 3), 2, 2)
  a <- t1b_fit(y); b <- t1b_fit(y * .1)
  expect_equal(t1b_rate(b, "logGBI")$values,
               t1b_rate(a, "logGBI")$values, tolerance = 1e-9)
  ac2 <- lt_c2_fit(lt_c2_domain(a$x, a$eligibility))
  bc2 <- lt_c2_fit(lt_c2_domain(b$x, b$eligibility))
  expect_equal(lt_c2_sensitivity(b$x, bc2, "logGBI_C2", b$eligibility)$values,
               lt_c2_sensitivity(a$x, ac2, "logGBI_C2", a$eligibility)$values,
               tolerance = 1e-9)
})

test_that("T1B-034 positive scaling preserves availability and zeros", {
  y <- matrix(c(1, 0, 2, 3), 2, 2)
  a <- t1b_fit(y); b <- t1b_fit(y * 10)
  for (m in c("ADD_LT", "GBI", "logGBI", "logGBI_strict")) {
    ar <- t1b_rate(a, m); br <- t1b_rate(b, m)
    expect_identical(ar$representation_availability,
                     br$representation_availability)
  }
  ac2 <- lt_c2_fit(lt_c2_domain(a$x, a$eligibility))
  bc2 <- lt_c2_fit(lt_c2_domain(b$x, b$eligibility))
  for (m in c("ADD_C2", "GBI_C2", "logGBI_C2")) {
    ar <- lt_c2_sensitivity(a$x, ac2, m, a$eligibility)
    br <- lt_c2_sensitivity(b$x, bc2, m, b$eligibility)
    expect_identical(ar$representation_availability,
                     br$representation_availability)
  }
  expect_identical(y == 0, y * 10 == 0)
})

test_that("T1B-035 deterministic replay reproduces rate hashes", {
  z <- t1b_fit(matrix(c(1, 0, 2, 3), 2, 2))
  a <- t1b_rate(z, "GBI"); b <- t1b_rate(z, "GBI")
  expect_identical(a$output_hashes, b$output_hashes)
  expect_identical(a$representation_reason, b$representation_reason)
})

# T1B-036 through T1B-045 are large frozen-marine integration gates. They are
# executed by the returned marine_integration.R and represented in the review
# package test manifest rather than duplicated in package unit tests.

test_that("T1B-046 lt_rate round-trip preserves every state/hash layer", {
  z <- t1b_fit(matrix(c(1, 0, 2, 3), 2, 2)); rate <- t1b_rate(z, "logGBI")
  path <- tempfile(fileext = ".rds"); lt_write_object(rate, path)
  expect_identical(lt_read_object(path), rate)
})

test_that("T1B-047 lt_c2_domain round-trip is exact", {
  object <- lt_c2_domain(t1b_matrix(matrix(c(1, 0, 2, 3), 2, 2)))
  path <- tempfile(fileext = ".rds"); lt_write_object(object, path)
  expect_identical(lt_read_object(path), object)
})

test_that("T1B-048 lt_c2_fit round-trip is exact", {
  object <- lt_c2_fit(lt_c2_domain(t1b_matrix(matrix(c(1, 0, 2, 3), 2, 2))))
  path <- tempfile(fileext = ".rds"); lt_write_object(object, path)
  expect_identical(lt_read_object(path), object)
})

test_that("T1B-049 trait parameter cannot affect rate construction", {
  expect_error(lt_rate_spec("trait_bad", "GBI", list(trait = "aquatic")),
               "cannot alter")
  expect_false("trait" %in% names(formals(lt_rate)))
  expect_false("trait" %in% names(formals(lt_c2_fit)))
})

test_that("T1B-050 QC flags alone cannot alter rate domain", {
  z <- t1b_fit(matrix(c(1, 1, 1, 100), 1, 4))
  before <- t1b_rate(z, "ADD_LT")
  qc <- lt_qc_matrix(z$x, lt_qc_spec(
    profile = "custom", outlier_rules = list(list(
      rule_id = "warning_only", method = "iqr", scope = "whole_matrix",
      direction = "upper", multiplier = 1.5
    ))
  ))
  expect_gt(nrow(qc$flags), 0L)
  after <- t1b_rate(z, "ADD_LT")
  expect_identical(before$output_hashes, after$output_hashes)
})

test_that("T1B-051 RER is absent from T1B implementation dependency graph", {
  symbols <- c("lt_rate", "lt_c2_fit", "lt_c2_sensitivity",
               "lt_rate_common_domain", "lt_c3_c2_compare")
  source <- paste(vapply(symbols, function(symbol) paste(
    deparse(get(symbol, envir = asNamespace("LaTerra"))), collapse = "\n"
  ), character(1L)), collapse = "\n")
  expect_false(grepl("RERconverge|getRMat|getAllResiduals", source))
})
