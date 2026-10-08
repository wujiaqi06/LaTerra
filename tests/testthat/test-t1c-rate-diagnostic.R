t1c_objects <- function() {
  y <- matrix(c(
    1, 0, 4,
    2, 5, 3,
    0, 0, 0
  ), 3, 3, byrow = TRUE,
  dimnames = list(c("g1", "g2", "DAOA"), c("terminal_1", "internal_1", "terminal_2")))
  state <- matrix("observed", 3, 3, dimnames = dimnames(y))
  x <- lt_matrix(y, state, fixture_payload(),
                 coordinate_provenance = fixture_coordinate_provenance())
  eligibility <- lt_measurement_eligibility(x)
  c3 <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "t1c_c3")
  c2 <- lt_c2_fit(lt_c2_domain(x, eligibility), run_id = "t1c_c2")
  base_rates <- setNames(lapply(c("ADD_LT", "GBI", "logGBI_strict"), function(method) {
    lt_rate(x, c3, lt_rate_spec(paste0("t1c_", method), method), eligibility)
  }), c("ADD_LT", "GBI", "logGBI_strict"))
  c3rates <- c(base_rates[c("ADD_LT", "GBI")],
               list(logGBI = base_rates$logGBI_strict,
                    logGBI_strict = base_rates$logGBI_strict))
  c2rates <- setNames(lapply(c("ADD_C2", "GBI_C2", "logGBI_C2"), function(method) {
    lt_c2_sensitivity(x, c2, method, eligibility)
  }), c("ADD_C2", "GBI_C2", "logGBI_C2"))
  c(list(x = x, eligibility = eligibility, c3 = c3, c2 = c2),
    c3rates, c2rates)
}

test_that("T1C-001 ADD identity is exact on the GBI domain", {
  z <- t1c_objects(); keep <- z$GBI$representation_availability
  mu <- matrix(NA_real_, 3, 3, dimnames = dimnames(z$x$values))
  mu[z$c3$edge_linear] <- z$c3$mu
  expect_equal(z$ADD_LT$values[keep],
               mu[keep] * (z$GBI$values[keep] - 1), tolerance = 1e-14)
})

test_that("T1C-002 positive logGBI equals natural log GBI", {
  z <- t1c_objects(); keep <- z$logGBI$representation_availability
  expect_equal(z$logGBI$values[keep], log(z$GBI$values[keep]), tolerance = 0)
})

test_that("T1C-003 positive GBI equals exp logGBI", {
  z <- t1c_objects(); keep <- z$logGBI$representation_availability
  expect_equal(z$GBI$values[keep], exp(z$logGBI$values[keep]), tolerance = 1e-15)
})

test_that("T1C-004 positive-common sign identities hold", {
  z <- t1c_objects(); keep <- z$logGBI$representation_availability
  expect_identical(sign(z$ADD_LT$values[keep]), sign(z$GBI$values[keep] - 1))
  expect_identical(sign(z$ADD_LT$values[keep]), sign(z$logGBI$values[keep]))
})

test_that("T1C-005 GBI and logGBI ranks are identical with deterministic ties", {
  z <- t1c_objects(); keep <- z$logGBI$representation_availability
  expect_identical(rank(z$GBI$values[keep], ties.method = "min"),
                   rank(z$logGBI$values[keep], ties.method = "min"))
  expect_equal(cor(z$GBI$values[keep], z$logGBI$values[keep],
                   method = "spearman"), 1, tolerance = 1e-15)
})

test_that("T1C-006 C3 ADD gene sums center at zero", {
  z <- t1c_objects()
  expect_lt(max(abs(rowSums(z$ADD_LT$values, na.rm = TRUE))), 1e-10)
})

test_that("T1C-007 C3 ADD branch sums center at zero", {
  z <- t1c_objects()
  expect_lt(max(abs(colSums(z$ADD_LT$values, na.rm = TRUE))), 1e-10)
})

test_that("T1C-008 C3 mu-weighted gene GBI means equal one", {
  z <- t1c_objects(); mu <- matrix(NA_real_, 3, 3, dimnames = dimnames(z$x$values))
  mu[z$c3$edge_linear] <- z$c3$mu
  got <- rowSums(mu * z$GBI$values, na.rm = TRUE) / rowSums(mu, na.rm = TRUE)
  expect_lt(max(abs(got[is.finite(got)] - 1)), 1e-10)
})

test_that("T1C-009 C3 mu-weighted branch GBI means equal one", {
  z <- t1c_objects(); mu <- matrix(NA_real_, 3, 3, dimnames = dimnames(z$x$values))
  mu[z$c3$edge_linear] <- z$c3$mu
  got <- colSums(mu * z$GBI$values, na.rm = TRUE) / colSums(mu, na.rm = TRUE)
  expect_lt(max(abs(got[is.finite(got)] - 1)), 1e-10)
})

test_that("T1C-010 C3 log centers are measured rather than hard-coded zero", {
  z <- t1c_objects()
  observed <- rowMeans(z$logGBI$values, na.rm = TRUE)
  expect_true(any(abs(observed[is.finite(observed)]) > 1e-4))
  d <- lt_diagnose_rate(z$logGBI)
  expect_equal(d$gene_summary$mean, unname(observed))
})

test_that("T1C-011 C2 positive-fit gene log means center at zero", {
  z <- t1c_objects(); m <- z$logGBI_C2$values
  got <- rowMeans(m, na.rm = TRUE)
  expect_lt(max(abs(got[is.finite(got)])), 1e-8)
})

test_that("T1C-012 C2 positive-fit branch log means center at zero", {
  z <- t1c_objects(); m <- z$logGBI_C2$values
  got <- colMeans(m, na.rm = TRUE)
  expect_lt(max(abs(got[is.finite(got)])), 1e-8)
})

test_that("T1C-013 native representation domains remain explicitly distinct", {
  z <- t1c_objects()
  expect_gt(sum(z$ADD_LT$representation_availability),
            sum(z$GBI$representation_availability))
  expect_gt(sum(z$GBI$representation_availability),
            sum(z$logGBI$representation_availability))
  expect_identical(z$logGBI$representation_availability,
                   z$logGBI_strict$representation_availability)
})

test_that("T1C-014 common positive count and keyed hash are exact", {
  z <- t1c_objects(); common <- lt_rate_common_domain(z$ADD_LT, z$GBI, z$logGBI)
  mask <- z$ADD_LT$representation_availability &
    z$GBI$representation_availability & z$logGBI$representation_availability
  expect_identical(common$common_cell_count, sum(mask))
  expect_identical(common$common_cell_sha256,
                   LaTerra:::.lt_rate_keyed_mask_hash(mask, z$x$gene_ids,
                                                      z$x$branch_ids))
})

test_that("T1C-015 zero-loss accounting is exact", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$logGBI)
  expect_identical(d$zero_loss$input_zero_log_undefined_count, 1L)
  expect_identical(d$zero_loss$baseline_zero_count, 3L)
  expect_identical(d$zero_loss$available_count, 5L)
  expect_identical(d$zero_loss$zero_encoded_count, 0L)
  strict <- lt_diagnose_rate(z$logGBI_strict)
  expect_identical(strict$zero_loss$input_zero_log_undefined_count, 1L)
  expect_identical(strict$zero_loss$zero_encoded_count, 0L)
  expect_identical(strict$zero_loss$available_count, 5L)
})

test_that("T1C-016 primary gene raw scale is arithmetic mean including zeros", {
  z <- t1c_objects()
  d <- lt_diagnose_rate(z$ADD_LT, z$x, z$c3, z$eligibility)
  expect_equal(d$gene_summary$raw_mean, unname(rowMeans(z$x$values)))
  expect_false(identical(d$gene_summary$raw_mean, apply(z$x$values, 1, median)))
})

test_that("T1C-017 DAOA zero-margin boundary stays visible", {
  z <- t1c_objects()
  d <- lt_diagnose_rate(z$logGBI, z$x, z$c3, z$eligibility)
  daoa <- d$gene_summary[d$gene_summary$scientific_key == "DAOA", ]
  expect_identical(daoa$raw_zero_fraction, 1)
  expect_identical(daoa$available_count, 0)
  expect_identical(daoa$baseline_mu_max, 0)
})

test_that("T1C-018 HQ median reference cannot change C3 mu or rate", {
  z <- t1c_objects(); before <- z$GBI$output_hashes$semantic_object_sha256
  d <- lt_diagnose_rate(z$GBI)
  d$metadata$branch_reference <- "HQ2275_raw_median"
  expect_invisible(validate_lt_rate_diagnostic(d))
  expect_identical(z$GBI$output_hashes$semantic_object_sha256, before)
})

test_that("T1C-019 raw mean reference cannot change C3 mu or rate", {
  z <- t1c_objects(); before <- z$ADD_LT$values
  d <- lt_diagnose_rate(z$ADD_LT)
  d$metadata$branch_reference <- "full_raw_mean"
  expect_invisible(validate_lt_rate_diagnostic(d))
  expect_identical(z$ADD_LT$values, before)
})

test_that("T1C-020 bounded extreme identities are exact and unique", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$ADD_LT, extreme_n = 3)
  expect_identical(anyDuplicated(d$extreme_ledger[c("rank_type", "rank")]), 0L)
  i <- d$extreme_ledger$linear_index
  expect_identical(d$extreme_ledger$value, z$ADD_LT$values[i])
  expect_true(all(d$extreme_ledger$gene_id %in% z$x$gene_ids))
  expect_true(all(d$extreme_ledger$branch_id %in% z$x$branch_ids))
})

test_that("T1C-021 diagnosis creates no exclusions", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$GBI)
  expect_true(d$provenance$diagnose_only)
  expect_identical(d$provenance$exclusion_count, 0L)
  expect_false(d$provenance$trimming)
  expect_identical(z$GBI$measurement_fit_eligibility$status,
                   z$eligibility$status)
})

test_that("T1C-022 no normal-theory outlier classification is present", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$ADD_LT)
  expect_false(d$provenance$normality_rule)
  expect_false(d$provenance$outlier_classification)
  expect_false(any(grepl("outlier", names(d$extreme_ledger), ignore.case = TRUE)))
})

test_that("T1C-023 C2 remains explicit sensitivity only", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$GBI_C2, z$x, z$c2, z$eligibility)
  expect_identical(d$baseline_estimand, "C2_positive_cell_log_OLS")
  expect_identical(d$compatibility$state, "CURRENT")
  expect_false(isTRUE(z$GBI_C2$metadata$production))
  expect_true(isTRUE(z$GBI_C2$metadata$sensitivity))
})

test_that("T1C-024 a trait object cannot enter rate diagnosis", {
  z <- t1c_objects()
  expect_false("trait" %in% names(formals(lt_diagnose_rate)))
  expect_error(do.call(lt_diagnose_rate, list(rate = z$GBI, trait = "aquatic")),
               "unused argument")
})

test_that("T1C-025 RER is absent from the diagnostic operator", {
  expect_false("rer" %in% tolower(names(formals(lt_diagnose_rate))))
  expect_false("lt_rer" %in% getNamespaceExports("LaTerra"))
  desc <- packageDescription("LaTerra")
  expect_false(grepl("RERconverge", paste(unclass(desc), collapse = " "),
                     fixed = TRUE))
})

test_that("T1C-026 rate diagnostic round trip is exact", {
  z <- t1c_objects(); d <- lt_diagnose_rate(z$logGBI, z$x, z$c3, z$eligibility)
  path <- tempfile(fileext = ".rds")
  lt_write_object(d, path)
  expect_identical(lt_read_object(path), d)
})

test_that("T1C-027 compatibility mismatch is visible and never repaired", {
  z <- t1c_objects(); edited <- z$x
  edited$values[1, 1] <- edited$values[1, 1] + .5
  expect_invisible(validate_lt_matrix(edited))
  expect_error(lt_diagnose_rate(z$ADD_LT, edited, z$c3, z$eligibility),
               "STALE_RECOMPUTE_REQUIRED.*no automatic repair")
})

test_that("T1C-028 ordinary runtime diagnosis needs no frozen-file SHA", {
  y <- matrix(c(1, 2, 3, 4), 2, 2,
              dimnames = list(c("user_g1", "user_g2"), c("user_b1", "user_b2")))
  state <- matrix("observed", 2, 2, dimnames = dimnames(y))
  x <- lt_matrix(y, state, fixture_payload(), coordinate_provenance = NULL)
  e <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, e), run_id = "ordinary_runtime_t1c")
  rate <- lt_rate(x, fit, lt_rate_spec("ordinary_runtime_gbi", "GBI"), e)
  d <- lt_diagnose_rate(rate, x, fit, e)
  expect_invisible(validate_lt_rate_diagnostic(d))
  expect_true(is.na(x$coordinate_provenance$reference_tree_sha256))
  expect_false("frozen_file_sha256" %in% names(d))
})
