t2a_matrix <- function(values, mask = matrix(TRUE, nrow(values), ncol(values))) {
  values <- as.matrix(values); storage.mode(values) <- "double"
  values[!mask] <- NA_real_
  if (is.null(rownames(values))) rownames(values) <- paste0("g", seq_len(nrow(values)))
  if (is.null(colnames(values))) colnames(values) <- paste0("b", seq_len(ncol(values)))
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  state[!mask] <- "NA_struct"
  lt_matrix(values, state, fixture_payload(),
            coordinate_provenance = fixture_coordinate_provenance())
}

t2a_rates <- function(values, mask = matrix(TRUE, nrow(values), ncol(values))) {
  x <- t2a_matrix(values, mask)
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "t2a_fixture")
  rates <- lapply(c("ADD_LT", "GBI", "logGBI_strict"), function(method) {
    lt_rate(x, fit, lt_rate_spec(paste0("t2a_", method), method), eligibility)
  })
  names(rates) <- c("ADD_LT", "GBI", "logGBI")
  list(x = x, eligibility = eligibility, fit = fit, rates = rates)
}

t2a_continuous_state <- function(branch_ids = paste0("b", 1:4),
                                 values = 0:3,
                                 availability = rep("available", 4),
                                 reason = rep(NA_character_, 4)) {
  lt_branch_state(
    "continuous_fixture", "continuous_trait", branch_ids, values,
    "continuous", availability, reason,
    provenance = list(source_type = "deterministic_fixture")
  )
}

t2a_binary_state <- function(branch_ids = paste0("b", 1:4),
                             values = c("reference", "reference", "focal", "focal"),
                             availability = rep("available", 4),
                             reason = rep(NA_character_, 4),
                             reference = "reference", focal = "focal") {
  lt_branch_state(
    "binary_fixture", "binary_trait", branch_ids, values, "binary",
    availability, reason,
    coding = list(reference_level = reference, focal_level = focal),
    provenance = list(source_type = "deterministic_fixture")
  )
}

t2a_exact_values <- function() {
  matrix(
    c(1, 3, 5, 7, 7, 5, 3, 1), 2, 4, byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("b", 1:4))
  )
}

test_that("T2A-001 lt_trait remains taxon-level and unchanged", {
  trait <- lt_trait("taxon_trait", c("taxon_a", "taxon_b"), c(0, 1), "binary")
  expect_identical(names(trait),
                   c("trait_id", "taxon_ids", "values", "type", "coding", "metadata"))
  expect_false("branch_ids" %in% names(trait))
  expect_invisible(validate_lt_trait(trait))
})

test_that("T2A-002 lt_branch_state exact constructor and validator", {
  x <- t2a_continuous_state()
  expect_s3_class(x, "lt_branch_state")
  expect_identical(x$branch_ids, paste0("b", 1:4))
  expect_invisible(validate_lt_branch_state(x))
})

test_that("T2A-003 available state requires a valid value", {
  expect_error(t2a_continuous_state(values = c(0, 1, NA, 3)),
               "Available branch states require")
  expect_error(t2a_continuous_state(values = c(0, 1, Inf, 3)), "finite")
})

test_that("T2A-004 unavailable state requires NA and reason", {
  expect_error(t2a_continuous_state(
    availability = c("available", "unavailable", "available", "available"),
    reason = c(NA, "missing", NA, NA)), "require an NA value")
  expect_error(t2a_continuous_state(
    values = c(0, NA, 2, 3),
    availability = c("available", "unavailable", "available", "available")),
    "non-empty reason")
})

test_that("T2A-005 unresolved state requires NA and reason", {
  x <- t2a_continuous_state(
    values = c(0, NA, 2, 3),
    availability = c("available", "unresolved", "available", "available"),
    reason = c(NA, "authority_unresolved", NA, NA))
  expect_invisible(validate_lt_branch_state(x))
})

test_that("T2A-006 duplicate branch IDs fail", {
  expect_error(t2a_continuous_state(branch_ids = c("b1", "b1", "b3", "b4")),
               "duplicate")
})

test_that("T2A-007 binary coding requires explicit focal/reference semantics", {
  expect_error(lt_branch_state(
    "bad", "trait", paste0("b", 1:4), c(0, 0, 1, 1), "binary",
    rep("available", 4), rep(NA_character_, 4)), "explicit reference_level")
})

test_that("T2A-008 binary lexical order cannot change effect direction", {
  z <- t2a_rates(t2a_exact_values())
  state <- t2a_binary_state(
    values = c("z_reference", "z_reference", "a_focal", "a_focal"),
    reference = "z_reference", focal = "a_focal")
  a <- lt_associate_univariate(z$rates$ADD_LT, state, mismatch = "error")
  expect_equal(a$estimates$beta, c(4, -4), tolerance = 0)
  expect_identical(a$effect_direction, "a_focal - z_reference")
})

test_that("T2A-009 continuous exact OLS slope and intercept fixture", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_equal(a$estimates$beta, c(2, -2), tolerance = 0)
  expect_equal(a$estimates$intercept, c(-3, 3), tolerance = 0)
})

test_that("T2A-010 continuous Pearson is exact", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_equal(a$descriptive_scores$pearson_r, c(1, -1), tolerance = 0)
  expect_equal(a$descriptive_scores$r_squared, c(1, 1), tolerance = 0)
})

test_that("T2A-011 continuous Spearman is exact", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_equal(a$descriptive_scores$spearman_rho, c(1, -1), tolerance = 0)
})

test_that("T2A-012 binary beta equals focal mean minus reference mean", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_binary_state(),
                               mismatch = "error")
  expect_equal(a$estimates$reference_mean, c(-2, 2), tolerance = 0)
  expect_equal(a$estimates$focal_mean, c(2, -2), tolerance = 0)
  expect_equal(a$estimates$beta,
               a$estimates$focal_mean - a$estimates$reference_mean,
               tolerance = 0)
})

test_that("T2A-013 binary point-biserial is exact", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_binary_state(),
                               mismatch = "error")
  expected <- vapply(seq_len(2), function(i) {
    cor(c(0, 0, 1, 1), z$rates$ADD_LT$values[i, ])
  }, numeric(1L))
  expect_equal(a$descriptive_scores$point_biserial_r, expected,
               tolerance = .Machine$double.eps)

  reference_values <- matrix(
    c(2, 7, 5, 11, 3, 13,
      17, 4, 19, 6, 23, 8,
      9, 29, 10, 31, 12, 37),
    3, 6, byrow = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("b", 1:6))
  )
  z_reference <- t2a_rates(reference_values)
  state_reference <- t2a_binary_state(
    branch_ids = paste0("b", 1:6),
    values = c("reference", "focal", "reference", "focal", "reference", "focal"),
    availability = rep("available", 6), reason = rep(NA_character_, 6)
  )
  actual_reference <- lt_associate_univariate(
    z_reference$rates$ADD_LT, state_reference, mismatch = "error"
  )$descriptive_scores$point_biserial_r
  expected_reference <- vapply(seq_len(3), function(i) {
    cor(c(0, 1, 0, 1, 0, 1), z_reference$rates$ADD_LT$values[i, ])
  }, numeric(1L))
  expect_equal(actual_reference, expected_reference, tolerance = 1e-14)
})

test_that("T2A-014 rate values are unchanged by association", {
  z <- t2a_rates(t2a_exact_values()); before <- serialize(z$rates$ADD_LT$values, NULL)
  lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(), mismatch = "error")
  expect_identical(serialize(z$rates$ADD_LT$values, NULL), before)
})

test_that("T2A-015 branch-state values are unchanged by association", {
  z <- t2a_rates(t2a_exact_values()); state <- t2a_continuous_state()
  before <- serialize(state$values, NULL)
  lt_associate_univariate(z$rates$ADD_LT, state, mismatch = "error")
  expect_identical(serialize(state$values, NULL), before)
})

test_that("T2A-016 Layers 1 through 4 remain unchanged", {
  z <- t2a_rates(t2a_exact_values()); rate <- z$rates$ADD_LT
  before <- serialize(rate[c("values", "coord_state", "representation_availability",
                             "representation_reason", "measurement_fit_eligibility")], NULL)
  lt_associate_univariate(rate, t2a_continuous_state(), mismatch = "error")
  expect_identical(serialize(rate[c("values", "coord_state",
    "representation_availability", "representation_reason",
    "measurement_fit_eligibility")], NULL), before)
})

test_that("T2A-017 Layer-5 reason ledger is exact", {
  mask <- matrix(TRUE, 2, 4, dimnames = dimnames(t2a_exact_values()))
  mask[1, 4] <- FALSE
  z <- t2a_rates(t2a_exact_values(), mask)
  state <- t2a_continuous_state(
    values = c(0, NA, 2, 3),
    availability = c("available", "unavailable", "available", "available"),
    reason = c(NA, "not_declared", NA, NA))
  d <- lt_association_domain(z$rates$ADD_LT, state,
                             branch_ids = c("b1", "b2", "b4"),
                             mismatch = "error")
  labels <- matrix(d$reason_levels[d$reason], nrow(d$reason), ncol(d$reason),
                   dimnames = dimnames(d$reason))
  expect_identical(unname(labels[1, ]),
                   c("eligible", "trait_unavailable", "branch_user_excluded",
                     "rate_unavailable"))
})

test_that("T2A-018 native association domain is exact", {
  z <- t2a_rates(t2a_exact_values())
  d <- lt_association_domain(z$rates$ADD_LT, t2a_continuous_state(),
                             mismatch = "error")
  expect_identical(d$domain_mode, "NATIVE_ASSOCIATION_DOMAIN")
  expect_identical(d$eligible, z$rates$ADD_LT$representation_availability)
})

test_that("T2A-019 three-rate common domain is exact keyed intersection", {
  values <- t2a_exact_values(); values[1, 1] <- 0
  z <- t2a_rates(values)
  d <- lt_common_association_domain(z$rates, t2a_continuous_state(),
                                    mismatch = "error")
  expected <- Reduce(`&`, lapply(z$rates, `[[`, "representation_availability"))
  expect_identical(d$eligible, expected)
  expect_identical(d$domain_mode, "COMMON_ASSOCIATION_DOMAIN")
})

test_that("T2A-020 common-domain matching never uses position alone", {
  z <- t2a_rates(t2a_exact_values())
  state <- t2a_continuous_state(
    branch_ids = c("b4", "b2", "b1", "b3"), values = c(3, 1, 0, 2))
  a <- lt_associate_univariate(z$rates$ADD_LT, state, mismatch = "error")
  expect_equal(a$estimates$beta, c(2, -2), tolerance = 0)
})

test_that("T2A-021 small branch-ID mismatch warns and uses exact intersection", {
  z <- t2a_rates(t2a_exact_values())
  state <- t2a_continuous_state(
    branch_ids = c("b1", "b2", "b3", "state_only"), values = 0:3)
  expect_warning(d <- lt_association_domain(z$rates$ADD_LT, state),
                 "1 rate-only and 1 state-only")
  expect_identical(d$matched_branch_count, 3L)
  expect_true(all(!d$eligible[, "b4"]))
})

test_that("T2A-022 strict mismatch mode errors", {
  z <- t2a_rates(t2a_exact_values())
  state <- t2a_continuous_state(
    branch_ids = c("b1", "b2", "b3", "state_only"), values = 0:3)
  expect_error(lt_association_domain(z$rates$ADD_LT, state, mismatch = "error"),
               "Exact branch-key mismatch")
})

test_that("T2A-023 ambiguous or duplicate mapping fails", {
  expect_error(t2a_binary_state(branch_ids = c("b1", "b2", "b2", "b4")),
               "duplicate")
})

test_that("T2A-024 no trait variation is explicit non-estimable status", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT,
    t2a_continuous_state(values = rep(1, 4)), mismatch = "error")
  expect_true(all(a$status == "NO_TRAIT_VARIATION"))
})

test_that("T2A-025 binary focal absent has explicit status", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT,
    t2a_binary_state(values = rep("reference", 4)), mismatch = "error")
  expect_true(all(a$status == "BINARY_FOCAL_ABSENT"))
})

test_that("T2A-026 binary reference absent has explicit status", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT,
    t2a_binary_state(values = rep("focal", 4)), mismatch = "error")
  expect_true(all(a$status == "BINARY_REFERENCE_ABSENT"))
})

test_that("T2A-027 all genes remain represented when non-estimable", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT,
    t2a_continuous_state(values = rep(1, 4)), mismatch = "error")
  expect_identical(a$gene_ids, rownames(t2a_exact_values()))
  expect_equal(nrow(a$estimates), 2L)
})

test_that("T2A-028 association performs no automatic rate scaling", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_equal(a$estimates$beta[1], cov(0:3, z$rates$ADD_LT$values[1, ]) /
                 var(0:3), tolerance = 0)
  expect_identical(a$provenance$rate_transformation, "none")
})

test_that("T2A-029 association performs no automatic branch weighting", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_true(a$provenance$unweighted)
  expect_false("weights" %in% names(a))
})

test_that("T2A-030 no trimming winsorization or pseudocount is applied", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  text <- paste(unlist(a$provenance), collapse = " ")
  expect_false(grepl("trim|winsor|pseudocount", text, ignore.case = TRUE))
})

test_that("T2A-031 no ordinary p-value is exposed as calibrated inference", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  released <- unlist(lapply(a[c("estimates", "descriptive_scores",
                                "domain_counts")], names))
  expect_false(any(grepl("p.value|p_value|fdr|q_value|confidence", released,
                         ignore.case = TRUE)))
})

test_that("T2A-032 inference status is point estimate only", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_identical(a$inference_status, "POINT_ESTIMATE_ONLY")
  expect_identical(a$significance_calibration, "NOT_PERFORMED")
})

test_that("T2A-033 association cannot alter C3 or rate hashes", {
  z <- t2a_rates(t2a_exact_values()); before <- z$rates$ADD_LT$output_hashes
  lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(), mismatch = "error")
  expect_identical(z$rates$ADD_LT$output_hashes, before)
  expect_invisible(validate_lt_rate(z$rates$ADD_LT))
})

test_that("T2A-034 trait cannot feed backward into rate construction", {
  expect_false(any(c("trait", "branch_state") %in% names(formals(lt_rate))))
  expect_true(isTRUE(lt_rate_spec("x", "ADD_LT")$trait_independent))
})

test_that("T2A-035 association is not prediction", {
  z <- t2a_rates(t2a_exact_values())
  a <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  expect_true(a$provenance$association_not_prediction)
  expect_false("prediction" %in% names(a))
})

test_that("T2A-036 categorical and ordinal association fail clearly", {
  z <- t2a_rates(t2a_exact_values())
  for (type in c("categorical", "ordinal")) {
    state <- lt_branch_state("x", "t", paste0("b", 1:4), letters[1:4], type,
                             rep("available", 4), rep(NA_character_, 4))
    expect_error(lt_associate_univariate(z$rates$ADD_LT, state),
                 "supports only binary or continuous")
  }
})

test_that("T2A-037 association-domain serialization is exact", {
  z <- t2a_rates(t2a_exact_values())
  x <- lt_association_domain(z$rates$ADD_LT, t2a_continuous_state(),
                             mismatch = "error")
  p <- tempfile(fileext = ".rds"); lt_write_object(x, p)
  expect_identical(lt_read_object(p), x)
})

test_that("T2A-038 branch-state serialization is exact", {
  x <- t2a_continuous_state(); p <- tempfile(fileext = ".rds")
  lt_write_object(x, p)
  expect_identical(lt_read_object(p), x)
})

test_that("T2A-039 association-result serialization is exact", {
  z <- t2a_rates(t2a_exact_values())
  x <- lt_associate_univariate(z$rates$ADD_LT, t2a_continuous_state(),
                               mismatch = "error")
  p <- tempfile(fileext = ".rds"); lt_write_object(x, p)
  expect_identical(lt_read_object(p), x)
})

test_that("T2A-040 ordinary runtime requires no frozen source-file SHA", {
  x <- lt_branch_state(
    "ordinary", "trait", c("b1", "b2"), c(0, 1), "continuous",
    c("available", "available"), rep(NA_character_, 2),
    provenance = list(source_type = "user_supplied", user_declared = TRUE))
  expect_invisible(validate_lt_branch_state(x))
  expect_false(any(grepl("sha", names(x$provenance), ignore.case = TRUE)))
})

test_that("T2A-041 RER dependency is absent", {
  imports <- read.dcf(system.file("DESCRIPTION", package = "LaTerra"),
                      fields = "Imports")
  expect_false(grepl("RER", imports, ignore.case = TRUE))
})

test_that("T2A-042 ASR dependency is absent", {
  expect_length(grep("asr|ancestral|stochastic", getNamespaceExports("LaTerra"),
                     ignore.case = TRUE, value = TRUE), 0L)
})

test_that("T2A-043 trait significance or FDR module is absent", {
  exports <- getNamespaceExports("LaTerra")
  expect_length(grep("pvalue|p_value|fdr|qvalue|significance", exports,
                     ignore.case = TRUE, value = TRUE), 0L)
})
