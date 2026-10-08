t2b_matrix <- function(values, mask = matrix(TRUE, nrow(values), ncol(values))) {
  values <- as.matrix(values); storage.mode(values) <- "double"
  values[!mask] <- NA_real_
  if (is.null(rownames(values))) rownames(values) <- paste0("g", seq_len(nrow(values)))
  if (is.null(colnames(values))) colnames(values) <- paste0("b", seq_len(ncol(values)))
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  state[!mask] <- "NA_struct"
  lt_matrix(values, state, fixture_payload(),
            coordinate_provenance = fixture_coordinate_provenance())
}

t2b_rates <- function(values, mask = matrix(TRUE, nrow(values), ncol(values))) {
  x <- t2b_matrix(values, mask)
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "t2b_fixture")
  rates <- lapply(c("ADD_LT", "GBI", "logGBI_strict"), function(method) {
    lt_rate(x, fit, lt_rate_spec(paste0("t2b_", method), method), eligibility)
  })
  names(rates) <- c("ADD_LT", "GBI", "logGBI")
  list(x = x, eligibility = eligibility, fit = fit, rates = rates)
}

skip_superseded_finite_logGBI <- function() {
  testthat::skip(paste(
    "Superseded 0.0.0.9013 default finite-logGBI assertion;",
    "BOUNDARYIMPL001 tests the boundary object and explicit legacy replay."
  ))
}

t2b_continuous_state <- function(branch_ids = paste0("b", 1:4),
                                 values = 0:3) {
  lt_branch_state(
    "continuous_fixture", "continuous_trait", branch_ids, values,
    "continuous", rep("available", length(branch_ids)),
    rep(NA_character_, length(branch_ids)),
    provenance = list(source_type = "deterministic_fixture")
  )
}

t2b_binary_state <- function(branch_ids = paste0("b", 1:4),
                             values = c("reference", "reference", "focal", "focal")) {
  lt_branch_state(
    "binary_fixture", "binary_trait", branch_ids, values, "binary",
    rep("available", length(branch_ids)),
    rep(NA_character_, length(branch_ids)),
    coding = list(reference_level = "reference", focal_level = "focal"),
    provenance = list(source_type = "deterministic_fixture")
  )
}

t2b_exact_values <- function() {
  matrix(
    c(1, 3, 5, 7, 7, 5, 3, 1), 2, 4, byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("b", 1:4))
  )
}

t2b_binary_null <- function(branch_ids = paste0("b", 1:4)) {
  canonical <- matrix(
    c(
      "reference", "focal", "reference", "focal",
      "focal", "reference", "reference", "focal",
      "focal", "reference", "focal", "reference"
    ),
    nrow = 4, ncol = 3,
    dimnames = list(paste0("b", 1:4), paste0("r", 1:3))
  )
  index <- match(branch_ids, rownames(canonical))
  values <- canonical[index, , drop = FALSE]
  rownames(values) <- branch_ids
  lt_branch_state_ensemble(
    "binary_null", "binary_trait", branch_ids, colnames(values), values,
    type = "binary",
    coding = list(reference_level = "reference", focal_level = "focal"),
    availability = matrix(
      "available", nrow(values), ncol(values), dimnames = dimnames(values)
    ),
    null_hypothesis = "deterministic external fixture",
    generator_provenance = list(generator_type = "fixed_fixture")
  )
}

t2b_binary_bundle <- function(values = t2b_exact_values(), mask = NULL) {
  if (is.null(mask)) mask <- matrix(TRUE, nrow(values), ncol(values),
                                    dimnames = dimnames(values))
  z <- t2b_rates(values, mask)
  state <- t2b_binary_state()
  domain <- lt_association_domain(z$rates$ADD_LT, state, mismatch = "error")
  association <- lt_associate_univariate(
    z$rates$ADD_LT, state, domain = domain
  )
  list(z = z, state = state, domain = domain, association = association,
       ensemble = t2b_binary_null())
}

t2b_calibrate <- function(bundle, alternative = "greater",
                          statistic = "beta") {
  lt_calibrate_association(
    bundle$z$rates$ADD_LT, bundle$state, bundle$association,
    bundle$domain, bundle$ensemble,
    lt_calibration_spec(statistic, alternative,
                        domain_mode = bundle$domain$domain_mode)
  )
}

t2b_common_bundle <- function() {
  values <- rbind(
    g0 = rep(0, 4),
    g1 = c(1, 3, 5, 7),
    g2 = c(7, 5, 3, 1)
  )
  colnames(values) <- paste0("b", 1:4)
  z <- t2b_rates(values)
  state <- t2b_binary_state()
  domain <- lt_common_association_domain(z$rates, state, mismatch = "error")
  associations <- lapply(z$rates, function(rate) {
    lt_associate_univariate(rate, state, domain = domain)
  })
  list(z = z, state = state, domain = domain,
       associations = associations, ensemble = t2b_binary_null())
}

test_that("T2B-001 null ensemble constructor and validator are exact", {
  x <- t2b_binary_null()
  expect_s3_class(x, "lt_branch_state_ensemble")
  expect_identical(dim(x$values), c(4L, 3L))
  expect_false(x$observed_realization_included)
  expect_invisible(validate_lt_branch_state_ensemble(x))
  expect_error(lt_branch_state_ensemble(
    "bad_observed", "binary_trait", x$branch_ids, x$replicate_ids,
    x$values, "binary", x$coding, x$availability, "invalid fixture",
    observed_realization_included = TRUE
  ), "must not be inserted")
  expect_error(lt_branch_state_ensemble(
    "bad_availability_length", "binary_trait", x$branch_ids,
    x$replicate_ids, x$values, "binary", x$coding,
    availability = c("available", "available"),
    null_hypothesis = "invalid availability fixture"
  ), "exactly one entry per branch")
})

test_that("T2B-002 observed and null types must match", {
  b <- t2b_binary_bundle()
  bad <- lt_branch_state_ensemble(
    "continuous_null", "binary_trait", paste0("b", 1:4), "r1",
    matrix(1:4, 4, 1, dimnames = list(paste0("b", 1:4), "r1")),
    "continuous", availability = rep("available", 4),
    null_hypothesis = "type mismatch fixture"
  )
  expect_error(lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, bad
  ), "types must match")
})

test_that("T2B-003 binary coding semantics are exact", {
  b <- t2b_binary_bundle()
  bad <- lt_branch_state_ensemble(
    "bad_coding", "binary_trait", paste0("b", 1:4), "r1",
    matrix(c("reference", "focal", "reference", "focal"), 4, 1,
           dimnames = list(paste0("b", 1:4), "r1")),
    "binary", coding = list(reference_level = "focal",
                              focal_level = "reference"),
    availability = rep("available", 4), null_hypothesis = "coding fixture"
  )
  expect_error(lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, bad
  ), "coding must match exactly")
})

test_that("T2B-004 null branches reorder by exact keys", {
  b <- t2b_binary_bundle()
  reordered <- t2b_binary_null(c("b4", "b2", "b1", "b3"))
  a <- lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, b$ensemble,
    lt_calibration_spec("beta", "greater")
  )
  z <- lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, reordered,
    lt_calibration_spec("beta", "greater")
  )
  expect_identical(z$p_empirical, a$p_empirical)
  expect_identical(z$exceedance_count, a$exceedance_count)
  short <- t2b_binary_null(c("b1", "b2", "b3"))
  expect_error(lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, short
  ), "same exact scientific-key set")
})

test_that("T2B-005 duplicate null branch IDs fail", {
  expect_error(t2b_binary_null(c("b1", "b1", "b3", "b4")), "duplicate")
})

test_that("T2B-006 replicate missingness cannot change a frozen gene domain", {
  values <- t2b_exact_values()
  mask <- matrix(TRUE, 2, 4, dimnames = dimnames(values)); mask[1, 4] <- FALSE
  b <- t2b_binary_bundle(values, mask)
  v <- b$ensemble$values; a <- b$ensemble$availability
  v["b4", "r1"] <- NA; a["b4", "r1"] <- "unavailable"
  e <- lt_branch_state_ensemble(
    "missing_one_branch", "binary_trait", rownames(v), colnames(v), v,
    "binary", list(reference_level = "reference", focal_level = "focal"),
    a, "missing-support fixture"
  )
  out <- lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, e,
    lt_calibration_spec("beta", "greater")
  )
  expect_identical(out$calibration_status,
                   c("CALIBRATED", "NULL_DOMAIN_MISMATCH"))
  expect_true(is.na(out$p_empirical[2]))
})

test_that("T2B-007 one coherent null world is applied to all genes", {
  out <- t2b_calibrate(t2b_binary_bundle(), "greater")
  expect_identical(out$p_empirical, c(0.25, 1))
  expect_true(out$provenance$one_null_world_applied_to_all_genes)
  expect_identical(nrow(out$null_replicate_status), 3L)
})

test_that("T2B-008 greater empirical p is exact", {
  out <- t2b_calibrate(t2b_binary_bundle(), "greater")
  expect_identical(out$exceedance_count, c(0L, 3L))
  expect_identical(out$p_empirical, c(0.25, 1))
})

test_that("T2B-009 less empirical p is exact", {
  out <- t2b_calibrate(t2b_binary_bundle(), "less")
  expect_identical(out$exceedance_count, c(3L, 0L))
  expect_identical(out$p_empirical, c(1, 0.25))
})

test_that("T2B-010 two-sided empirical p is exact", {
  out <- t2b_calibrate(t2b_binary_bundle(), "two_sided")
  expect_identical(out$exceedance_count, c(0L, 0L))
  expect_identical(out$p_empirical, c(0.25, 0.25))
})

test_that("T2B-011 empirical tails include exact ties", {
  b <- t2b_binary_bundle()
  v <- cbind(tie = c("reference", "reference", "focal", "focal"),
             b$ensemble$values)
  e <- lt_branch_state_ensemble(
    "tie_null", "binary_trait", rownames(v), colnames(v), v, "binary",
    list(reference_level = "reference", focal_level = "focal"),
    matrix("available", 4, 4, dimnames = dimnames(v)), "tie fixture"
  )
  out <- lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, e,
    lt_calibration_spec("beta", "greater")
  )
  expect_identical(out$exceedance_count, c(1L, 4L))
  expect_identical(out$p_empirical, c(0.4, 1))
})

test_that("T2B-012 empirical p is never zero", {
  out <- t2b_calibrate(t2b_binary_bundle(), "two_sided")
  expect_true(all(out$p_empirical > 0))
  expect_identical(min(out$p_empirical), 1 / 4)
})

test_that("T2B-013 p resolution is exact", {
  out <- t2b_calibrate(t2b_binary_bundle())
  expect_identical(out$p_resolution, rep(1 / 4, 2))
})

test_that("T2B-014 Monte-Carlo precision is exact", {
  out <- t2b_calibrate(t2b_binary_bundle())
  expect_identical(out$monte_carlo_se,
                   sqrt(out$p_empirical * (1 - out$p_empirical) / 4))
})

test_that("T2B-015 binary non-estimable null makes calibration incomplete", {
  b <- t2b_binary_bundle()
  v <- b$ensemble$values; v[, "r2"] <- "reference"
  e <- lt_branch_state_ensemble(
    "binary_bad_replicate", "binary_trait", rownames(v), colnames(v), v,
    "binary", list(reference_level = "reference", focal_level = "focal"),
    matrix("available", 4, 3, dimnames = dimnames(v)), "one-group fixture"
  )
  out <- lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, e,
    lt_calibration_spec("beta", "greater")
  )
  expect_true(all(out$calibration_status == "NULL_CALIBRATION_INCOMPLETE"))
  expect_true(all(is.na(out$p_empirical)))
  expect_true(any(out$null_replicate_status$status ==
                    "NULL_REPLICATE_NONESTIMABLE"))
})

test_that("T2B-016 continuous non-estimable null makes calibration incomplete", {
  z <- t2b_rates(t2b_exact_values()); state <- t2b_continuous_state()
  d <- lt_association_domain(z$rates$ADD_LT, state, mismatch = "error")
  a <- lt_associate_univariate(z$rates$ADD_LT, state, domain = d)
  v <- cbind(r1 = 3:0, r2 = c(0, 1, 0, 1), r3 = rep(1, 4))
  rownames(v) <- paste0("b", 1:4)
  e <- lt_branch_state_ensemble(
    "continuous_bad_replicate", "continuous_trait", rownames(v),
    colnames(v), v, "continuous",
    availability = matrix("available", 4, 3, dimnames = dimnames(v)),
    null_hypothesis = "no-variation fixture"
  )
  out <- lt_calibrate_association(
    z$rates$ADD_LT, state, a, d, e,
    lt_calibration_spec("beta", "greater")
  )
  expect_true(all(out$calibration_status == "NULL_CALIBRATION_INCOMPLETE"))
  expect_true(all(out$null_nonestimable_count == 1L))
})

test_that("T2B-017 observed non-estimable genes remain represented", {
  z <- t2b_rates(t2b_exact_values())
  state <- t2b_continuous_state(values = rep(1, 4))
  d <- lt_association_domain(z$rates$ADD_LT, state, mismatch = "error")
  a <- lt_associate_univariate(z$rates$ADD_LT, state, domain = d)
  v <- cbind(r1 = 0:3, r2 = 3:0); rownames(v) <- paste0("b", 1:4)
  e <- lt_branch_state_ensemble(
    "continuous_null", "continuous_trait", rownames(v), colnames(v), v,
    "continuous", availability = matrix("available", 4, 2,
                                         dimnames = dimnames(v)),
    null_hypothesis = "observed-nonestimable fixture"
  )
  out <- lt_calibrate_association(z$rates$ADD_LT, state, a, d, e)
  expect_identical(out$gene_ids, z$rates$ADD_LT$ordered_gene_ledger)
  expect_true(all(out$calibration_status ==
                    "OBSERVED_ASSOCIATION_NONESTIMABLE"))
  expect_true(all(is.na(out$p_empirical)))
})

test_that("T2B-018 native calibration identity is exact", {
  b <- t2b_binary_bundle(); out <- t2b_calibrate(b)
  expect_identical(out$domain_mode, "NATIVE_ASSOCIATION_DOMAIN")
  expect_identical(out$domain_identity, b$domain$association_domain_identity)
  expect_identical(out$rate_identity, b$association$rate_identity)
})

test_that("T2B-019 production common association domain is exact", {
  b <- t2b_common_bundle()
  native_add <- lt_association_domain(b$z$rates$ADD_LT, b$state,
                                      mismatch = "error")
  expect_identical(b$domain$eligible,
                   Reduce(`&`, lapply(b$z$rates, `[[`,
                                      "representation_availability")))
  expect_gt(sum(native_add$eligible), sum(b$domain$eligible))
  expect_identical(b$domain$domain_mode, "COMMON_ASSOCIATION_DOMAIN")
})

test_that("T2B-020 the same null worlds apply across representations", {
  b <- t2b_common_bundle()
  out <- Map(function(rate, association) {
    lt_calibrate_association(
      rate, b$state, association, b$domain, b$ensemble,
      lt_calibration_spec("beta", "two_sided",
                          domain_mode = "COMMON_ASSOCIATION_DOMAIN")
    )
  }, b$z$rates, b$associations)
  expect_length(unique(vapply(out, `[[`, character(1),
                              "null_ensemble_identity")), 1L)
  expect_true(all(vapply(out, function(x)
    identical(x$null_replicate_status$replicate_id,
              b$ensemble$replicate_ids), logical(1))))
})

test_that("T2B-021 default logGBI zero cells remain available", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  z <- t2b_rates(values); rate <- z$rates$logGBI
  expect_true(any(rate$zero_encoded %in% TRUE, na.rm = TRUE))
  expect_true(all(rate$representation_availability[rate$zero_encoded %in% TRUE]))
})

test_that("T2B-022 default logGBI zero cells keep certified sentinel values", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  z <- t2b_rates(values); rate <- z$rates$logGBI
  sentinel <- rate$representation_provenance$zero_sentinel
  before <- rate$values[rate$zero_encoded %in% TRUE]
  state <- t2b_binary_state(); d <- lt_association_domain(rate, state,
                                                          mismatch = "error")
  a <- lt_associate_univariate(rate, state, domain = d)
  lt_calibrate_association(rate, state, a, d, t2b_binary_null())
  expect_identical(before, rate$values[rate$zero_encoded %in% TRUE])
  expect_true(all(before == sentinel))
})

test_that("T2B-023 zero-encoded cells never become numeric-zero placeholders", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  rate <- t2b_rates(values)$rates$logGBI
  expect_true(all(rate$values[rate$zero_encoded %in% TRUE] != 0))
})

test_that("T2B-024 default logGBI zero mask is unchanged", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  z <- t2b_rates(values); rate <- z$rates$logGBI
  before <- serialize(rate$zero_encoded, NULL)
  state <- t2b_binary_state(); d <- lt_association_domain(rate, state,
                                                          mismatch = "error")
  a <- lt_associate_univariate(rate, state, domain = d)
  lt_calibrate_association(rate, state, a, d, t2b_binary_null())
  expect_identical(serialize(rate$zero_encoded, NULL), before)
})

test_that("T2B-025 logGBI sentinel is never recomputed from null traits", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  z <- t2b_rates(values); rate <- z$rates$logGBI
  state <- t2b_binary_state(); d <- lt_association_domain(rate, state,
                                                          mismatch = "error")
  a <- lt_associate_univariate(rate, state, domain = d)
  out <- lt_calibrate_association(rate, state, a, d, t2b_binary_null())
  expect_false(out$provenance$zero_sentinel_recomputed)
  expect_true(out$provenance$zero_encoded_cells_consumed_as_is)
  expect_false(any(grepl("m_ref|sentinel", names(out), ignore.case = TRUE)))
})

test_that("T2B-026 logGBI_strict remains distinguishable", {
  skip_superseded_finite_logGBI()
  values <- t2b_exact_values(); values[1, 1] <- 0
  z <- t2b_rates(values)
  strict <- lt_rate(
    z$x, z$fit, lt_rate_spec("strict_fixture", "logGBI_strict"), z$eligibility
  )
  expect_identical(strict$representation, "logGBI_strict")
  expect_false(identical(strict$ordered_cell_identity$available_cell_sha256,
                         z$rates$logGBI$ordered_cell_identity$available_cell_sha256))
})

test_that("T2B-027 calibration cannot change C3", {
  b <- t2b_binary_bundle(); before <- serialize(b$z$fit, NULL)
  t2b_calibrate(b)
  expect_identical(serialize(b$z$fit, NULL), before)
})

test_that("T2B-028 calibration cannot change ADD", {
  b <- t2b_binary_bundle(); before <- b$z$rates$ADD_LT$output_hashes
  t2b_calibrate(b)
  expect_identical(b$z$rates$ADD_LT$output_hashes, before)
})

test_that("T2B-029 calibration cannot change GBI", {
  b <- t2b_common_bundle(); rate <- b$z$rates$GBI
  before <- rate$output_hashes
  lt_calibrate_association(
    rate, b$state, b$associations$GBI, b$domain, b$ensemble,
    lt_calibration_spec("beta", "two_sided",
                        domain_mode = "COMMON_ASSOCIATION_DOMAIN")
  )
  expect_identical(rate$output_hashes, before)
})

test_that("T2B-030 calibration cannot change default logGBI", {
  skip_superseded_finite_logGBI()
  b <- t2b_common_bundle(); rate <- b$z$rates$logGBI
  before <- rate$output_hashes
  lt_calibrate_association(
    rate, b$state, b$associations$logGBI, b$domain, b$ensemble,
    lt_calibration_spec("beta", "two_sided",
                        domain_mode = "COMMON_ASSOCIATION_DOMAIN")
  )
  expect_identical(rate$output_hashes, before)
})

test_that("T2B-031 calibration performs no null generation", {
  out <- t2b_calibrate(t2b_binary_bundle())
  expect_identical(out$provenance$null_generation, "not_performed")
  expect_false(out$provenance$hidden_rng)
})

test_that("T2B-032 no naive branch-permutation default exists", {
  code <- paste(deparse(body(lt_calibrate_association)), collapse = " ")
  expect_false(grepl("sample\\(|shuffle|permute", code, ignore.case = TRUE))
})

test_that("T2B-033 BH adjustment is exact", {
  x <- t2b_calibrate(t2b_binary_bundle())
  out <- lt_adjust_calibration(x, "BH")
  expect_identical(out$adjusted_p,
                   p.adjust(x$p_empirical, method = "BH"))
})

test_that("T2B-034 BY adjustment is exact", {
  x <- t2b_calibrate(t2b_binary_bundle())
  out <- lt_adjust_calibration(x, "BY")
  expect_identical(out$adjusted_p,
                   p.adjust(x$p_empirical, method = "BY"))
})

test_that("T2B-035 Holm adjustment is exact and is FWER", {
  x <- t2b_calibrate(t2b_binary_bundle())
  out <- lt_adjust_calibration(x, "Holm")
  expect_identical(out$adjusted_p,
                   p.adjust(x$p_empirical, method = "holm"))
  expect_identical(out$procedure_class, "FWER_PROCEDURE")
})

test_that("T2B-036 Bonferroni adjustment is exact and is FWER", {
  x <- t2b_calibrate(t2b_binary_bundle())
  out <- lt_adjust_calibration(x, "Bonferroni")
  expect_identical(out$adjusted_p,
                   p.adjust(x$p_empirical, method = "bonferroni"))
  expect_identical(out$procedure_class, "FWER_PROCEDURE")
})

test_that("T2B-037 multiple-testing family boundaries are explicit", {
  x <- t2b_calibrate(t2b_binary_bundle())
  out <- lt_adjust_calibration(x, "none")
  expect_identical(out$n_tests, 2L)
  expect_match(out$provenance$family_definition,
               "one trait x one rate representation")
  expect_false(out$provenance$automatic)
})

test_that("T2B-038 deterministic replay is exact", {
  b <- t2b_binary_bundle()
  a <- t2b_calibrate(b); z <- t2b_calibrate(b)
  expect_identical(a, z)
  expect_identical(a$calibration_identity, z$calibration_identity)
})

test_that("T2B-039 ensemble calibration and adjustment serialize exactly", {
  b <- t2b_binary_bundle(); c <- t2b_calibrate(b)
  objects <- list(b$ensemble, lt_calibration_spec(), c,
                  lt_adjust_calibration(c, "BY"))
  for (x in objects) {
    p <- tempfile(fileext = ".rds")
    lt_write_object(x, p)
    expect_identical(lt_read_object(p), x)
  }
})

test_that("T2B-040 ordinary runtime needs no frozen raw-file SHA", {
  x <- t2b_binary_null()
  expect_false(any(grepl("sha", names(x$generator_provenance),
                         ignore.case = TRUE)))
  x$metadata$display_note <- "user-editable"
  expect_invisible(validate_lt_branch_state_ensemble(x))
  p <- tempfile(fileext = ".rds"); lt_write_object(x, p)
  expect_error(lt_write_object(x, p), "already exists")
})

test_that("T2B-041 RER remains absent", {
  imports <- read.dcf(system.file("DESCRIPTION", package = "LaTerra"),
                      fields = "Imports")
  expect_false(grepl("RER", imports, ignore.case = TRUE))
})

test_that("T2B-042 prediction remains absent", {
  exports <- getNamespaceExports("LaTerra")
  expect_length(grep("predict|prediction", exports, ignore.case = TRUE,
                     value = TRUE), 0L)
})

test_that("T2B-043 frozen 9013 is superseded only by explicit boundary contract", {
  skip_superseded_finite_logGBI()
  expect_identical(as.character(packageVersion("LaTerra")), "0.0.0.9014")
  values <- t2b_exact_values(); values[1, 1] <- 0
  rate <- t2b_rates(values)$rates$logGBI
  expect_true(any(rate$zero_encoded %in% TRUE, na.rm = TRUE))
  expect_match(rate$representation_provenance$zero_policy,
               "global_positive_log_min_minus_10")
})

test_that("T2B-SPEC-001 calibration specification and result fields are explicit", {
  spec <- lt_calibration_spec("beta", "greater", "none",
                              "NATIVE_ASSOCIATION_DOMAIN")
  expect_s3_class(spec, "lt_calibration_spec")
  expect_error(lt_calibration_spec("hidden_statistic"), "Unknown")
  out <- t2b_calibrate(t2b_binary_bundle())
  required <- c(
    "calibration_id", "association_identity", "rate_identity",
    "branch_state_identity", "null_ensemble_identity", "calibration_spec",
    "statistic", "alternative", "gene_ids", "observed_statistic",
    "null_replicate_count", "exceedance_count", "p_empirical",
    "p_resolution", "monte_carlo_se", "calibration_status", "reason",
    "domain_identity", "provenance", "metadata"
  )
  expect_true(all(required %in% names(out)))
  bad <- out
  bad$null_replicate_status$nonestimable_gene_count <-
    as.double(bad$null_replicate_status$nonestimable_gene_count)
  expect_error(validate_lt_calibration(bad), "status ledger is malformed")
})

test_that("T2B-SPEC-002 statistic vocabularies are type-specific", {
  b <- t2b_binary_bundle()
  expect_error(lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, b$ensemble,
    lt_calibration_spec("pearson_r", "two_sided")
  ), "not defined")
  expect_s3_class(lt_calibrate_association(
    b$z$rates$ADD_LT, b$state, b$association, b$domain, b$ensemble,
    lt_calibration_spec("point_biserial_r", "two_sided")
  ), "lt_calibration")

  z <- t2b_rates(t2b_exact_values()); state <- t2b_continuous_state()
  d <- lt_association_domain(z$rates$ADD_LT, state, mismatch = "error")
  a <- lt_associate_univariate(z$rates$ADD_LT, state, domain = d)
  v <- cbind(r1 = 3:0, r2 = c(0, 1, 3, 2), r3 = c(2, 0, 1, 3))
  rownames(v) <- paste0("b", 1:4)
  e <- lt_branch_state_ensemble(
    "continuous_complete", "continuous_trait", rownames(v), colnames(v), v,
    "continuous", availability = matrix("available", 4, 3,
                                         dimnames = dimnames(v)),
    null_hypothesis = "continuous statistic fixture"
  )
  for (statistic in c("beta", "pearson_r", "spearman_rho")) {
    expect_s3_class(lt_calibrate_association(
      z$rates$ADD_LT, state, a, d, e,
      lt_calibration_spec(statistic, "two_sided")
    ), "lt_calibration")
  }
  expect_error(lt_calibrate_association(
    z$rates$ADD_LT, state, a, d, e,
    lt_calibration_spec("point_biserial_r", "two_sided")
  ), "not defined")
})

test_that("T2B-SPEC-003 continuous null values are consumed without preprocessing", {
  values <- matrix(c(10.5, -2.25, 7, 100), 2, 2,
                   dimnames = list(c("b1", "b2"), c("r1", "r2")))
  x <- lt_branch_state_ensemble(
    "raw_continuous", "continuous_trait", rownames(values), colnames(values),
    values, "continuous",
    availability = matrix("available", 2, 2, dimnames = dimnames(values)),
    null_hypothesis = "raw-value preservation fixture"
  )
  expect_identical(x$values, values)
  expect_false(any(c("center", "scale", "transform", "impute") %in%
                     names(x$generator_provenance)))
})
