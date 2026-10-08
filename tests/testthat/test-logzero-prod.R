lzprod_matrix <- function(values) {
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  if (is.null(rownames(values))) rownames(values) <- paste0("lz_g", seq_len(nrow(values)))
  if (is.null(colnames(values))) colnames(values) <- paste0("lz_b", seq_len(ncol(values)))
  state <- matrix("observed", nrow(values), ncol(values), dimnames = dimnames(values))
  lt_matrix(values, state, fixture_payload(), coordinate_provenance = NULL)
}

lzprod_objects <- function(values = matrix(
  c(1, 0, 2, 0, 3, 4), 2, 3, byrow = TRUE
)) {
  testthat::skip(paste(
    "Superseded 0.0.0.9013 finite-sentinel fixture;",
    "BOUNDARYIMPL001 conformance tests exercise lt_log_relative and explicit replay."
  ))
  x <- lzprod_matrix(values)
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "lzprod_fixture")
  make <- function(method) lt_rate(
    x, fit, lt_rate_spec(paste0("lzprod_", method), method), eligibility
  )
  list(
    x = x, eligibility = eligibility, fit = fit,
    ADD_LT = make("ADD_LT"), GBI = make("GBI"),
    logGBI = make("logGBI"), logGBI_strict = make("logGBI_strict")
  )
}

test_that("LZPROD-001 default logGBI accepts Y zero with positive mu", {
  z <- lzprod_objects()
  edge <- z$fit$edge_linear
  mask <- matrix(FALSE, nrow(z$x$values), ncol(z$x$values), dimnames = dimnames(z$x$values))
  mask[edge[z$fit$observed_value == 0 & z$fit$mu > 0]] <- TRUE
  expect_true(any(mask))
  expect_true(all(z$logGBI$representation_availability[mask]))
  expect_true(all(is.finite(z$logGBI$values[mask])))
})

test_that("LZPROD-002 sentinel is the positive strict minimum minus ten", {
  z <- lzprod_objects()
  expected <- min(z$logGBI_strict$values, na.rm = TRUE) - 10
  expect_identical(z$logGBI$representation_provenance$zero_sentinel, expected)
  expect_identical(unique(z$logGBI$values[z$logGBI$zero_encoded %in% TRUE]), expected)
})

test_that("LZPROD-003 natural logarithm is used on positive cells", {
  z <- lzprod_objects()
  positive <- z$GBI$representation_availability & z$GBI$values > 0
  expect_equal(z$logGBI$values[positive], log(z$GBI$values[positive]), tolerance = 0)
  expect_identical(z$logGBI$representation_provenance$log_base, "natural")
})

test_that("LZPROD-004 one global sentinel is used per object", {
  z <- lzprod_objects()
  encoded <- z$logGBI$zero_encoded %in% TRUE
  expect_gt(sum(encoded), 1L)
  expect_length(unique(z$logGBI$values[encoded]), 1L)
  expect_identical(z$logGBI$representation_provenance$m_ref_source_domain_count,
                   as.integer(sum(!encoded & z$logGBI$representation_availability)))
})

test_that("LZPROD-005 positive values are exact-identical to strict logGBI", {
  z <- lzprod_objects()
  positive <- z$logGBI_strict$representation_availability
  expect_identical(z$logGBI$values[positive], z$logGBI_strict$values[positive])
})

test_that("LZPROD-006 zero encoding preserves source Y equal zero", {
  z <- lzprod_objects(); before <- serialize(z$x$values, NULL)
  encoded <- z$logGBI$zero_encoded %in% TRUE
  expect_true(all(z$x$values[encoded] == 0))
  expect_identical(serialize(z$x$values, NULL), before)
})

test_that("LZPROD-007 zero encoding preserves GBI equal zero", {
  z <- lzprod_objects(); encoded <- z$logGBI$zero_encoded %in% TRUE
  expect_identical(z$GBI$values[encoded], rep(0, sum(encoded)))
})

test_that("LZPROD-008 zero_encoded mask is exact and first class", {
  z <- lzprod_objects()
  mu <- matrix(NA_real_, nrow(z$x$values), ncol(z$x$values), dimnames = dimnames(z$x$values))
  mu[z$fit$edge_linear] <- z$fit$mu
  expected <- z$x$values == 0 & mu > 0
  expected[!z$logGBI$representation_availability] <- NA
  expect_identical(z$logGBI$zero_encoded, expected)
  expect_identical(z$logGBI$output_hashes$zero_encoded_sha256,
                   digest::digest(z$logGBI$zero_encoded, algo = "sha256"))
})

test_that("LZPROD-009 baseline-zero cells remain unavailable", {
  z <- lzprod_objects(matrix(c(1, 2, 0, 0), 2, 2, byrow = TRUE))
  boundary <- matrix(FALSE, 2, 2, dimnames = dimnames(z$x$values))
  boundary[z$fit$edge_linear[z$fit$mu == 0]] <- TRUE
  expect_true(any(boundary))
  expect_true(all(is.na(z$logGBI$values[boundary])))
  expect_true(all(is.na(z$logGBI$zero_encoded[boundary])))
  expect_true(all(lt_rate_reason(z$logGBI)[boundary] == "baseline_zero"))
})

test_that("LZPROD-010 positive Y with zero mu remains fatal", {
  z <- lzprod_objects(matrix(c(1, 2, 3, 4), 2, 2))
  z$fit$mu[1L] <- 0
  expect_error(lt_rate(
    z$x, z$fit, lt_rate_spec("fatal", "logGBI"), z$eligibility
  ), "positive observed payload|fatal_positive_y_zero_mu")
})

test_that("LZPROD-011 zero-state encoding is not a pseudocount", {
  z <- lzprod_objects()
  expect_identical(z$logGBI$representation_provenance$pseudocount, "none")
  expect_false(any(c("pseudocount", "floor") %in% names(z$logGBI$rate_spec$parameters)))
  expect_identical(z$GBI$values[z$GBI$values == 0],
                   rep(0, sum(z$GBI$values == 0, na.rm = TRUE)))
})

test_that("LZPROD-012 production logGBI stores no Inf minus Inf or NaN", {
  z <- lzprod_objects()
  expect_false(any(is.infinite(z$logGBI$values), na.rm = TRUE))
  expect_false(any(is.nan(z$logGBI$values)))
  expect_false(any(is.infinite(z$logGBI_strict$values), na.rm = TRUE))
  expect_false(any(is.nan(z$logGBI_strict$values)))
})

test_that("LZPROD-013 default logGBI and GBI have exact keyed domains", {
  z <- lzprod_objects()
  expect_identical(z$logGBI$representation_availability,
                   z$GBI$representation_availability)
  expect_identical(z$logGBI$ordered_cell_identity$available_cell_sha256,
                   z$GBI$ordered_cell_identity$available_cell_sha256)
})

# LZPROD-014 through LZPROD-016 are frozen marine count gates executed by
# 09_SCRIPTS/freeze_loggbi_prod.R and independent_validation.R.

test_that("LZPROD-017 logGBI_strict preserves positive-only behavior", {
  z <- lzprod_objects()
  zero <- z$x$values == 0
  expect_true(all(is.na(z$logGBI_strict$values[zero])))
  expect_true(all(lt_rate_reason(z$logGBI_strict)[zero] ==
                    "input_zero_log_undefined"))
  positive <- z$logGBI_strict$representation_availability
  expect_identical(z$logGBI_strict$values[positive],
                   log(z$GBI$values[positive]))
})

# LZPROD-018 is the frozen marine strict numeric-hash gate.

test_that("LZPROD-019 old serialized strict object is not reinterpreted", {
  path <- system.file("extdata", "legacy_logGBI_9010.rds", package = "LaTerra")
  expect_true(nzchar(path))
  legacy <- lt_read_object(path)
  expect_identical(legacy$representation, "logGBI")
  expect_null(legacy$zero_encoded)
  expect_identical(
    LaTerra:::.lt_rate_log_contract(legacy),
    "LEGACY_DEVELOPMENT_LOGGBI_STRICT_NOT_REINTERPRETED"
  )
})

# LZPROD-020 is the frozen marine association-domain identity gate.

test_that("LZPROD-021 association result exposes no p-value or FDR", {
  values <- matrix(c(1, 0, 2, 3, 0, 4, 2, 1), 2, 4, byrow = TRUE,
                   dimnames = list(c("g1", "g2"), paste0("b", 1:4)))
  z <- lzprod_objects(values)
  state <- lt_branch_state(
    "lz_binary", "lz_trait", paste0("b", 1:4),
    c("reference", "reference", "focal", "focal"), "binary",
    rep("available", 4), rep(NA_character_, 4),
    coding = list(reference_level = "reference", focal_level = "focal"),
    provenance = list(source_type = "deterministic_fixture")
  )
  association <- lt_associate_univariate(z$logGBI, state, mismatch = "error")
  exposed <- unlist(lapply(
    association[c("estimates", "descriptive_scores", "domain_counts")], names
  ))
  expect_false(any(grepl("p.value|p_value|fdr|q_value", exposed, ignore.case = TRUE)))
})

test_that("LZPROD-022 exact minus-ten depth lies between K 1e4 and 1e6", {
  z <- lzprod_objects()
  min_log <- z$logGBI$representation_provenance$L_min_plus
  exact <- z$logGBI$representation_provenance$zero_sentinel
  expect_lt(exact, min_log - log(1e4))
  expect_lt(min_log - log(1e6), exact)
  expect_gt(exp(10), 1e4)
  expect_lt(exp(10), 1e6)
})

# LZPROD-023 is the frozen marine C3/ADD/GBI regression barrier.

test_that("LZPROD-024 deterministic replay is exact", {
  a <- lzprod_objects()$logGBI
  b <- lzprod_objects()$logGBI
  expect_identical(a, b)
  expect_identical(a$output_hashes, b$output_hashes)
})

test_that("LZPROD-025 serialization preserves zero mask and provenance", {
  object <- lzprod_objects()$logGBI
  path <- tempfile(fileext = ".rds")
  lt_write_object(object, path)
  observed <- lt_read_object(path)
  expect_identical(observed, object)
  expect_identical(observed$zero_encoded, object$zero_encoded)
  expect_identical(observed$representation_provenance,
                   object$representation_provenance)
})

test_that("LZPROD-026 ordinary runtime does not require input-file SHA", {
  z <- lzprod_objects()
  expect_true(is.na(z$x$coordinate_provenance$reference_tree_sha256))
  expect_invisible(validate_lt_rate(z$logGBI))
  expect_false("frozen_file_sha256" %in% names(z$logGBI))
})

test_that("LZPROD-027 prediction metadata requires training freeze", {
  z <- lzprod_objects()
  expect_identical(
    z$logGBI$representation_provenance$prediction_use_status,
    "ASSOCIATION_OK_TRAINING_FREEZE_REQUIRED_FOR_HELDOUT_PREDICTION"
  )
  expect_false("prediction" %in% names(formals(lt_rate)))
})
