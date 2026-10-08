boundary_fixture <- function() {
  values <- matrix(
    c(
      0, 1, 2, 0,
      1, 0, 3, 4,
      2, 3, 0, 1
    ), 3, 4, byrow = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("b", 1:4))
  )
  coord <- matrix("observed", 3, 4, dimnames = dimnames(values))
  x <- lt_matrix(
    values, coord, fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "boundary_fixture")
  gbi <- lt_rate(
    x, fit, lt_rate_spec("boundary_GBI", "GBI"), eligibility
  )
  state <- lt_branch_state(
    "boundary_binary", "boundary_trait", colnames(values),
    c("reference", "reference", "focal", "focal"), "binary",
    rep("available", 4), rep(NA_character_, 4),
    coding = list(reference_level = "reference", focal_level = "focal"),
    provenance = list(source_type = "fixed_fixture")
  )
  null_values <- matrix(
    c(
      "reference", "focal", "reference",
      "focal", "reference", "reference",
      "reference", "reference", "focal",
      "focal", "focal", "focal"
    ), 4, 3, byrow = TRUE,
    dimnames = list(colnames(values), paste0("r", 1:3))
  )
  null <- lt_branch_state_ensemble(
    "boundary_null", "boundary_trait", colnames(values),
    colnames(null_values), null_values, "binary",
    coding = list(reference_level = "reference", focal_level = "focal"),
    availability = matrix(
      "available", 4, 3, dimnames = dimnames(null_values)
    ),
    null_hypothesis = "fixed unit-test null",
    generator_provenance = list(generator_type = "fixed_fixture")
  )
  list(x = x, eligibility = eligibility, fit = fit, gbi = gbi,
       state = state, null = null)
}

rehash_gbi <- function(gbi) {
  gbi$output_hashes <- LaTerra:::.lt_rate_authoritative_hashes(gbi)
  gbi
}

test_that("BOUNDARY-001 primary state/reason truth table is exhaustive", {
  z <- boundary_fixture()
  x <- lt_log_relative(z$gbi)
  expect_identical(sort(unique(as.vector(ratio_state(x)))),
                   c("EXACT_ZERO_BOUNDARY", "POSITIVE"))
  expect_true(all(unavailability_reason(x) == "not_applicable"))

  reasons <- setdiff(LaTerra:::.lt_log_relative_reason_levels(),
                     "not_applicable")
  for (reason in reasons) {
    g <- z$gbi
    g$values[1, 1] <- NA_real_
    g$representation_availability[1, 1] <- FALSE
    code <- match(reason, g$representation_reason_levels)
    g$representation_reason[1, 1] <- code
    g$value_reason[1, 1] <- code
    g$zero_origin_state[1, 1] <- as.raw(0L)
    g$ordered_cell_identity$available_cell_sha256 <-
      LaTerra:::.lt_rate_keyed_mask_hash(
        g$representation_availability, g$ordered_gene_ledger,
        g$ordered_branch_ledger
      )
    g <- rehash_gbi(g)
    y <- lt_log_relative(g)
    expect_identical(ratio_state(y)[1, 1], "UNAVAILABLE")
    expect_identical(unavailability_reason(y)[1, 1], reason)
  }
})

test_that("BOUNDARY-002 exact_zero is mandatory three-valued", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  g <- z$gbi
  positive_cell <- which(g$values > 0, arr.ind = TRUE)[1, ]
  zero_cell <- which(g$values == 0, arr.ind = TRUE)[1, ]
  g$values[positive_cell[1], positive_cell[2]] <- NA_real_
  g$representation_availability[positive_cell[1], positive_cell[2]] <- FALSE
  code <- match("numeric_indeterminate", g$representation_reason_levels)
  g$representation_reason[positive_cell[1], positive_cell[2]] <- code
  g$value_reason[positive_cell[1], positive_cell[2]] <- code
  g$zero_origin_state[positive_cell[1], positive_cell[2]] <- as.raw(0L)
  g$ordered_cell_identity$available_cell_sha256 <-
    LaTerra:::.lt_rate_keyed_mask_hash(
      g$representation_availability, g$ordered_gene_ledger,
      g$ordered_branch_ledger
    )
  g <- rehash_gbi(g); x <- lt_log_relative(g)
  ez <- exact_zero(x)
  expect_false(ez[which(ratio_state(x) == "POSITIVE", arr.ind = TRUE)[1, 1],
                  which(ratio_state(x) == "POSITIVE", arr.ind = TRUE)[1, 2]])
  expect_true(ez[zero_cell[1], zero_cell[2]])
  expect_true(is.na(ez[positive_cell[1], positive_cell[2]]))
})

test_that("BOUNDARY-003 finite positive logs occur only on POSITIVE", {
  x <- lt_log_relative(boundary_fixture()$gbi)
  state <- ratio_state(x)
  expect_true(all(is.finite(positive_log_values(x)[state == "POSITIVE"])))
  expect_true(all(is.na(positive_log_values(x)[state != "POSITIVE"])))
  expect_identical(ratio_availability(x), state != "UNAVAILABLE")
  expect_identical(positive_log_availability(x), state == "POSITIVE")
})

test_that("BOUNDARY-004 exact zero requires certified numerator origin", {
  z <- boundary_fixture(); g <- z$gbi
  cell <- which(g$values == 0, arr.ind = TRUE)[1, ]
  g$zero_origin_state[cell[1], cell[2]] <- as.raw(0L)
  g <- rehash_gbi(g)
  x <- lt_log_relative(g)
  expect_identical(ratio_state(x)[cell[1], cell[2]], "UNAVAILABLE")
  expect_identical(unavailability_reason(x)[cell[1], cell[2]],
                   "numeric_indeterminate")
})

test_that("BOUNDARY-005 floating underflow zero fails closed", {
  z <- boundary_fixture(); g <- z$gbi
  cell <- which(g$values > 0, arr.ind = TRUE)[1, ]
  g$values[cell[1], cell[2]] <- 0
  expect_identical(as.integer(g$zero_origin_state[cell[1], cell[2]]), 1L)
  g <- rehash_gbi(g)
  x <- lt_log_relative(g)
  expect_identical(ratio_state(x)[cell[1], cell[2]], "UNAVAILABLE")
  expect_identical(unavailability_reason(x)[cell[1], cell[2]],
                   "numeric_indeterminate")
})

test_that("BOUNDARY-006 vocabularies and axes are identity protected", {
  x <- lt_log_relative(boundary_fixture()$gbi)
  bad <- x; bad$ratio_state_levels[1] <- "OTHER"
  expect_error(validate_lt_log_relative(bad), "vocabulary")
  bad <- x; bad$ordered_branch_ledger <- rev(bad$ordered_branch_ledger)
  expect_error(validate_lt_log_relative(bad), "aligned|axis")
})

test_that("BOUNDARY-007 semantic and audit identities are separate", {
  x <- lt_log_relative(boundary_fixture()$gbi)
  labels <- x; labels$gene_display_labels[] <- "renamed"
  expect_invisible(validate_lt_log_relative(labels))
  expect_identical(labels$semantic_identity, x$semantic_identity)
  expect_identical(labels$audit_identity, x$audit_identity)

  migrated <- x
  migrated$migration_provenance <- list(note = "explicit migration evidence")
  migrated$audit_identity <- LaTerra:::.lt_log_relative_identity_hash(
    LaTerra:::.lt_log_relative_audit_payload(migrated),
    migrated$canonicalization_version
  )
  expect_invisible(validate_lt_log_relative(migrated))
  expect_identical(migrated$semantic_identity, x$semantic_identity)
  expect_false(identical(migrated$audit_identity, x$audit_identity))
})

test_that("BOUNDARY-008 stale dependencies reject at the boundary", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  expect_identical(lt_log_relative_dependency_status(x, z$gbi)$status,
                   "CURRENT")
  changed <- z$gbi
  changed$values[changed$values > 0][1] <-
    changed$values[changed$values > 0][1] + 1e-6
  changed <- rehash_gbi(changed)
  expect_identical(lt_log_relative_dependency_status(x, changed)$status,
                   "STALE_DEPENDENCY")
  expect_error(lt_log_relative_dependency_status(x, changed, error = TRUE),
               class = "lt_error_stale_log_relative_dependency")
})

test_that("BOUNDARY-009 old constructor and config receive classed guidance", {
  expect_error(lt_rate_spec("old", "logGBI"),
               class = "lt_error_logGBI_constructor_moved")
  config <- yaml::read_yaml(system.file(
    "recipes", "marine_2026_legacy.yaml", package = "LaTerra"
  ))
  config$rate$method <- "logGBI"
  expect_error(validate_lt_config(config),
               class = "lt_error_logGBI_constructor_moved")
})

test_that("BOUNDARY-010 internal bypass and numeric coercions fail closed", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  expect_error(
    LaTerra:::.lt_boundary_score(
      x, z$state$values, rep(TRUE, 4), z$state
    ), class = "lt_error_log_relative_internal_bypass"
  )
  expect_error(as.numeric(x), class = "lt_error_log_relative_numeric_coercion")
  expect_error(as.matrix(x), class = "lt_error_log_relative_numeric_coercion")
  expect_error(mean(x), class = "lt_error_log_relative_numeric_coercion")
  expect_error(rowMeans(x))
  expect_error(stats::t.test(x), class = "lt_error_log_relative_numeric_coercion")
  expect_error(lt_associate_univariate(x, z$state),
               class = "lt_error_log_relative_requires_component_association")
  expect_identical(positive_log_values(x), x$positive_log_values)
})

test_that("BOUNDARY-011 source-free and ambiguous migration refuse", {
  z <- boundary_fixture()
  expect_error(lt_migrate_logGBI(z$gbi),
               class = "lt_error_logGBI_migration_source_required")
  expect_error(lt_migrate_logGBI(z$gbi, z$gbi),
               class = "lt_error_unclassified_historical_logGBI")
  old <- lt_rate(
    z$x, z$fit, lt_rate_spec("boundary_old_strict", "logGBI_strict"),
    z$eligibility
  )
  positive <- which(z$gbi$values > 0)[1]
  old$values[positive] <- min(log(z$gbi$values[z$gbi$values > 0])) - 10
  old$output_hashes <- LaTerra:::.lt_rate_authoritative_hashes(old)
  expect_error(
    lt_migrate_logGBI(old, z$gbi),
    class = "lt_error_log_relative_source_mismatch"
  )
})

test_that("BOUNDARY-012 source m_ref is preserved without log-exp recovery", {
  x <- lt_log_relative(boundary_fixture()$gbi)
  encoded <- logGBI_encoded(x)
  expect_identical(encoded$m_ref,
                   x$representation_provenance$original_source_m_ref)
  expect_identical(encoded$sentinel, log(encoded$m_ref) - 10)
  expect_false(encoded$default_continuous_inference_authorized)
})

test_that("BOUNDARY-013 all-zero boundary is valid and legacy encoding is not", {
  z <- boundary_fixture(); g <- z$gbi
  available <- g$representation_availability
  g$values[available] <- 0
  g$zero_origin_state[available] <- as.raw(2L)
  g <- rehash_gbi(g)
  x <- lt_log_relative(g)
  expect_invisible(validate_lt_log_relative(x))
  expect_false(any(positive_log_availability(x)))
  expect_error(logGBI_encoded(x),
               class = "lt_error_logGBI_encoded_nonestimable")
})

test_that("BOUNDARY-014 association keeps components separate and no joint p", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  a <- lt_associate_log_relative(x, z$state, mismatch = "error")
  expect_s3_class(a, "lt_boundary_association")
  expect_identical(names(a$components), c("ZERO_MASS", "POSITIVE_LOG"))
  expect_false("joint_p" %in% names(a))
  expect_identical(a$joint_inference_status, "NOT_DEFINED_BY_CONTRACT")
  expect_identical(component(a, "ZERO_MASS")$estimand,
                   "delta_z_focal_minus_reference")
})

test_that("BOUNDARY-015 null calibration is plus-one and separately adjusted", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  a <- lt_associate_log_relative(x, z$state, mismatch = "error")
  cal <- lt_calibrate_log_relative(
    x, z$state, a, list(NULL_A = z$null, NULL_B = z$null)
  )
  expect_length(cal$calibrations, 4L)
  expect_false("joint_p" %in% names(cal))
  for (one in cal$calibrations) {
    ok <- one$calibration_status == "CALIBRATED"
    expect_identical(
      one$p_empirical[ok],
      (1 + one$exceedance_count[ok]) /
        (one$null_replicate_count[ok] + 1)
    )
  }
  adj <- lt_adjust_log_relative(cal, c("BH", "BY"))
  expect_length(adj$adjustments, 8L)
  families <- vapply(adj$adjustments, `[[`, character(1L), "family_identity")
  expect_length(unique(families), 8L)
})

test_that("BOUNDARY-016 legacy view materializes diagnostically but refuses inference", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  encoded <- logGBI_encoded(x)
  expect_identical(as.matrix(encoded), encoded$values)
  expect_invisible(validate_lt_logGBI_encoded_source(encoded, x))
  expect_error(lt_associate_univariate(encoded, z$state),
               class = "lt_error_legacy_encoded_inference_not_supported")
  expect_error(
    lt_associate_univariate(
      encoded, z$state, encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ), class = "lt_error_legacy_encoded_inference_not_supported"
  )
})

test_that("BOUNDARY-017 legacy calibration and adjustment entry points refuse", {
  z <- boundary_fixture(); x <- lt_log_relative(z$gbi)
  encoded <- logGBI_encoded(x)
  expect_error(
    lt_calibrate_logGBI_encoded(
      encoded, z$state, NULL, list(NULL_A = z$null)
    ), class = "lt_error_legacy_encoded_inference_not_supported"
  )
  expect_error(
    lt_calibrate_logGBI_encoded(
      encoded, z$state, NULL, list(NULL_A = z$null),
      encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ), class = "lt_error_legacy_encoded_inference_not_supported"
  )
  expect_error(
    lt_adjust_logGBI_encoded(list()),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
})

test_that("BOUNDARY-018 boundary objects round trip without migration", {
  x <- lt_log_relative(boundary_fixture()$gbi)
  path <- tempfile(fileext = ".rds")
  lt_write_object(x, path)
  restored <- lt_read_object(path)
  expect_identical(restored, x)
  expect_identical(restored$semantic_identity, x$semantic_identity)
  expect_identical(restored$audit_identity, x$audit_identity)
})
