fix001_matrix <- function(values = matrix(c(1, 2, 3, 5, 7, 11), 2, 3)) {
  values <- as.matrix(values)
  storage.mode(values) <- "double"
  dimnames(values) <- list(c("gene_key_1", "gene_key_2"),
                           c("branch_key_1", "branch_key_2", "branch_key_3"))
  state <- matrix("observed", nrow(values), ncol(values),
                  dimnames = dimnames(values))
  lt_matrix(values, state, fixture_payload(), coordinate_provenance = NULL)
}

fix001_objects <- function() {
  x <- fix001_matrix()
  eligibility <- lt_measurement_eligibility(x)
  c3 <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "fix001_c3")
  add <- lt_rate(x, c3, lt_rate_spec("fix001_add", "ADD_LT"), eligibility)
  gbi <- lt_rate(x, c3, lt_rate_spec("fix001_gbi", "GBI"), eligibility)
  loggbi <- lt_rate(
    x, c3, lt_rate_spec("fix001_loggbi", "logGBI_strict"), eligibility
  )
  c2 <- lt_c2_fit(lt_c2_domain(x, eligibility), run_id = "fix001_c2")
  c2_rate <- lt_c2_sensitivity(x, c2, "GBI_C2", eligibility)
  common <- lt_rate_common_domain(add, gbi, loggbi)
  list(
    x = x, eligibility = eligibility, c3 = c3, add = add, gbi = gbi,
    loggbi = loggbi, c2 = c2, c2_rate = c2_rate, common = common
  )
}

test_that("FIX001-A embedded Layer-2 status mutation fails closed", {
  z <- fix001_objects(); bad <- z$add
  bad$measurement_fit_eligibility$status[1, 1] <- "ineligible"
  bad$measurement_fit_eligibility$reason[1, 1] <- "user_predeclared_excluded"
  expect_error(validate_lt_rate(bad), "Layer-2.*identity")
})

test_that("FIX001-B embedded Layer-2 reason mutation fails closed", {
  z <- fix001_objects(); bad <- z$add
  bad$measurement_fit_eligibility$reason[1, 1] <- "measurement_qc_excluded"
  expect_error(validate_lt_rate(bad), "Layer-2.*identity")
})

test_that("FIX001-C through G rate semantic masquerades fail closed", {
  z <- fix001_objects()

  bad <- z$add; bad$representation <- "GBI"
  expect_error(validate_lt_rate(bad), "masquerade|semantic identity")

  bad <- z$add; bad$baseline_estimand <- "C2_positive_cell_log_OLS"
  expect_error(validate_lt_rate(bad), "masquerade|semantic identity")

  bad <- z$add; bad$baseline_fit_id <- "different_fit"
  expect_error(validate_lt_rate(bad), "semantic identity")

  bad <- z$add; bad$baseline_mu_sha256 <- strrep("0", 64)
  expect_error(validate_lt_rate(bad), "semantic identity")

  bad <- z$add; bad$rate_spec$method <- "GBI"
  expect_error(validate_lt_rate(bad), "specification hash|masquerade")
})

test_that("FIX001-H and I C2 factor mutations fail closed", {
  z <- fix001_objects()
  bad <- z$c2; bad$gene_factor[1] <- bad$gene_factor[1] * 2
  expect_error(validate_lt_c2_fit(bad), "effect/factor semantic identity")

  bad <- z$c2; bad$branch_factor[1] <- bad$branch_factor[1] * 2
  expect_error(validate_lt_c2_fit(bad), "effect/factor semantic identity")
})

test_that("FIX001-J and K C2 estimand run control and certificate are bound", {
  z <- fix001_objects()
  bad <- z$c2; bad$sensitivity_estimand <- "other_estimand"
  expect_error(validate_lt_c2_fit(bad), "estimand identity")

  bad <- z$c2; bad$run_id <- "other_run"
  expect_error(validate_lt_c2_fit(bad), "semantic identity")

  bad <- z$c2; bad$control$tolerance <- bad$control$tolerance * 2
  expect_error(validate_lt_c2_fit(bad), "normalized frozen control|certificate")

  bad <- z$c2; bad$certificates[[1]]$component_id <- "other_component"
  expect_error(validate_lt_c2_fit(bad), "certificate component")
})

test_that("FIX001-L through N common-domain ledgers are self-validating", {
  z <- fix001_objects()
  bad <- z$common; bad$common_linear <- bad$common_linear[-1]
  expect_error(validate_lt_rate_common_domain(bad), "ledger/count/hash")

  bad <- z$common; bad$common_cell_sha256 <- strrep("0", 64)
  expect_error(validate_lt_rate_common_domain(bad), "keyed semantic identity")

  bad <- z$common; bad$ordered_gene_ledger <- rev(bad$ordered_gene_ledger)
  expect_error(validate_lt_rate_common_domain(bad), "keyed semantic identity")

  bad <- z$common; bad$ordered_branch_ledger <- rev(bad$ordered_branch_ledger)
  expect_error(validate_lt_rate_common_domain(bad), "keyed semantic identity")
})

test_that("FIX001 rate C2 fit and common domain round trips validate", {
  z <- fix001_objects()
  for (object in list(z$add, z$c2, z$common)) {
    path <- tempfile(fileext = ".rds")
    lt_write_object(object, path)
    observed <- lt_read_object(path)
    expect_identical(observed, object)
  }
})

test_that("runtime usability accepts an ordinary user matrix without prior SHA", {
  x <- fix001_matrix()
  expect_invisible(validate_lt_matrix(x))
  expect_true(is.na(x$coordinate_provenance$reference_tree_sha256))
  expect_true(is.na(x$coordinate_provenance$ordered_taxon_ledger_sha256))
  expect_match(x$coordinate_provenance$audit_identity_status, "UNFROZEN")

  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(lt_c3_domain(x, eligibility), run_id = "ordinary_user")
  rate <- lt_rate(x, fit, lt_rate_spec("ordinary_add", "ADD_LT"), eligibility)
  expect_invisible(validate_lt_rate(rate))
})

test_that("display labels and non-authoritative metadata remain editable", {
  z <- fix001_objects()
  matrix_fingerprint <- lt_audit_fingerprint(z$x)
  edited <- lt_set_display_labels(
    z$x, c("display Gene A", "display Gene B"),
    c("display Branch A", "display Branch B", "display Branch C")
  )
  edited$metadata$user_note <- "ordinary mutable annotation"
  edited$coordinate_provenance$source_method <- "user relabelled import note"
  expect_invisible(validate_lt_matrix(edited))
  expect_false(lt_same_frozen_snapshot(edited, matrix_fingerprint))
  status <- lt_rate_dependency_status(z$add, edited, z$c3, z$eligibility)
  expect_identical(status$state, "CURRENT")

  rate_fingerprint <- lt_audit_fingerprint(z$add)
  semantic_before <- z$add$output_hashes$semantic_object_sha256
  edited_rate <- lt_set_display_labels(
    z$add, c("A", "B"), c("I", "II", "III")
  )
  edited_rate$metadata$user_note <- "display-only annotation"
  expect_invisible(validate_lt_rate(edited_rate))
  expect_identical(
    edited_rate$output_hashes$semantic_object_sha256, semantic_before
  )
  expect_false(lt_same_frozen_snapshot(edited_rate, rate_fingerprint))

  edited_common <- lt_set_display_labels(
    z$common, c("Common A", "Common B"), c("X", "Y", "Z")
  )
  expect_invisible(validate_lt_rate_common_domain(edited_common))
  expect_identical(
    edited_common$semantic_common_domain_sha256,
    z$common$semantic_common_domain_sha256
  )
})

test_that("rate-spec metadata does not change its scientific identity", {
  spec <- lt_rate_spec("editable_note", "GBI", metadata = list(note = "v1"))
  identity <- spec$rate_spec_id
  spec$metadata$note <- "v2"
  spec$metadata$display_group <- "user choice"
  expect_invisible(validate_lt_rate_spec(spec))
  expect_identical(spec$rate_spec_id, identity)
})

test_that("authoritative input edits stale dependencies without invalidating input", {
  z <- fix001_objects()
  edited <- z$x
  edited$values[1, 1] <- edited$values[1, 1] + 0.5
  expect_invisible(validate_lt_matrix(edited))
  status <- lt_rate_dependency_status(z$add, edited, z$c3, z$eligibility)
  expect_identical(status$state, "STALE_RECOMPUTE_REQUIRED")
  expect_match(status$reasons, "source_matrix")
  expect_error(
    validate_lt_rate_dependency(z$add, edited, z$c3, z$eligibility),
    "STALE_RECOMPUTE_REQUIRED"
  )
})

test_that("scientific-key edits require recomputation not global invalidity", {
  z <- fix001_objects()
  edited <- lt_matrix(
    z$x$values, z$x$coord_state, z$x$payload,
    gene_ids = c("new_gene_key_1", "new_gene_key_2"),
    branch_ids = z$x$branch_ids
  )
  eligibility <- lt_measurement_eligibility(edited)
  expect_invisible(validate_lt_matrix(edited))
  status <- lt_rate_dependency_status(z$add, edited, z$c3, eligibility)
  expect_identical(status$state, "STALE_RECOMPUTE_REQUIRED")
})

test_that("audit snapshot mismatch is distinct from structural invalidity", {
  z <- fix001_objects()
  fingerprint <- lt_audit_fingerprint(z$add)
  edited <- z$add
  edited$metadata$review_note <- "changed after snapshot"
  expect_invisible(validate_lt_rate(edited))
  expect_false(lt_same_frozen_snapshot(edited, fingerprint))
  expect_error(
    validate_lt_frozen_snapshot(edited, fingerprint),
    "NOT_THE_SAME_FROZEN_SNAPSHOT.*structurally valid"
  )
})
