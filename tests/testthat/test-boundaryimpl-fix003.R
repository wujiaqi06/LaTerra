fix003_fixture <- function() {
  values <- matrix(
    c(0, 1, 2, 0, 1, 0, 3, 4, 2, 3, 0, 1),
    3, 4, byrow = TRUE,
    dimnames = list(paste0("f3_g", 1:3), paste0("f3_b", 1:4))
  )
  coord <- matrix("observed", 3, 4, dimnames = dimnames(values))
  x <- lt_matrix(
    values, coord, fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(
    lt_c3_domain(x, eligibility), run_id = "boundaryimpl_fix003_fixture"
  )
  gbi <- lt_rate(
    x, fit, lt_rate_spec("boundaryimpl_fix003_GBI", "GBI"), eligibility
  )
  state <- lt_branch_state(
    "boundaryimpl_fix003_binary", "boundaryimpl_fix003_trait",
    colnames(values), c("reference", "reference", "focal", "focal"),
    "binary", rep("available", 4), rep(NA_character_, 4),
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
    dimnames = list(colnames(values), paste0("f3_r", 1:3))
  )
  null <- lt_branch_state_ensemble(
    "boundaryimpl_fix003_null", "boundaryimpl_fix003_trait",
    colnames(values), colnames(null_values), null_values, "binary",
    coding = list(reference_level = "reference", focal_level = "focal"),
    availability = matrix(
      "available", 4, 3, dimnames = dimnames(null_values)
    ),
    null_hypothesis = "fixed FIX003 unit-test null",
    generator_provenance = list(generator_type = "fixed_fixture")
  )
  list(x = x, eligibility = eligibility, fit = fit, gbi = gbi,
       state = state, null = null)
}

fix003_source_root <- function() {
  test_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  candidates <- c(test_root, file.path(test_root, "00_pkg_src", "LaTerra"))
  is_source_root <- vapply(candidates, function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      dir.exists(file.path(path, "R")) &&
      dir.exists(file.path(path, "inst", "schema"))
  }, logical(1L))
  if (!any(is_source_root)) stop("Cannot locate package source root.")
  normalizePath(candidates[which(is_source_root)[[1L]]], mustWork = TRUE)
}

fix003_rehash_encoded <- function(x, rehash_source_claims = TRUE) {
  state <- as.integer(x$ratio_state)
  positive <- matrix(
    state == 1L, nrow(x$values), ncol(x$values),
    dimnames = dimnames(x$values)
  )
  zero <- matrix(
    state == 2L, nrow(x$values), ncol(x$values),
    dimnames = dimnames(x$values)
  )
  available <- positive | zero
  x$positive_log_values_identity <-
    LaTerra:::.lt_legacy_encoded_positive_identity(
      x$values, x$ratio_state, x$ordered_gene_ledger,
      x$ordered_branch_ledger
    )
  x$positive_domain_identity <- LaTerra:::.lt_rate_keyed_mask_hash(
    positive, x$ordered_gene_ledger, x$ordered_branch_ledger
  )
  x$exact_zero_identity <- LaTerra:::.lt_rate_keyed_mask_hash(
    zero, x$ordered_gene_ledger, x$ordered_branch_ledger
  )
  x$view_domain_identity <- LaTerra:::.lt_rate_keyed_mask_hash(
    available, x$ordered_gene_ledger, x$ordered_branch_ledger
  )
  x$view_ordered_axes_identity <- LaTerra:::.lt_hash(list(
    gene_ids = x$ordered_gene_ledger,
    branch_ids = x$ordered_branch_ledger
  ))
  if (rehash_source_claims) {
    x$source_semantic_identity <- LaTerra:::.lt_hash(list(
      attacker_payload = x$values, attacker_state = x$ratio_state
    ))
    x$source_audit_identity <- LaTerra:::.lt_hash(list(
      attacker_semantic = x$source_semantic_identity,
      attacker_m_ref = x$m_ref
    ))
    x$source_positive_log_values_sha256 <- LaTerra:::.lt_hash(x$values)
    x$source_ratio_state_sha256 <- LaTerra:::.lt_hash(x$ratio_state)
    x$source_positive_domain_identity <- x$positive_domain_identity
    x$source_exact_zero_identity <- x$exact_zero_identity
    x$source_object_identity <-
      LaTerra:::.lt_legacy_source_object_identity(x)
  }
  x$legacy_mode_identity <- LaTerra:::.lt_hash(list(
    mode = "legacy_global_k10",
    contract = x$encoding_contract_version,
    m_ref = x$m_ref,
    L_min_plus = x$L_min_plus,
    sentinel = x$sentinel
  ))
  x$m_ref_provenance <- LaTerra:::.lt_legacy_m_ref_provenance(x)
  x$m_ref_provenance_identity <- LaTerra:::.lt_hash(x$m_ref_provenance)
  x$encoded_identity <- LaTerra:::.lt_hash(
    LaTerra:::.lt_legacy_encoded_payload(x)
  )
  x
}

fix003_rehash_log_relative <- function(x) {
  x$semantic_identity <- LaTerra:::.lt_log_relative_identity_hash(
    LaTerra:::.lt_log_relative_semantic_payload(x),
    x$canonicalization_version
  )
  x$audit_identity <- LaTerra:::.lt_log_relative_identity_hash(
    LaTerra:::.lt_log_relative_audit_payload(x),
    x$canonicalization_version
  )
  x
}

fix003_replay_authority <- function(encoded, authority_id = "fixed_test_authority") {
  list(
    authority_id = authority_id,
    purpose = "HISTORICAL_REPLAY_VALIDATION",
    source_semantic_identity = encoded$source_semantic_identity,
    source_audit_identity = encoded$source_audit_identity,
    source_object_identity = encoded$source_object_identity,
    encoded_identity = encoded$encoded_identity,
    mode = encoded$mode,
    m_ref = encoded$m_ref,
    L_min_plus = encoded$L_min_plus,
    sentinel = encoded$sentinel,
    positive_log_values_identity = encoded$positive_log_values_identity,
    positive_domain_identity = encoded$positive_domain_identity,
    exact_zero_identity = encoded$exact_zero_identity,
    view_domain_identity = encoded$view_domain_identity,
    view_ordered_axes_identity = encoded$view_ordered_axes_identity
  )
}

test_that("FIX003-001 no encoded authorization token or seal implementation remains", {
  root <- fix003_source_root()
  r_files <- list.files(file.path(root, "R"), pattern = "[.]R$",
                        full.names = TRUE)
  text <- paste(unlist(lapply(r_files, readLines, warn = FALSE)),
                collapse = "\n")
  expect_false(grepl("lt_encoded_authorization_seal", text, fixed = TRUE))
  expect_false(grepl("lt_authorize_encoded_inference", text, fixed = TRUE))
  expect_false(grepl("lt_encoded_inference_authorization", text, fixed = TRUE))
})

test_that("FIX003-002 every ordinary and internal encoded inference route refuses", {
  z <- fix003_fixture()
  encoded <- logGBI_encoded(lt_log_relative(z$gbi))
  refusal <- "lt_error_legacy_encoded_inference_not_supported"

  expect_error(lt_association_domain(encoded, z$state), class = refusal)
  expect_error(
    lt_common_association_domain(list(encoded, z$gbi), z$state),
    class = refusal
  )
  expect_error(LaTerra:::.lt_rate_association_identity(encoded),
               class = refusal)
  expect_error(
    LaTerra:::.lt_domain_rates(
      list(encoded), z$state, NULL, "error", "forbidden", "forbidden"
    ), class = refusal
  )
  expect_error(lt_associate_univariate(encoded, z$state), class = refusal)
  expect_error(
    lt_associate_univariate(
      encoded, z$state, encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ), class = refusal
  )
  expect_error(
    LaTerra:::.lt_associate_legacy_encoded(
      encoded, z$state, "forged", list(record = "forged")
    ), class = refusal
  )
  expect_error(
    lt_calibrate_association(encoded, NULL, NULL, NULL, NULL),
    class = refusal
  )
  expect_error(
    LaTerra:::.lt_null_score(encoded, NULL, NULL, "binary", NULL, "beta"),
    class = refusal
  )
  expect_error(
    LaTerra:::.lt_calibrate_legacy_encoded_null(
      encoded, z$state, NULL, z$null, "NULL_A", "two_sided", "forged",
      list(record = "forged")
    ), class = refusal
  )
  expect_error(
    lt_calibrate_logGBI_encoded(
      encoded, z$state, NULL, list(NULL_A = z$null),
      encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ), class = refusal
  )
  expect_error(lt_adjust_logGBI_encoded(list()), class = refusal)
  expect_error(stats::t.test(encoded), class = refusal)
})

test_that("FIX003-003 structural and source-bound validation are separate", {
  z <- fix003_fixture()
  source <- lt_log_relative(z$gbi)
  encoded <- logGBI_encoded(source)
  expect_invisible(validate_lt_logGBI_encoded(encoded))
  expect_invisible(validate_lt_logGBI_encoded_source(encoded, source))

  positive_tamper <- encoded
  positive_index <- which(as.integer(positive_tamper$ratio_state) == 1L)[[1L]]
  positive_tamper$values[positive_index] <-
    positive_tamper$values[positive_index] + 0.125
  positive_tamper <- fix003_rehash_encoded(positive_tamper)
  expect_invisible(validate_lt_logGBI_encoded(positive_tamper))
  expect_error(
    validate_lt_logGBI_encoded_source(positive_tamper, source),
    class = "lt_error_legacy_encoded_source_mismatch"
  )

  mref_tamper <- encoded
  mref_tamper$m_ref <- encoded$m_ref / 2
  mref_tamper$L_min_plus <- log(mref_tamper$m_ref)
  mref_tamper$sentinel <- mref_tamper$L_min_plus - 10
  mref_tamper$values[as.integer(mref_tamper$ratio_state) == 2L] <-
    mref_tamper$sentinel
  mref_tamper <- fix003_rehash_encoded(mref_tamper)
  expect_invisible(validate_lt_logGBI_encoded(mref_tamper))
  expect_error(
    validate_lt_logGBI_encoded_source(mref_tamper, source),
    class = "lt_error_legacy_encoded_source_mismatch"
  )

  state_tamper <- encoded
  zero_index <- which(as.integer(state_tamper$ratio_state) == 2L)[[1L]]
  state_tamper$ratio_state[zero_index] <- as.raw(1L)
  state_tamper$values[zero_index] <- state_tamper$L_min_plus + 1
  state_tamper <- fix003_rehash_encoded(state_tamper)
  expect_invisible(validate_lt_logGBI_encoded(state_tamper))
  expect_error(
    validate_lt_logGBI_encoded_source(state_tamper, source),
    class = "lt_error_legacy_encoded_source_mismatch"
  )
})

test_that("FIX003-004 subset validation preserves complete-source m_ref", {
  source <- lt_log_relative(fix003_fixture()$gbi)
  state <- ratio_state(source)
  positive <- which(state == "POSITIVE", arr.ind = TRUE)
  minimum <- positive[which.min(
    positive_log_values(source)[state == "POSITIVE"]
  ), , drop = FALSE]
  genes <- source$ordered_gene_ledger[-minimum[1L, 1L]]
  if (!length(genes)) genes <- rev(source$ordered_gene_ledger)
  branches <- rev(source$ordered_branch_ledger)
  encoded <- logGBI_encoded(source, gene_ids = genes, branch_ids = branches)

  expect_invisible(validate_lt_logGBI_encoded_source(encoded, source))
  expect_identical(
    encoded$m_ref, source$representation_provenance$original_source_m_ref
  )
  expect_identical(
    encoded$sentinel,
    source$representation_provenance$legacy_global_k10_sentinel
  )
  expect_false(encoded$provenance$m_ref_recomputed_after_subset)
})

test_that("FIX003-005 replay validation needs fixed external authority", {
  source <- lt_log_relative(fix003_fixture()$gbi)
  encoded <- logGBI_encoded(source)
  authority <- fix003_replay_authority(encoded)
  result <- validate_legacy_logGBI_replay(encoded, source, authority)
  expect_identical(result$status, "EXACT_FROZEN_REPLAY_MATCH")
  expect_identical(result$purpose, "HISTORICAL_REPLAY_VALIDATION")
  expect_false(result$inference_authorized)

  changed_source <- source
  positive <- which(as.integer(changed_source$ratio_state) == 1L)
  nonminimum <- positive[which.max(changed_source$positive_log_values[positive])]
  changed_source$positive_log_values[nonminimum] <-
    changed_source$positive_log_values[nonminimum] + 0.25
  changed_source <- fix003_rehash_log_relative(changed_source)
  expect_invisible(validate_lt_log_relative(changed_source))
  changed_encoded <- logGBI_encoded(changed_source)
  expect_invisible(validate_lt_logGBI_encoded_source(
    changed_encoded, changed_source
  ))
  expect_error(
    validate_legacy_logGBI_replay(
      changed_encoded, changed_source, authority
    ), class = "lt_error_legacy_replay_authority_mismatch"
  )
})

test_that("FIX003-006 replay authority cannot be replaced by arbitrary metadata", {
  source <- lt_log_relative(fix003_fixture()$gbi)
  encoded <- logGBI_encoded(source)
  bad <- fix003_replay_authority(encoded)
  bad$purpose <- "USER_AUTHORIZED"
  expect_error(
    validate_legacy_logGBI_replay(encoded, source, bad),
    class = "lt_error_legacy_replay_authority_invalid"
  )
  expect_error(
    validate_legacy_logGBI_replay(encoded, source, list()),
    class = "lt_error_legacy_replay_authority_invalid"
  )
})

test_that("FIX003-007 documentation scanner rejects current acknowledgement unlock claims", {
  current <- tempfile(pattern = "README_", fileext = ".md")
  writeLines(
    "Acknowledgement authorizes finite-sentinel encoded inference.", current
  )
  expect_error(
    LaTerra:::.lt_assert_loggbi_documentation_conforms(
      dirname(current), current_version = "0.0.0.9017", paths = current
    ), class = "lt_error_logGBI_documentation_nonconformant"
  )
  historical <- tempfile(pattern = "NEWS_", fileext = ".md")
  writeLines(c(
    "# LaTerra 0.0.0.9016",
    "Acknowledgement authorizes finite-sentinel encoded inference."
  ), historical)
  expect_equal(nrow(LaTerra:::.lt_scan_loggbi_documentation(
    dirname(historical), current_version = "0.0.0.9017",
    paths = historical
  )), 0L)
})
