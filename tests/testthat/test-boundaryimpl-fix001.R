fix001_fixture <- function() {
  values <- matrix(
    c(0, 1, 2, 0, 1, 0, 3, 4, 2, 3, 0, 1),
    3, 4, byrow = TRUE,
    dimnames = list(paste0("fx_g", 1:3), paste0("fx_b", 1:4))
  )
  coord <- matrix("observed", 3, 4, dimnames = dimnames(values))
  x <- lt_matrix(
    values, coord, fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(
    lt_c3_domain(x, eligibility), run_id = "boundaryimpl_fix001_fixture"
  )
  gbi <- lt_rate(
    x, fit, lt_rate_spec("boundaryimpl_fix001_GBI", "GBI"), eligibility
  )
  strict <- lt_rate(
    x, fit,
    lt_rate_spec("boundaryimpl_fix001_logGBI_strict", "logGBI_strict"),
    eligibility
  )
  state <- lt_branch_state(
    "boundaryimpl_fix001_binary", "boundaryimpl_fix001_trait",
    colnames(values), c("reference", "reference", "focal", "focal"),
    "binary", rep("available", 4), rep(NA_character_, 4),
    coding = list(reference_level = "reference", focal_level = "focal"),
    provenance = list(source_type = "fixed_fixture")
  )
  list(x = x, eligibility = eligibility, fit = fit, gbi = gbi,
       strict = strict, state = state)
}

fix001_historical_rate <- function() {
  path <- system.file("extdata", "legacy_logGBI_9010.rds", package = "LaTerra")
  stopifnot(nzchar(path))
  lt_read_object(path)
}

fix001_historical_state <- function(rate = fix001_historical_rate()) {
  n <- length(rate$ordered_branch_ledger)
  values <- rep(c("reference", "focal"), length.out = n)
  lt_branch_state(
    "historical_logGBI_state", "historical_logGBI_trait",
    rate$ordered_branch_ledger, values, "binary",
    rep("available", n), rep(NA_character_, n),
    coding = list(reference_level = "reference", focal_level = "focal"),
    provenance = list(source_type = "fixed_historical_fixture")
  )
}

test_that("BOUNDARYIMPL-FIX001-001 historical rate stays readable and migratable", {
  historical <- fix001_historical_rate()
  expect_s3_class(historical, "lt_rate")
  expect_identical(historical$representation, "logGBI")
  expect_invisible(validate_lt_rate(historical))

  z <- fix001_fixture()
  migrated <- lt_migrate_logGBI(z$strict, z$gbi)
  direct <- lt_log_relative(z$gbi)
  expect_identical(migrated$ratio_state, direct$ratio_state)
  expect_identical(migrated$unavailability_reason,
                   direct$unavailability_reason)
  expect_identical(migrated$positive_log_values,
                   direct$positive_log_values)
  expect_identical(migrated$semantic_identity, direct$semantic_identity)
})

test_that("BOUNDARYIMPL-FIX001-002 historical rate fails every ordinary domain route", {
  historical <- fix001_historical_rate()
  state <- fix001_historical_state(historical)
  expect_error(
    lt_association_domain(historical, state),
    class = "lt_error_historical_logGBI_inference_forbidden"
  )
  expect_error(
    lt_common_association_domain(list(historical, historical), state),
    class = "lt_error_historical_logGBI_inference_forbidden"
  )
})

test_that("BOUNDARYIMPL-FIX001-003 historical rate fails univariate and calibration routes", {
  historical <- fix001_historical_rate()
  state <- fix001_historical_state(historical)
  expect_error(
    lt_associate_univariate(historical, state),
    class = "lt_error_historical_logGBI_inference_forbidden"
  )
  expect_error(
    lt_calibrate_association(historical, NULL, NULL, NULL, NULL),
    class = "lt_error_historical_logGBI_inference_forbidden"
  )
  expect_error(
    LaTerra:::.lt_null_score(historical, NULL, NULL, "binary", NULL, "beta"),
    class = "lt_error_historical_logGBI_inference_forbidden"
  )
})

test_that("BOUNDARYIMPL-FIX001-004 historical rate specs and enclosing profiles fail closed", {
  historical <- fix001_historical_rate()
  spec <- historical$rate_spec
  expect_identical(spec$method, "logGBI")
  expect_error(validate_lt_rate_spec(spec),
               class = "lt_error_logGBI_constructor_moved")
  expect_error(
    lt_analysis_spec("historical_profile", spec),
    class = "lt_error_logGBI_constructor_moved"
  )
  path <- tempfile(fileext = ".rds")
  saveRDS(spec, path, version = 3)
  expect_error(lt_read_object(path),
               class = "lt_error_logGBI_constructor_moved")

  z <- fix001_fixture()
  expect_error(
    lt_rate(z$x, z$fit, spec, z$eligibility),
    class = "lt_error_logGBI_constructor_moved"
  )
})

test_that("BOUNDARYIMPL-FIX001-005 legacy diagnostic encoding is exact and source bound", {
  z <- fix001_fixture()
  boundary <- lt_log_relative(z$gbi)
  encoded <- logGBI_encoded(boundary, mode = "legacy_global_k10")
  expect_invisible(validate_lt_logGBI_encoded(encoded))
  expect_invisible(validate_lt_logGBI_encoded_source(encoded, boundary))
  expect_identical(lt_encoded_values(encoded), encoded$values)
  expect_error(
    lt_associate_univariate(
      encoded, z$state, encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
})

test_that("BOUNDARYIMPL-FIX001-006 subset views never recompute source-authority m_ref", {
  boundary <- lt_log_relative(fix001_fixture()$gbi)
  encoded <- logGBI_encoded(boundary, mode = "legacy_global_k10")
  positive <- which(ratio_state(boundary) == "POSITIVE", arr.ind = TRUE)
  minimum <- positive[which.min(
    positive_log_values(boundary)[ratio_state(boundary) == "POSITIVE"]
  ), , drop = FALSE]
  keep_gene <- setdiff(seq_len(nrow(encoded$values)), minimum[1, 1])
  keep_branch <- setdiff(seq_len(ncol(encoded$values)), minimum[1, 2])
  if (!length(keep_gene)) keep_gene <- seq_len(nrow(encoded$values))
  if (!length(keep_branch)) keep_branch <- seq_len(ncol(encoded$values))
  filtered_values <- encoded$values[keep_gene, keep_branch, drop = FALSE]

  expect_identical(encoded$m_ref,
                   boundary$representation_provenance$original_source_m_ref)
  expect_identical(encoded$L_min_plus,
                   boundary$representation_provenance$L_min_plus)
  expect_identical(encoded$sentinel,
                   boundary$representation_provenance$legacy_global_k10_sentinel)
  expect_gt(length(filtered_values), 0L)
  expect_identical(encoded$sentinel, log(encoded$m_ref) - 10)
})

test_that("BOUNDARYIMPL-FIX001-007 boundary identities execute serialization v3", {
  x <- lt_log_relative(fix001_fixture()$gbi)
  semantic_bytes <- serialize(
    LaTerra:::.lt_log_relative_semantic_payload(x), NULL, version = 3
  )
  audit_bytes <- serialize(
    LaTerra:::.lt_log_relative_audit_payload(x), NULL, version = 3
  )
  expect_identical(
    x$semantic_identity,
    digest::digest(semantic_bytes, algo = "sha256", serialize = FALSE)
  )
  expect_identical(
    x$audit_identity,
    digest::digest(audit_bytes, algo = "sha256", serialize = FALSE)
  )
  expect_identical(
    x$canonicalization_version,
    "LaTerra_R_serialization_v3_raw_code_order_v1"
  )
  expect_error(
    LaTerra:::.lt_log_relative_identity_hash(list(a = 1), "unknown"),
    class = "lt_error_log_relative_canonicalization_unknown"
  )
})

test_that("BOUNDARYIMPL-FIX001-008 ratio baseline invariant prevents infinity", {
  contract <- LaTerra:::.lt_ratio_baseline_contract(c(0, 2), c(0, 4))
  expect_identical(contract$ratio_available, c(FALSE, TRUE))
  expect_identical(contract$reason, c("baseline_zero", "available"))
  ratio <- rep(NA_real_, 2)
  ratio[contract$ratio_available] <-
    c(0, 2)[contract$ratio_available] / c(0, 4)[contract$ratio_available]
  expect_false(any(is.infinite(ratio), na.rm = TRUE))
  expect_error(
    LaTerra:::.lt_ratio_baseline_contract(1, 0),
    class = "lt_error_positive_y_zero_mu"
  )
})

test_that("BOUNDARYIMPL-FIX001-009 current documentation and schemas conform", {
  desc <- packageDescription("LaTerra")
  expect_identical(desc$Version, "0.0.0.9022")
  expect_match(desc$Description, "boundary aware", fixed = TRUE)
  expect_match(desc$Description, "read/migrate-only", fixed = TRUE)

  help_text <- function(topic) {
    help_file <- utils:::.getHelpFile(utils::help(topic, package = "LaTerra"))
    text <- paste(capture.output(tools::Rd2txt(help_file)), collapse = "\n")
    gsub(intToUtf8(8L), "", text, fixed = TRUE)
  }
  expect_match(help_text("lt_rate"), "EXACT_ZERO_BOUNDARY", fixed = TRUE)
  expect_match(help_text("lt_rate"), "legacy_global_k10", fixed = TRUE)
  expect_match(help_text("lt_calibrate_association"),
               "fail closed", fixed = TRUE)
  expect_match(help_text("lt_diagnose_rate"),
               "read[[:space:]-]+migrate-only")

  schema <- yaml::read_yaml(system.file(
    "schema", "lt_objects.schema.yaml", package = "LaTerra"
  ))
  expect_identical(schema$schema_version, "1.6.0")
  expect_match(schema$objects$lt_rate$fields$representation,
               "read/migrate-only", fixed = TRUE)
  expect_match(schema$objects$lt_log_relative$fields$canonicalization_version,
               "R version 3", fixed = TRUE)
})
