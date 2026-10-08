fix002_fixture <- function() {
  values <- matrix(
    c(0, 1, 2, 0, 1, 0, 3, 4, 2, 3, 0, 1),
    3, 4, byrow = TRUE,
    dimnames = list(paste0("f2_g", 1:3), paste0("f2_b", 1:4))
  )
  coord <- matrix("observed", 3, 4, dimnames = dimnames(values))
  x <- lt_matrix(
    values, coord, fixture_payload(),
    coordinate_provenance = fixture_coordinate_provenance()
  )
  eligibility <- lt_measurement_eligibility(x)
  fit <- lt_c3_fit(
    lt_c3_domain(x, eligibility), run_id = "boundaryimpl_fix002_fixture"
  )
  gbi <- lt_rate(
    x, fit, lt_rate_spec("boundaryimpl_fix002_GBI", "GBI"), eligibility
  )
  state <- lt_branch_state(
    "boundaryimpl_fix002_binary", "boundaryimpl_fix002_trait",
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
    dimnames = list(colnames(values), paste0("f2_r", 1:3))
  )
  null <- lt_branch_state_ensemble(
    "boundaryimpl_fix002_null", "boundaryimpl_fix002_trait",
    colnames(values), colnames(null_values), null_values, "binary",
    coding = list(reference_level = "reference", focal_level = "focal"),
    availability = matrix(
      "available", 4, 3, dimnames = dimnames(null_values)
    ),
    null_hypothesis = "fixed FIX002 unit-test null",
    generator_provenance = list(generator_type = "fixed_fixture")
  )
  list(x = x, eligibility = eligibility, fit = fit, gbi = gbi,
       state = state, null = null)
}

fix002_invalid_config <- function() {
  path <- system.file("recipes", "marine_2026_legacy.yaml", package = "LaTerra")
  stopifnot(nzchar(path))
  config <- lt_read_config(path, validate = FALSE)
  config$rate$method <- "logGBI"
  config
}

fix002_rehash_encoded <- function(x) {
  x$encoded_identity <- LaTerra:::.lt_hash(
    LaTerra:::.lt_legacy_encoded_payload(x)
  )
  x
}

fix002_source_root <- function() {
  test_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
  candidates <- c(
    test_root,
    file.path(test_root, "00_pkg_src", "LaTerra")
  )
  is_source_root <- vapply(candidates, function(path) {
    file.exists(file.path(path, "DESCRIPTION")) &&
      dir.exists(file.path(path, "R")) &&
      dir.exists(file.path(path, "man")) &&
      dir.exists(file.path(path, "inst", "schema")) &&
      dir.exists(file.path(path, "inst", "recipes"))
  }, logical(1))
  if (!any(is_source_root)) {
    stop("Cannot locate the frozen package source root for conformance tests.")
  }
  normalizePath(candidates[which(is_source_root)[[1L]]], mustWork = TRUE)
}

test_that("FIX002-001 lt_run validates every attached config route", {
  invalid <- fix002_invalid_config()
  expect_error(
    lt_run("hand_built_invalid", config = invalid),
    "lt_migrate_logGBI",
    class = "lt_error_logGBI_constructor_moved"
  )

  config_rds <- tempfile(fileext = ".rds")
  saveRDS(invalid, config_rds, version = 3)
  restored_config <- readRDS(config_rds)
  expect_error(
    lt_run("deserialized_invalid", config = restored_config),
    class = "lt_error_logGBI_constructor_moved"
  )

  profile_yaml <- tempfile(fileext = ".yaml")
  lt_write_config(invalid, profile_yaml, validate = FALSE)
  expect_error(
    lt_read_config(profile_yaml),
    class = "lt_error_logGBI_constructor_moved"
  )
  profile <- lt_read_config(profile_yaml, validate = FALSE)
  expect_error(
    lt_run("profile_invalid", config = profile),
    class = "lt_error_logGBI_constructor_moved"
  )
})

test_that("FIX002-002 invalid configured runs cannot be written or resurrected", {
  forged <- lt_run("forged_invalid_run")
  forged$config <- fix002_invalid_config()
  expect_error(
    validate_lt_run(forged),
    class = "lt_error_logGBI_constructor_moved"
  )
  path <- tempfile(fileext = ".rds")
  expect_error(
    lt_write_object(forged, path),
    class = "lt_error_logGBI_constructor_moved"
  )
  saveRDS(forged, path, version = 3)
  expect_error(
    lt_read_object(path),
    class = "lt_error_logGBI_constructor_moved"
  )
})

test_that("FIX002-003 encoded inference is unconditionally refused", {
  z <- fix002_fixture()
  encoded <- logGBI_encoded(lt_log_relative(z$gbi))
  expect_error(
    lt_associate_univariate(
      encoded, z$state, encoded_mode = "legacy_global_k10"
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
  expect_error(
    LaTerra:::.lt_associate_legacy_encoded(
      encoded, z$state, "direct_without_authorization"
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
  expect_error(
    LaTerra:::.lt_associate_legacy_encoded(
      encoded, z$state, "forged_authorization", list(forged = TRUE)
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
  expect_error(
    LaTerra:::.lt_calibrate_legacy_encoded_null(
      encoded, z$state, NULL, z$null, "NULL_A", "two_sided",
      "direct_calibration_without_authorization"
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )
})

test_that("FIX002-004 acknowledgement and mode cannot unlock calibration", {
  z <- fix002_fixture()
  boundary <- lt_log_relative(z$gbi)
  encoded <- logGBI_encoded(boundary)
  expect_invisible(validate_lt_logGBI_encoded_source(encoded, boundary))
  expect_error(
    lt_calibrate_logGBI_encoded(
      encoded, z$state, NULL, list(NULL_A = z$null),
      encoded_mode = "legacy_global_k10",
      acknowledge_encoded_zero = TRUE
    ),
    class = "lt_error_legacy_encoded_inference_not_supported"
  )

  r_files <- list.files(file.path(fix002_source_root(), "R"),
                        pattern = "[.]R$", full.names = TRUE)
  lines <- unlist(lapply(r_files, readLines, warn = FALSE), use.names = FALSE)
  expect_false(any(grepl("lt_encoded_authorization_seal", lines, fixed = TRUE)))
  expect_false(any(grepl("lt_authorize_encoded_inference", lines, fixed = TRUE)))
})

test_that("FIX002-005 legacy encoded mathematical corruption fails closed", {
  encoded <- logGBI_encoded(lt_log_relative(fix002_fixture()$gbi))
  expect_invisible(validate_lt_logGBI_encoded(encoded))

  corruptions <- list(
    m_ref = function(x) { x$m_ref <- x$m_ref * 2; x },
    L_min_plus = function(x) { x$L_min_plus <- x$L_min_plus + 1; x },
    sentinel = function(x) { x$sentinel <- x$sentinel - 1; x },
    zero_value = function(x) {
      i <- which(as.integer(x$ratio_state) == 2L)[[1L]]
      x$values[i] <- x$sentinel + 1
      x
    },
    positive_value = function(x) {
      i <- which(as.integer(x$ratio_state) == 1L)[[1L]]
      x$values[i] <- x$values[i] + 1
      x
    },
    positive_identity = function(x) {
      x$positive_domain_identity <- strrep("a", 64L)
      x
    },
    zero_identity = function(x) {
      x$exact_zero_identity <- strrep("b", 64L)
      x
    },
    source_identity = function(x) {
      x$source_object_identity <- strrep("c", 64L)
      x
    },
    legacy_mode_identity = function(x) {
      x$legacy_mode_identity <- strrep("d", 64L)
      x
    },
    m_ref_provenance_identity = function(x) {
      x$m_ref_provenance_identity <- strrep("e", 64L)
      x
    }
  )
  for (mutate in corruptions) {
    bad <- fix002_rehash_encoded(mutate(encoded))
    expect_error(validate_lt_logGBI_encoded(bad))
  }
})

test_that("FIX002-006 unavailable cells cannot impersonate encoded zero", {
  encoded <- logGBI_encoded(lt_log_relative(fix002_fixture()$gbi))
  bad <- encoded
  i <- which(as.integer(bad$ratio_state) == 1L)[[1L]]
  bad$ratio_state[i] <- as.raw(3L)
  bad$values[i] <- bad$sentinel
  bad <- fix002_rehash_encoded(bad)
  expect_error(validate_lt_logGBI_encoded(bad), "unavailable|state")
})

test_that("FIX002-007 real public subset preserves full-source m_ref", {
  z <- fix002_fixture()
  boundary <- lt_log_relative(z$gbi)
  positive <- which(ratio_state(boundary) == "POSITIVE", arr.ind = TRUE)
  min_row <- positive[which.min(positive_log_values(boundary)[
    ratio_state(boundary) == "POSITIVE"
  ]), 1L]
  keep_genes <- boundary$ordered_gene_ledger[-min_row]
  encoded <- logGBI_encoded(
    boundary, mode = "legacy_global_k10", gene_ids = keep_genes
  )

  expect_identical(encoded$m_ref,
                   boundary$representation_provenance$original_source_m_ref)
  expect_identical(encoded$L_min_plus,
                   boundary$representation_provenance$L_min_plus)
  expect_identical(
    encoded$sentinel,
    boundary$representation_provenance$legacy_global_k10_sentinel
  )
  state <- as.integer(encoded$ratio_state)
  expect_true(min(encoded$values[state == 1L]) > encoded$L_min_plus)
  expect_true(all(encoded$values[state == 2L] == encoded$sentinel))
  expect_true(encoded$provenance$subset_applied)
  expect_false(encoded$provenance$m_ref_recomputed_after_subset)

  expect_invisible(validate_lt_logGBI_encoded_source(encoded, boundary))
  path <- tempfile(fileext = ".rds")
  lt_write_object(encoded, path)
  expect_identical(lt_read_object(path), encoded)
})

test_that("FIX002-008 negative documentation scanner covers current materials", {
  root <- fix002_source_root()
  paths <- LaTerra:::.lt_loggbi_doc_files(root)
  relative <- substring(paths, nchar(root) + 2L)
  expect_true("DESCRIPTION" %in% relative)
  expect_true(any(grepl("^README", relative)))
  expect_true(any(grepl("^NEWS", relative)))
  expect_true(any(grepl("^man/", relative)))
  expect_true(any(grepl("^inst/doc/", relative)))
  expect_true(any(grepl("^inst/schema/", relative)))
  expect_true(any(grepl("^inst/recipes/", relative)))
  expect_invisible(LaTerra:::.lt_assert_loggbi_documentation_conforms(
    root, current_version = "0.0.0.9017"
  ))
})

test_that("FIX002-009 scanner catches current claims but preserves versioned history", {
  current <- tempfile(pattern = "README_", fileext = ".md")
  writeLines(
    "Default production logGBI uses a finite zero-state sentinel.", current
  )
  expect_error(
    LaTerra:::.lt_assert_loggbi_documentation_conforms(
      dirname(current), current_version = "0.0.0.9017", paths = current
    ),
    class = "lt_error_logGBI_documentation_nonconformant"
  )

  historical <- tempfile(pattern = "NEWS_", fileext = ".md")
  writeLines(c(
    "# LaTerra 0.0.0.9012",
    "Default production logGBI uses a finite zero-state sentinel."
  ), historical)
  expect_equal(nrow(LaTerra:::.lt_scan_loggbi_documentation(
    dirname(historical), current_version = "0.0.0.9017",
    paths = historical
  )), 0L)

  current_news <- tempfile(pattern = "NEWS_", fileext = ".md")
  writeLines(c(
    "# LaTerra 0.0.0.9017",
    "Finite zero-state logGBI is the current primary representation."
  ), current_news)
  expect_gt(nrow(LaTerra:::.lt_scan_loggbi_documentation(
    dirname(current_news), current_version = "0.0.0.9017",
    paths = current_news
  )), 0L)
})
