test_that("lt_matrix validates orthogonal matrix fields", {
  x <- small_lt_matrix()

  expect_s3_class(x, "lt_matrix")
  expect_identical(x$values[1, 1], 0)
  expect_identical(x$coord_state[1, 1], "observed")
  expect_true(all(is.na(x$value_reason)))
  expect_identical(x$payload$type, "branch_length")
  expect_invisible(validate_lt_matrix(x))
  expect_output(print(x), "2 genes x 2 branches")
})

test_that("duplicate gene IDs are rejected", {
  x <- small_lt_matrix()
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      value_reason = x$value_reason,
      gene_ids = c("gene_a", "gene_a"),
      branch_ids = x$branch_ids,
      coordinate_provenance = x$coordinate_provenance
    ),
    "duplicate identifiers"
  )
})

test_that("duplicate branch IDs are rejected", {
  x <- small_lt_matrix()
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      value_reason = x$value_reason,
      gene_ids = x$gene_ids,
      branch_ids = c("B1", "B1"),
      coordinate_provenance = x$coordinate_provenance
    ),
    "duplicate identifiers"
  )
})

test_that("matrix coordinate state and value reason dimensions must match", {
  x <- small_lt_matrix()
  bad_state <- x$coord_state[, 1, drop = FALSE]
  expect_error(
    lt_matrix(
      x$values, bad_state, x$payload,
      coordinate_provenance = x$coordinate_provenance
    ),
    "dimension mismatch"
  )

  bad_reason <- x$value_reason[, 1, drop = FALSE]
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      value_reason = bad_reason,
      coordinate_provenance = x$coordinate_provenance
    ),
    "dimension mismatch"
  )
})

test_that("coordinate-state vocabulary remains closed", {
  x <- small_lt_matrix()
  bad_state <- x$coord_state
  bad_state[2, 2] <- "rate_trimmed"
  expect_error(
    lt_matrix(
      x$values, bad_state, x$payload,
      value_reason = x$value_reason,
      coordinate_provenance = x$coordinate_provenance
    ),
    "Unknown coordinate-state labels"
  )
})

test_that("coordinate truth and current-payload absence remain separate", {
  x <- small_lt_matrix()
  values <- x$values
  reasons <- x$value_reason
  values[1, 2] <- NA_real_
  reasons[1, 2] <- "legacy_gbi_upper_tail_trim"

  derived <- lt_matrix(
    values, x$coord_state, fixture_payload("derived", "rate_representation"),
    value_reason = reasons,
    coordinate_provenance = x$coordinate_provenance
  )
  expect_identical(derived$coord_state[1, 2], "observed")
  expect_identical(derived$value_reason[1, 2], "legacy_gbi_upper_tail_trim")

  reasons[1, 2] <- NA_character_
  expect_error(
    lt_matrix(
      values, x$coord_state, fixture_payload("derived", "rate_representation"),
      value_reason = reasons,
      coordinate_provenance = x$coordinate_provenance
    ),
    "require a non-empty"
  )
})

test_that("value reasons cannot annotate finite or coordinate-absent cells", {
  x <- small_lt_matrix()
  reasons <- x$value_reason
  reasons[1, 1] <- "not_applicable"
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      value_reason = reasons,
      coordinate_provenance = x$coordinate_provenance
    ),
    "finite values"
  )

  reasons <- x$value_reason
  reasons[2, 1] <- "rate_trimmed"
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      value_reason = reasons,
      coordinate_provenance = x$coordinate_provenance
    ),
    "Coordinate-absent"
  )
})

test_that("non-observed coordinates cannot contain values", {
  x <- small_lt_matrix()
  bad_values <- x$values
  bad_values[2, 1] <- 1
  expect_error(
    lt_matrix(
      bad_values, x$coord_state, x$payload,
      value_reason = x$value_reason,
      coordinate_provenance = x$coordinate_provenance
    ),
    "non-observed coordinate state"
  )
})

test_that("payload distinguishes imported and derived rates", {
  x <- small_lt_matrix()
  imported_rate <- fixture_payload("imported", "rate_representation")
  imported_rate$spec_id <- "external_rate_spec/1.0.0"
  imported <- lt_matrix(
    x$values, x$coord_state, imported_rate,
    coordinate_provenance = x$coordinate_provenance
  )
  expect_identical(imported$payload$origin, "imported")

  derived_rate <- fixture_payload("derived", "rate_representation")
  derived <- lt_matrix(
    x$values, x$coord_state, derived_rate,
    coordinate_provenance = x$coordinate_provenance
  )
  expect_identical(derived$payload$origin, "derived")
  expect_false(identical(imported$payload$origin, derived$payload$origin))

  missing_payload_field <- x$payload[names(x$payload) != "units"]
  expect_error(
    lt_matrix(
      x$values, x$coord_state, missing_payload_field,
      coordinate_provenance = x$coordinate_provenance
    ),
    "payload.*missing"
  )

  bad_derived <- derived_rate
  bad_derived$spec_id <- NA_character_
  expect_error(
    lt_matrix(
      x$values, x$coord_state, bad_derived,
      coordinate_provenance = x$coordinate_provenance
    ),
    "derived payload requires"
  )
})

test_that("coordinate provenance supports unfrozen runtime objects", {
  x <- small_lt_matrix()
  incomplete <- x$coordinate_provenance
  incomplete$reference_tree_sha256 <- NULL
  runtime <- lt_matrix(
    x$values, x$coord_state, x$payload,
    coordinate_provenance = incomplete
  )
  expect_invisible(validate_lt_matrix(runtime))
  expect_true(is.na(runtime$coordinate_provenance$reference_tree_sha256))
  expect_identical(
    runtime$coordinate_provenance$audit_identity_status,
    "UNFROZEN_RUNTIME_OBJECT"
  )

  invalid <- x$coordinate_provenance
  invalid$ordered_branch_ledger_sha256 <- "not-a-hash"
  expect_error(
    lt_matrix(
      x$values, x$coord_state, x$payload,
      coordinate_provenance = invalid
    ),
    "SHA-256"
  )
})

test_that("NaN and infinite payload values are prohibited", {
  x <- small_lt_matrix()
  for (bad in c(NaN, Inf, -Inf)) {
    values <- x$values
    values[1, 1] <- bad
    expect_error(
      lt_matrix(
        values, x$coord_state, x$payload,
        coordinate_provenance = x$coordinate_provenance
      ),
      "NaN or infinite"
    )
  }
})
