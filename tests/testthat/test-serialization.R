test_that("S3 objects round-trip without changing ledgers or masks", {
  x <- small_lt_matrix()
  path <- tempfile(fileext = ".rds")

  lt_write_object(x, path)
  restored <- lt_read_object(path)

  expect_s3_class(restored, "lt_matrix")
  expect_identical(restored, x)
})

test_that("lt_run provides the complete provenance container", {
  run <- lt_run("draft_run", matrix = small_lt_matrix())

  expect_s3_class(run, "lt_run")
  expect_identical(names(run$provenance), names(lt_empty_provenance()))
  expect_identical(names(lt_manifest_schema()), names(lt_empty_provenance()))
  expect_identical(
    names(run$provenance$fold_ledger),
    c("fold_id", "taxon_id", "role", "group")
  )
  expect_identical(
    names(run$provenance$seed_ledger),
    c("scope", "seed", "rng_kind", "normal_kind", "sample_kind", "r_version")
  )
  expect_true(all(c("specification_ledger", "filter_ledger") %in% names(run$provenance)))
  expect_false("scientific_certification" %in% names(run$provenance))

  path <- tempfile(fileext = ".rds")
  lt_write_object(run, path)
  expect_identical(lt_read_object(path), run)
})

test_that("provenance validator enforces taxon-level fold structure", {
  provenance <- lt_empty_provenance()
  provenance$fold_ledger <- data.frame(
    fold_id = "fold_1",
    taxon_id = c("Taxon_a", "Taxon_b"),
    role = c("test", "train"),
    group = c("Genus_a", "Genus_b"),
    stringsAsFactors = FALSE
  )
  expect_invisible(validate_lt_provenance(provenance))

  provenance$fold_ledger$role[1] <- "validation"
  expect_error(validate_lt_provenance(provenance), "train, test, or excluded")
})

test_that("provenance validator rejects obsolete fold ledgers", {
  provenance <- lt_empty_provenance()
  provenance$fold_ledger <- data.frame(
    fold_id = character(), held_out_group = character(),
    stringsAsFactors = FALSE
  )
  expect_error(validate_lt_provenance(provenance), "fold_ledger")
})
