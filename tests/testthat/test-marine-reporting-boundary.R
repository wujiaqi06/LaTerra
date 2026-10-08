reporting_fixture <- function(source_id = "OUT_TABLE_S4", sheet = "Sheet1") {
  x <- LaTerra:::.lt_marine_reporting_layout(source_id, sheet)
  dynamic <- c("EXISTING_ADMITTED_INPUT_JOIN", "RECOMPUTE_FROM_THIS_RUN", "THIS_RUN_KEY_OR_MANIFEST_JOIN")
  mask <- matrix(x$roles %in% dynamic, nrow(x$roles))
  rows <- which(rowSums(mask) > 0L)
  headers <- unlist(x$schema$headers, use.names = FALSE)
  cols <- headers[colSums(mask) > 0L]
  keys <- unlist(x$schema$key_fields, use.names = FALSE)
  current <- as.data.frame(setNames(rep(list(seq_along(rows) + 1000), length(cols)), cols),
    check.names = FALSE, stringsAsFactors = FALSE)
  for (j in seq_along(keys)) current[[keys[j]]] <- vapply(rows, function(r)
    as.character(x$schema$ordered_keys[[r - 1L]][[j]]), "")
  current
}

test_that("trusted reporting projection covers exactly the original eighteen schemas", {
  a <- LaTerra:::.lt_marine_reporting_authority()
  expect_length(a$schemas, 18L)
  expect_length(a$references, 8169L)
  expect_equal(nrow(a$fields), 496L)
  total <- 0L
  for (schema in a$schemas) {
    x <- LaTerra:::.lt_marine_reporting_layout(schema$source_id, schema$sheet)
    expect_false(anyNA(x$roles))
    expect_identical(unname(x$values[1, ]), schema$headers)
    total <- total + length(x$roles)
  }
  expect_equal(total, 53108L)
  expect_true(all(vapply(a$references, function(x) x$role %in%
    c("FROZEN_REPORTING_REFERENCE", "HISTORICAL_REFERENCE_IMPORTED_NOT_RECOMPUTED"), TRUE)))
  expect_error(LaTerra:::.lt_marine_reporting_layout("OUT_TABLE_S4", "renamed"),
    class = "lt_error_marine_reporting_schema")
})

test_that("current key coverage precedes historical ordering and never selects rows", {
  current <- reporting_fixture()
  expect_equal(nrow(current), 22L)
  assemble <- function(d) LaTerra:::.lt_marine_assemble_reporting_sheet("OUT_TABLE_S4", "Sheet1", d)
  x <- assemble(current)
  reversed <- assemble(current[22:1, , drop = FALSE])
  expect_identical(x$values, reversed$values)
  expect_true(x$provenance$key_coverage_before_archived_order)
  expect_false(x$provenance$scientific_or_target_pass_claimed)
  expect_error(assemble(current[-1, ]), class = "lt_error_marine_reporting_keys")
  expect_error(assemble(rbind(current, current[1, ])), class = "lt_error_marine_reporting_keys")
  bad <- current; bad$metric[1] <- "invented_metric"
  expect_error(assemble(bad), class = "lt_error_marine_reporting_keys")
  bad <- current; bad$metric[1] <- NA_character_
  expect_error(assemble(bad), class = "lt_error_marine_reporting_keys")
  bad <- current; bad$metric[1] <- " "
  expect_error(assemble(bad), class = "lt_error_marine_reporting_keys")
  bad <- current; bad$metric[1] <- paste0(bad$metric[1], "\r")
  expect_error(assemble(bad), class = "lt_error_marine_reporting_keys")
})

test_that("computed fields cannot be imported and explicit blanks differ from zero", {
  current <- reporting_fixture()
  assemble <- function(d) LaTerra:::.lt_marine_assemble_reporting_sheet("OUT_TABLE_S4", "Sheet1", d)
  x <- assemble(current)
  expect_equal(x$values[[2, 3]], 1001)
  changed <- current; changed$value[1] <- 0
  expect_identical(assemble(changed)$values[[2, 3]], 0)
  changed$value[1] <- NA_real_
  expect_null(assemble(changed)$values[[2, 3]])
  changed$value[1] <- Inf
  expect_error(assemble(changed), class = "lt_error_marine_reporting_payload")
  changed <- current; changed$notes <- "caller-authored reference"
  expect_error(assemble(changed), class = "lt_error_marine_reporting_fields")
  expect_error(assemble(current[, names(current) != "value"]), class = "lt_error_marine_reporting_fields")
  expect_identical(x$provenance$historical_reference_qualifier,
    "historical_reference_imported_not_recomputed")
  expect_equal(sum(x$roles == "HISTORICAL_REFERENCE_IMPORTED_NOT_RECOMPUTED"), 48L)
  expect_false(x$provenance$original_workbook_read)
})

test_that("reference-only sheets accept no caller result payload", {
  for (pair in list(c("OUT_TABLE_S5", "column_description"),
    c("OUT_TABLE_S6", "overview"), c("OUT_TABLE_S6", "column_description"))) {
    x <- LaTerra:::.lt_marine_assemble_reporting_sheet(pair[1], pair[2])
    expect_true(all(x$assigned))
    expect_error(LaTerra:::.lt_marine_assemble_reporting_sheet(pair[1], pair[2], data.frame(x = 1)),
      class = "lt_error_marine_reporting_fields")
  }
})

test_that("coordinated projection or role rehash cannot expand reporting authority", {
  asset <- system.file("reporting", "marine_TGT021_projection_v1.rds", package = "LaTerra")
  d <- readRDS(asset)
  temp <- tempfile(fileext = ".rds"); on.exit(unlink(temp))
  d$references[[1]]$value <- "edited reference"
  d$source_sha256[] <- "caller rehash"
  saveRDS(d, temp, version = 3, compress = "xz")
  expect_error(LaTerra:::.lt_marine_reporting_authority(temp), class = "lt_error_marine_authority_mismatch")
  d <- readRDS(asset); d$fields$role[d$fields$role == "RECOMPUTE_FROM_THIS_RUN"] <- "FROZEN_REPORTING_REFERENCE"
  d$authority_sha256 <- digest::digest(d, algo = "sha256")
  saveRDS(d, temp, version = 3, compress = "xz")
  expect_error(LaTerra:::.lt_marine_reporting_authority(temp), class = "lt_error_marine_authority_mismatch")
  expect_identical(names(formals(LaTerra:::.lt_marine_assemble_reporting_sheet)), c("source_id", "sheet", "current"))
})

test_that("row display ordinals follow verified models and never define their keys", {
  current <- reporting_fixture("OUT_TABLE_S6", "Fig5A_fold_summary")
  x <- LaTerra:::.lt_marine_assemble_reporting_sheet("OUT_TABLE_S6", "Fig5A_fold_summary", current[8:1, ])
  expect_identical(unlist(x$values[2:9, 1], use.names = FALSE), 1:8)
  expect_error(LaTerra:::.lt_marine_assemble_reporting_sheet("OUT_TABLE_S6", "Fig5A_fold_summary", current[-1, ]),
    class = "lt_error_marine_reporting_keys")
})
