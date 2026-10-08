test_that("terminal manifest retains all twenty-two required targets without self-certification", {
  names <- paste0("s", c(1, 2, 4, 5, 6), ".xlsx")
  targets <- LaTerra:::.lt_marine_target_artifacts(names)
  expect_identical(unique(targets$target_id), sprintf("TGT-%03d", 1:22))
  expect_true(all(targets$required))
  expect_true(all(targets$independent_acceptance == "NOT_CLAIMED"))
  expect_identical(unique(targets$source_root[targets$target_id == "TGT-021"]), "reporting")
  expect_true(all(c("nc_fingerprints.tsv", "focal_projections.tsv") %in%
    targets$path[targets$target_id == "TGT-019"]))
  expect_false(any(grepl("^TGT-02[34]$", targets$target_id)))
})

test_that("composite reporting run preserves scientific ledgers and distinct dependency roots", {
  tmp <- tempfile("lt_reporting_manifest_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  sci <- file.path(tmp, "science"); out <- file.path(tmp, "report")
  dir.create(sci); dir.create(out)
  sci <- normalizePath(sci); out <- normalizePath(out)
  names <- paste0("s", c(1, 2, 4, 5, 6), ".xlsx")
  targets <- LaTerra:::.lt_marine_target_artifacts(names)
  for (i in seq_len(nrow(targets))) {
    base <- if (targets$source_root[i] == "scientific") sci else out
    p <- file.path(base, targets$path[i]); dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
    if (!file.exists(p)) writeLines("synthetic target evidence only", p)
  }
  old <- lt_run("science-fixture", status = "completed")
  old$provenance$eligibility_masks <- list(M1 = list(directory = "M1", layers = "keyed fixture"))
  old$provenance$configuration_sha256 <- strrep("a", 64)
  saveRDS(old, file.path(sci, "run.rds"), version = 3)
  files <- list.files(sci, recursive = TRUE)
  hashes <- vapply(file.path(sci, files), LaTerra:::.lt_marine_sha, "")
  writeLines(paste(hashes, files, sep = "  "), file.path(sci, "SHA256SUMS"))
  reader <- list(run = old, manifest_sha256 = LaTerra:::.lt_marine_sha(file.path(sci, "SHA256SUMS")))
  a <- LaTerra:::.lt_marine_reporting_authority()
  reporting <- list(reporting_authority = a[c("authority_id", "authority_sha256", "source_sha256")],
    clarification = LaTerra:::.lt_marine_reporting_clarification(), reporting_package_version = "test-version",
    precision_clarification = LaTerra:::.lt_marine_reporting_precision_clarification(),
    reporting_installed_runtime_files = data.frame(path = "fixture", sha256 = strrep("b", 64)),
    reporting_environment = "test fixture", command = "terminal-reporting-fixture")
  before <- serialize(old, NULL, version = 3)
  x <- LaTerra:::.lt_marine_complete_reporting_run(reader, reporting, sci, out, names)
  expect_s3_class(x, "lt_run")
  expect_identical(serialize(old, NULL, version = 3), before)
  expect_identical(x$provenance$eligibility_masks$layers, old$provenance$eligibility_masks)
  expect_identical(x$results$dependency_roots, list(scientific = sci, reporting = out))
  expect_identical(x$results$scientific_source_configuration_sha256, old$provenance$configuration_sha256)
  expect_identical(x$provenance$configuration_sha256, LaTerra:::.lt_marine_sha(file.path(out, "reporting_recipe.rds")))
  expect_identical(x$results$target_status$target_id, sprintf("TGT-%03d", 1:22))
  expect_true(all(x$results$target_status$independent_acceptance == "NOT_CLAIMED"))
  expect_identical(x, readRDS(file.path(out, "run.rds")))
  expect_false(file.path(out, "run.rds") %in% x$provenance$output_hashes$path)
  expect_true(all(file.exists(x$provenance$output_hashes$path)))
  expect_true(LaTerra:::.lt_marine_reporting_clarification()$id %in% x$provenance$specification_ledger$specification_id)
  # A later edit to any linked source target must be caught, including targets
  # that do not supply an individual worksheet field.
  writeLines("changed scientific output", file.path(sci, "nc_fingerprints.tsv"))
  out2 <- file.path(tmp, "report2"); dir.create(out2)
  expect_error(LaTerra:::.lt_marine_complete_reporting_run(reader, reporting, sci, out2, names),
    class = "lt_error_marine_authority_mismatch")
  expect_false(file.exists(file.path(out2, "run.rds")))
})
