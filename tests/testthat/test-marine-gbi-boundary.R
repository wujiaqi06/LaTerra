test_that("historical computation view and read-back preserve original axes and raw values", {
  raw <- list(values = rbind(g1 = c(0, 1.2345678901234567, 2.5, 9.7),
    g2 = c(1.1, NA, 2.1, 6), g3 = c(0.7, 3.1, 1.8, 4.2)))
  colnames(raw$values) <- c("oldA", "oldB", "oldC", "oldD")
  raw$gene_ids <- rownames(raw$values); raw$branch_ids <- colnames(raw$values)
  mapping <- data.frame(new_branch_label = paste0("native", 1:4),
    old_branch_label = c("oldC", "oldA", "oldD", "oldB"))
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive = TRUE))
  before <- serialize(raw, NULL, version = 3)
  b <- LaTerra:::.lt_marine_gbi_history_kernel(raw, mapping, out)
  expect_identical(serialize(raw, NULL, version = 3), before)
  expect_identical(dimnames(b$gbi), dimnames(raw$values))
  expect_identical(dimnames(b$trim_mask), dimnames(raw$values))
  expect_identical(dimnames(b$trimmed_values), dimnames(raw$values))
  expect_identical(b$branch_effect$branch, raw$branch_ids)
  expect_identical(b$gbi_boundary$native_from_public, c(3L, 1L, 4L, 2L))
  expect_identical(b$gbi_boundary$public_from_native, c(2L, 4L, 1L, 3L))
  expect_identical(b$gbi_boundary$native_from_public[b$gbi_boundary$public_from_native], 1:4)
  expect_false(b$gbi_boundary$source_admission_claimed)
  expect_false(b$gbi_boundary$intermediate_full_payload_duplicated)
  read_back <- LaTerra:::.lt_marine_read_generated_gbi(file.path(out, b$gbi_boundary$consumed_file), raw$gene_ids, raw$branch_ids)
  expect_identical(b$gbi, read_back)
  expect_identical(b$gbi_boundary$consumed_values_sha256, LaTerra:::.lt_marine_payload_sha(read_back))
  expect_error(LaTerra:::.lt_marine_consumed_gbi(b, out), class = "lt_error_marine_gbi_dependency")
  expect_error(LaTerra:::.lt_marine_gbi_history_kernel(raw, mapping, out), class = "lt_error_marine_output_exists")
  # The generic operator still evaluates the user-supplied order in memory.
  generic <- LaTerra:::.lt_marine_baseline(raw)
  expect_null(generic$gbi_boundary)
  manual <- raw$values
  for (i in 1:3) manual[i, which(manual[i, ] >= quantile(manual[i, ], .975, na.rm=TRUE))] <- NA_real_
  avg <- function(x) if (all(is.na(x))) NA_real_ else mean(x[!is.na(x)])
  expected <- sweep(sweep(manual, 2L, apply(manual, 2L, avg), "/"), 1L, apply(manual, 1L, avg), "/")
  expected[!is.finite(expected)] <- NA_real_
  expect_identical(generic$gbi, expected)
})

test_that("literal text permutation neither parses nor rounds numeric tokens", {
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive = TRUE))
  source <- file.path(out, "native.tsv"); dest <- file.path(out, "old.tsv")
  writeLines(c("gene\tN1\tN2\tN3", "g1\t1.2345678901234567\t0\tNA",
    "g2\t1e-120\t3.1415926535897931\t2.0000000000000001"), source)
  order <- data.frame(native_branch_id = c("N1","N2","N3"), public_branch_id = c("B3","B1","B2"))
  LaTerra:::.lt_marine_reorder_gbi_tokens(source, dest, order, c("B1","B2","B3"))
  expect_identical(readLines(dest), c("gene\tB1\tB2\tB3", "g1\t0\tNA\t1.2345678901234567",
    "g2\t3.1415926535897931\t2.0000000000000001\t1e-120"))
  expect_error(LaTerra:::.lt_marine_reorder_gbi_tokens(source, dest, order, c("B1","B2","B3")), class="lt_error_marine_output_exists")
  expect_error(LaTerra:::.lt_marine_reorder_gbi_tokens(source, file.path(out,"wrong.tsv"), order, c("B1","B2","missing")), class="lt_error_marine_gbi_order")
  v <- LaTerra:::.lt_marine_read_generated_gbi(dest, c("g1","g2"), c("B1","B2","B3"))
  expect_identical(v["g1","B1"], 0)
  expect_true(is.na(v["g1","B2"]))
  expect_identical(v["g1","B3"], as.numeric("1.2345678901234567"))
  expect_error(LaTerra:::.lt_marine_read_generated_gbi(dest,c("g2","g1"),c("B1","B2","B3")),class="lt_error_marine_gbi_reader")
  bad <- file.path(out,"bad.tsv"); writeLines(c("gene\tB1", "g1\tsurprise"),bad)
  expect_error(LaTerra:::.lt_marine_read_generated_gbi(bad,"g1","B1"),class="lt_error_marine_gbi_reader")
})

test_that("canonical Marine intake cannot be activated by a caller-rehashed source", {
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive = TRUE))
  raw <- file.path(out,"raw.tsv"); map <- file.path(out,"map.tsv")
  writeLines(c("gene\tB1", "g1\t1"),raw)
  writeLines(c("new_branch_label\told_branch_label", "B1\tB1"),map)
  fake <- list(matrix_path=raw,crosswalk_path=map,
    raw_sha256=digest::digest(file=raw,algo="sha256"),crosswalk_sha256=digest::digest(file=map,algo="sha256"))
  expect_error(LaTerra:::.lt_marine_validate_gbi_intake(fake),class="lt_error_marine_authority_mismatch")
  expect_error(LaTerra:::.lt_marine_gbi_intake(raw,map),class="lt_error_marine_authority_mismatch")
  expect_false(any(c("profile","hash","raw","baseline","fit","authorization") %in% names(formals(lt_run_marine))))
  expect_false(any(c("profile","hash","gbi","baseline","fit","authorization") %in% names(formals(lt_run_marine_m2))))
  a <- LaTerra:::.lt_marine_gbi_authority()
  expect_identical(digest::digest(file=system.file(a$producer_installed_path,package="LaTerra"),algo="sha256"),a$producer_sha256)
  expect_identical(lt_marine_profile()$execution_path_addendum$authority_sha256,a$sha256)
  expect_false(lt_marine_profile()$execution_path_addendum$generic_S2_affected)
  expect_false(any(grepl("expected|archived_GBI",names(a))))
})
