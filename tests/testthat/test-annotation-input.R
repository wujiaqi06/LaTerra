annotation_fixture <- function() {
  data.frame(gene_id = c("custom-gene-Z", "custom-gene-A"),
    candidate_module = c("specific label with spaces", "REVIEW_UNANNOTATED"),
    annotation_status = c("declared", "NA"),
    module_annotation_source = c("/does/not/exist/a.tsv", "historical text"),
    stringsAsFactors = FALSE)
}

test_that("keyed annotation parsing works with other genes and no snapshot hash", {
  f <- tempfile(fileext = ".tsv"); on.exit(unlink(f))
  x <- annotation_fixture()
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = TRUE)
  read <- LaTerra:::.lt_read_keyed_character_fields
  expect_identical(read(f, names(x), "gene_id"), x)
  x$gene_id <- rev(x$gene_id); x$candidate_module[1] <- "edited module"
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = TRUE)
  expect_identical(read(f, names(x), "gene_id"), x)
  expect_identical(read(f, c("gene_id", "annotation_status"), "gene_id"),
    x[c("gene_id", "annotation_status")])
})

test_that("structural annotation validation fails on ambiguous keys and columns", {
  f <- tempfile(); on.exit(unlink(f)); x <- annotation_fixture()
  read <- LaTerra:::.lt_read_keyed_character_fields
  put <- function(y) utils::write.table(y, f, sep = "\t", row.names = FALSE)
  put(x)
  expect_error(read(f, c("gene_id", "missing"), "gene_id"), "Missing annotation")
  expect_error(read(f, names(x), "absent"), "key")
  x$gene_id[2] <- x$gene_id[1]; put(x)
  expect_error(read(f, names(x), "gene_id"), "unique|duplicate")
  x$gene_id[2] <- ""; put(x)
  expect_error(read(f, names(x), "gene_id"), "empty|non-empty")
})

test_that("snapshot boundary is explicit and never opens historical source paths", {
  f <- tempfile(); on.exit(unlink(f)); x <- annotation_fixture()
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = TRUE)
  spec <- list(delivery_sha256 = digest::digest(file = f, algo = "sha256"),
    full_source_sha256 = paste(rep("a", 64), collapse = ""),
    authority_id = "test-declared-input", authority_sha256 = paste(rep("b", 64), collapse = ""))
  read <- LaTerra:::.lt_marine_turnover_annotation_input
  ans <- read(f, spec)
  expect_identical(ans$annotation, x)
  expect_true(ans$provenance$source_paths_are_inert_strings)
  expect_false(ans$provenance$original_full_source_required_at_runtime)
  expect_false(ans$provenance$feature_selection_authority)
  expect_identical(ans$provenance$full_source_sha256, spec$full_source_sha256)
  expect_error(read(f, list()), class = "lt_error_marine_annotation_authority")
  bad <- spec; bad$full_source_sha256 <- "not-a-hash"
  expect_error(read(f, bad), class = "lt_error_marine_annotation_authority")
  bad <- spec; bad$delivery_sha256 <- paste(rep("c", 64), collapse = "")
  expect_error(read(f, bad), class = "lt_error_marine_authority_mismatch")
  x$coefficient <- c(1, 2)
  utils::write.table(x, f, sep = "\t", row.names = FALSE, quote = TRUE)
  spec$delivery_sha256 <- digest::digest(file = f, algo = "sha256")
  expect_error(read(f, spec), class = "lt_error_marine_annotation_fields")
})
