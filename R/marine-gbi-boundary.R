# ADDENDUM006: historical execution order/serialization for Marine replay only.
# Generic S2 operators continue to use their supplied order and in-memory values.

.lt_marine_gbi_authority <- function() {
  list(id = "TARGET001-AUTHORITY-ADDENDUM006-GBI-ORDER-SERIALIZATION001",
    sha256 = "cebaafc6fdfb019d36d9fcb0960f171a765f65cb6a767f06c9a156a5d9ec99c5",
    raw_sha256 = "3ec73fa19d2706d6ff86f397e1235f1303f6364f7cfb422189d0bd3917551a7e",
    # Canonical input payload identity, not an output oracle or caller assertion.
    raw_values_sha256 = "6cbc499299e94601f96f9e7f17c849e5bb41e6f316c79687c9d6a64b60c46a9d",
    crosswalk_sha256 = "6737625cd8e6dbeb7a1f3d8b247ef71d9efe836ecfa4d124c25773e0662ca14c",
    producer_sha256 = "1655a0f5dd84a66d7b80d4c3dddc10c483ce319dac02713075e414b8a27aa316",
    producer_installed_path = "authority/s2_gbi_source.R",
    mapper_sha256 = "486f2e9ee85aae7bf62e57fcbcfe5ecf74db81e1a1865fe1e3b6908c80149a16",
    fingerprint_source_sha256 = "f1893550f1f395c8b4a4cf194e06fbc4220c3aadeef4f7f5b5b9a40a704bdc47",
    full_data_source_sha256 = "debc1792f941ef7b34b587acd5f590b84471666e95730bf7bf90dd7f84a6d53e")
}

.lt_marine_payload_sha <- function(x) {
  digest::digest(serialize(x, NULL, version = 3), algo = "sha256", serialize = FALSE)
}

.lt_marine_gbi_intake <- function(matrix_path, crosswalk_path) {
  authority <- .lt_marine_gbi_authority()
  .lt_marine_verify(matrix_path, authority$raw_sha256, "canonical_original_raw")
  .lt_marine_verify(crosswalk_path, authority$crosswalk_sha256, "canonical_computation_order_crosswalk")
  source <- system.file(authority$producer_installed_path, package = "LaTerra")
  .lt_marine_verify(source, authority$producer_sha256, "canonical_inert_GBI_source")
  x <- list(raw = .lt_marine_parse_matrix(matrix_path),
    mapping = .lt_marine_read_tsv(crosswalk_path),
    matrix_path = matrix_path, crosswalk_path = crosswalk_path)
  .lt_marine_validate_gbi_intake(x)
  x
}

.lt_marine_validate_gbi_intake <- function(x) {
  # Rebind to exact compiled authority, not a mutable context's reported hashes.
  a <- .lt_marine_gbi_authority()
  .lt_marine_verify(x$matrix_path, a$raw_sha256, "canonical_original_raw")
  .lt_marine_verify(x$crosswalk_path, a$crosswalk_sha256, "canonical_computation_order_crosswalk")
  .lt_marine_verify(system.file(a$producer_installed_path, package = "LaTerra"),
    a$producer_sha256, "canonical_inert_GBI_source")
  if (!identical(.lt_marine_payload_sha(x$raw$values), a$raw_values_sha256) ||
      !identical(x$raw$gene_ids, rownames(x$raw$values)) ||
      !identical(x$raw$branch_ids, colnames(x$raw$values))) {
    .lt_marine_abort("Marine raw payload or ordered axes differ from the canonical input.", "gbi_source")
  }
  canonical <- .lt_marine_read_tsv(x$crosswalk_path)
  if (!identical(x$mapping, canonical) ||
      !identical(canonical$new_branch_label, paste0("B", seq_len(601L))) ||
      anyDuplicated(canonical$old_branch_label) ||
      !setequal(canonical$old_branch_label, x$raw$branch_ids)) {
    .lt_marine_abort("Marine computation order must be the complete canonical inverse keyed view.", "gbi_order")
  }
  invisible(TRUE)
}

.lt_marine_gbi_order <- function(raw, mapping) {
  .lt_assert_unique_ids(raw$gene_ids, "raw gene IDs")
  .lt_assert_unique_ids(raw$branch_ids, "public branch IDs")
  .lt_assert_unique_ids(mapping$new_branch_label, "native branch labels")
  .lt_assert_unique_ids(mapping$old_branch_label, "mapped public branch IDs")
  if (!identical(dimnames(raw$values), list(raw$gene_ids, raw$branch_ids)) ||
      !setequal(raw$branch_ids, mapping$old_branch_label)) {
    .lt_marine_abort("Computational view must bijectively preserve the supplied scientific axes.", "gbi_order")
  }
  data.frame(native_position = seq_along(mapping$new_branch_label),
    native_branch_id = mapping$new_branch_label,
    public_position = match(mapping$old_branch_label, raw$branch_ids),
    public_branch_id = mapping$old_branch_label, stringsAsFactors = FALSE)
}

.lt_marine_reorder_gbi_tokens <- function(src, dst, order, public_ids) {
  if (file.exists(dst)) .lt_marine_abort("Generated GBI destination exists; no overwrite.", "output_exists")
  inp <- file(src, "rt"); on.exit(close(inp), add = TRUE)
  dest <- file(dst, "wt"); on.exit(close(dest), add = TRUE)
  header <- strsplit(readLines(inp, n = 1L), "\t", fixed = TRUE)[[1L]]
  if (!identical(header, c("gene", order$native_branch_id))) {
    .lt_marine_abort("Native GBI header does not match its recorded computation order.", "gbi_order")
  }
  .lt_assert_unique_ids(public_ids, "restored public branch IDs")
  idx <- match(order$native_branch_id[match(public_ids, order$public_branch_id)], header)
  if (length(idx) != nrow(order) || anyNA(idx) || anyDuplicated(idx)) {
    .lt_marine_abort("Restored GBI order is not a complete inverse permutation.", "gbi_order")
  }
  writeLines(paste(c(header[1L], public_ids), collapse = "\t"), dest)
  repeat {
    lines <- readLines(inp, n = 256L)
    if (!length(lines)) break
    tokens <- strsplit(lines, "\t", fixed = TRUE)
    if (any(lengths(tokens) != length(header))) {
      .lt_marine_abort("Malformed generated GBI token row; no partial reader fallback.", "gbi_reader")
    }
    writeLines(vapply(tokens, function(v) paste(c(v[1L], v[idx]), collapse = "\t"), ""), dest)
  }
  invisible(dst)
}

.lt_marine_read_generated_gbi <- function(path, genes, branches) {
  # Literal original downstream reader. No precision selection or old output read.
  dat <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!identical(names(dat), c("gene", branches)) || !identical(dat$gene, genes)) {
    .lt_marine_abort("Generated/read-back GBI axes differ from the actual input.", "gbi_reader")
  }
  rownames(dat) <- dat$gene; dat$gene <- NULL
  valid <- vapply(dat, function(x) is.numeric(x) || (is.logical(x) && all(is.na(x))), TRUE)
  if (!all(valid)) .lt_marine_abort("Generated GBI has an unexpected nonnumeric token.", "gbi_reader")
  value <- as.matrix(dat)
  storage.mode(value) <- "numeric"
  if (any(is.nan(value)) || any(!is.na(value) & !is.finite(value))) {
    .lt_marine_abort("Generated GBI contains a nonfinite payload.", "gbi_reader")
  }
  value
}

.lt_marine_gbi_history_kernel <- function(raw, mapping, output_dir) {
  # Low-level testable transformation; does not assert canonical source admission.
  order <- .lt_marine_gbi_order(raw, mapping)
  native_from_public <- order$public_position
  public_from_native <- match(raw$branch_ids, order$public_branch_id)
  native <- list(values = raw$values[, native_from_public, drop = FALSE],
    gene_ids = raw$gene_ids, branch_ids = order$native_branch_id)
  colnames(native$values) <- native$branch_ids
  baseline <- .lt_marine_baseline(native)
  native_path <- file.path(output_dir, "GBI.generated_native.tsv")
  consumed_path <- file.path(output_dir, "GBI.generated_consumed.oldlabels.tsv")
  if (any(file.exists(c(native_path, consumed_path)))) {
    .lt_marine_abort("GBI boundary outputs already exist; no overwrite.", "output_exists")
  }
  intermediate_native_sha <- .lt_marine_payload_sha(baseline$gbi)
  # Same call/defaults as original producer; no digits or independent rounding.
  utils::write.table(data.frame(gene = raw$gene_ids, baseline$gbi, check.names = FALSE),
    file = native_path, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
  .lt_marine_reorder_gbi_tokens(native_path, consumed_path, order, raw$branch_ids)
  restore <- function(m) {
    m <- m[, public_from_native, drop = FALSE]
    colnames(m) <- raw$branch_ids
    m
  }
  intermediate_old_sha <- .lt_marine_payload_sha(restore(baseline$gbi))
  baseline$gbi <- .lt_marine_read_generated_gbi(consumed_path, raw$gene_ids, raw$branch_ids)
  baseline$trimmed_values <- restore(baseline$trimmed_values)
  baseline$trim_mask <- restore(baseline$trim_mask)
  baseline$branch_effect <- baseline$branch_effect[public_from_native, , drop = FALSE]
  baseline$branch_effect$branch <- raw$branch_ids
  rownames(baseline$branch_effect) <- NULL
  baseline$gbi_boundary <- list(
    consumption = "this_run_generated_literal_writer_token_reorder_reader",
    intermediate_native_values_sha256 = intermediate_native_sha,
    intermediate_restored_values_sha256 = intermediate_old_sha,
    consumed_values_sha256 = .lt_marine_payload_sha(baseline$gbi),
    native_file = basename(native_path), native_file_sha256 = .lt_marine_sha(native_path),
    consumed_file = basename(consumed_path), consumed_file_sha256 = .lt_marine_sha(consumed_path),
    computational_order = order,
    native_from_public = native_from_public, public_from_native = public_from_native,
    order_sha256 = .lt_marine_payload_sha(order),
    gene_axes_sha256 = .lt_marine_payload_sha(raw$gene_ids),
    public_axes_sha256 = .lt_marine_payload_sha(raw$branch_ids),
    writer = "original_CODE_GBI_write.table_quote_FALSE_sep_tab_row.names_FALSE_na_NA",
    reorder = "written_text_tokens_only_no_numeric_parser",
    reader = "original_read.delim_check.names_FALSE_gene_rownames_as.matrix_storage.mode_numeric",
    intermediate_full_payload_duplicated = FALSE, source_admission_claimed = FALSE)
  baseline
}

.lt_marine_historical_baseline <- function(intake, output_dir) {
  .lt_marine_validate_gbi_intake(intake)
  baseline <- .lt_marine_gbi_history_kernel(intake$raw, intake$mapping, output_dir)
  baseline$gbi_boundary$source_authority <- .lt_marine_gbi_authority()
  baseline$gbi_boundary$source_admission_claimed <- TRUE
  saveRDS(baseline$gbi_boundary, file.path(output_dir, "GBI_consumption_boundary.rds"), version = 3)
  .lt_marine_write_tsv(baseline$gbi_boundary$computational_order,
    file.path(output_dir, "GBI_computational_order.tsv"))
  baseline
}

.lt_marine_consumed_gbi <- function(baseline, output_dir) {
  b <- baseline$gbi_boundary
  if (!is.list(b) || !isTRUE(b$source_admission_claimed) ||
      !identical(b$source_authority, .lt_marine_gbi_authority()) ||
      !identical(b$consumption, "this_run_generated_literal_writer_token_reorder_reader") ||
      !identical(b$consumed_file, "GBI.generated_consumed.oldlabels.tsv")) {
    .lt_marine_abort("This Marine dependency lacks the admitted generated/read-back GBI boundary; recompute from raw.", "gbi_dependency")
  }
  file <- file.path(output_dir, b$consumed_file)
  .lt_marine_verify(file, b$consumed_file_sha256, "this_run_generated_GBI_not_archived_input")
  value <- .lt_marine_read_generated_gbi(file, rownames(baseline$gbi), colnames(baseline$gbi))
  if (!identical(.lt_marine_payload_sha(value), b$consumed_values_sha256) ||
      !identical(value, baseline$gbi)) {
    .lt_marine_abort("GBI cache and generated consumed bits differ; recompute stale dependencies.", "gbi_dependency")
  }
  value
}
