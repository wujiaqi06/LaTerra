# Terminal reporting boundary for the explicit historical Marine profile only.
# Ordinary data objects never pass through this snapshot-identity requirement.
# No original workbook, fit, code evidence, or all-values JSON is a dependency.
.lt_marine_reporting_authority <- function(path = system.file("reporting",
    "marine_TGT021_projection_v1.rds", package = "LaTerra")) {
  .lt_marine_verify(path, "c9d5927108dd2653129789364544a0b3a44f172f017cc8b3a0fe165b20c31bd3",
    "trusted_TGT021_reporting_projection_not_caller_authority")
  readRDS(path)
}

.lt_marine_report_cell <- function(cell) {
  if (length(cell) != 1L || is.na(cell) || !grepl("^[A-Z]+[1-9][0-9]*$", cell))
    .lt_marine_abort("Invalid reporting cell address.", "reporting_schema")
  letters <- utf8ToInt(sub("[0-9]+$", "", cell)) - 64L
  column <- Reduce(function(a, b) 26L * a + b, letters, init = 0L)
  c(row = as.integer(sub("^[A-Z]+", "", cell)), col = as.integer(column))
}

.lt_marine_report_range <- function(address) {
  parts <- strsplit(address, ":", fixed = TRUE)[[1L]]
  if (!length(parts) %in% 1:2)
    .lt_marine_abort("Invalid reporting range.", "reporting_schema")
  a <- .lt_marine_report_cell(parts[1L]); b <- .lt_marine_report_cell(utils::tail(parts, 1L))
  if (any(b < a)) .lt_marine_abort("Reversed reporting range.", "reporting_schema")
  list(rows = seq.int(a[[1L]], b[[1L]]), cols = seq.int(a[[2L]], b[[2L]]))
}

.lt_marine_report_keys <- function(rows) {
  vapply(rows, function(row) {
    vals <- vapply(row, function(x) {
      if (length(x) != 1L || is.na(x) || !is.atomic(x))
        .lt_marine_abort("Missing or ambiguous reporting key.", "reporting_keys")
      v <- as.character(x)
      if (!nzchar(trimws(v)) || grepl("\r", v, fixed = TRUE))
        .lt_marine_abort("Blank or ambiguous reporting key.", "reporting_keys")
      v
    }, "")
    paste(vals, collapse = "\r")
  }, "")
}

.lt_marine_reporting_layout <- function(source_id, sheet) {
  authority <- .lt_marine_reporting_authority()
  matches <- vapply(authority$schemas, function(x)
    identical(x$source_id, source_id) && identical(x$sheet, sheet), TRUE)
  if (sum(matches) != 1L)
    .lt_marine_abort("Unknown or ambiguous historical reporting sheet.", "reporting_schema")
  schema <- authority$schemas[[which(matches)]]
  box <- .lt_marine_report_range(schema$range)
  nr <- max(box$rows); nc <- max(box$cols)
  if (min(box$rows) != 1L || min(box$cols) != 1L || length(schema$headers) != nc ||
      length(schema$ordered_keys) != nr - 1L || anyDuplicated(unlist(schema$headers)))
    .lt_marine_abort("Malformed reporting schema rectangle or headers.", "reporting_schema")
  values <- matrix(vector("list", nr * nc), nr, nc)
  roles <- matrix(NA_character_, nr, nc)
  assigned <- matrix(FALSE, nr, nc)
  fields <- authority$fields[authority$fields$source_id == source_id &
    authority$fields$sheet == sheet, , drop = FALSE]
  for (i in seq_len(nrow(fields))) {
    r <- .lt_marine_report_range(fields$range[i])
    if (fields$source_sha256[i] != schema$source_sha256 || length(r$cols) != 1L ||
        max(r$rows) > nr || max(r$cols) > nc ||
        fields$field[i] != schema$headers[[r$cols]] || any(!is.na(roles[r$rows, r$cols])))
      .lt_marine_abort("Reporting field/role/source binding conflict.", "reporting_schema")
    roles[r$rows, r$cols] <- fields$role[i]
  }
  permitted <- c("EXACT_HEADER_SCHEMA", "EXISTING_ADMITTED_INPUT_JOIN",
    "FROZEN_REPORTING_REFERENCE", "RECOMPUTE_FROM_THIS_RUN", "THIS_RUN_KEY_OR_MANIFEST_JOIN",
    "SCHEMA_ORDER_METADATA", "HISTORICAL_REFERENCE_IMPORTED_NOT_RECOMPUTED")
  if (anyNA(roles) || any(!roles %in% permitted) || any(roles[1L, ] != "EXACT_HEADER_SCHEMA"))
    .lt_marine_abort("Unclassified or invalid reporting cells.", "reporting_schema")
  values[1L, ] <- schema$headers; assigned[1L, ] <- TRUE
  refs <- authority$references[vapply(authority$references, function(x)
    identical(x$source_id, source_id) && identical(x$sheet, sheet), TRUE)]
  for (ref in refs) {
    p <- .lt_marine_report_cell(ref$cell); r <- p[[1L]]; c <- p[[2L]]
    if (r > nr || c > nc || ref$source_sha256 != schema$source_sha256 ||
        !ref$role %in% c("FROZEN_REPORTING_REFERENCE", "HISTORICAL_REFERENCE_IMPORTED_NOT_RECOMPUTED") ||
        roles[r, c] != ref$role || assigned[r, c])
      .lt_marine_abort("Unclassified, duplicate, or source-mismatched reference cell.", "reporting_schema")
    values[r, c] <- list(ref$value); assigned[r, c] <- TRUE
  }
  # This layout is disposable; entry points reload the trusted projection rather
  # than accepting an edited layout or its self-supplied fingerprint as authority.
  list(schema = schema, values = values, roles = roles, assigned = assigned,
    authority_id = authority$authority_id, authority_sha256 = authority$authority_sha256,
    source_sha256 = authority$source_sha256)
}

.lt_marine_assemble_reporting_sheet <- function(source_id, sheet, current = NULL) {
  x <- .lt_marine_reporting_layout(source_id, sheet)
  dynamic <- c("EXISTING_ADMITTED_INPUT_JOIN", "RECOMPUTE_FROM_THIS_RUN", "THIS_RUN_KEY_OR_MANIFEST_JOIN")
  needed <- which(rowSums(matrix(x$roles %in% dynamic, nrow(x$roles))) > 0L)
  headers <- unlist(x$schema$headers, use.names = FALSE)
  keys <- unlist(x$schema$key_fields, use.names = FALSE)
  key_cols <- match(keys, headers)
  expected_keys <- .lt_marine_report_keys(x$schema$ordered_keys)
  if (anyNA(key_cols) || anyDuplicated(expected_keys))
    .lt_marine_abort("Ambiguous frozen reporting keys.", "reporting_keys")
  if (length(needed)) {
    required_cols <- headers[which(colSums(matrix(x$roles %in% dynamic, nrow(x$roles))) > 0L)]
    if (!is.data.frame(current) || anyDuplicated(names(current)) ||
        !setequal(names(current), union(required_cols, keys)))
      .lt_marine_abort("Supply exactly the required this-run fields and scientific keys; no imported result columns.",
        "reporting_fields")
    actual_keys <- .lt_marine_report_keys(lapply(seq_len(nrow(current)), function(i)
      as.list(current[i, keys, drop = FALSE])))
    required_keys <- expected_keys[needed - 1L]
    if (anyDuplicated(actual_keys) || !setequal(actual_keys, required_keys))
      .lt_marine_abort("This-run reporting keys are missing, extra, duplicated, or ambiguous; archived order cannot select them.",
        "reporting_keys")
    order <- match(required_keys, actual_keys)
    for (i in seq_along(needed)) {
      r <- needed[i]
      for (c in which(x$roles[r, ] %in% dynamic)) {
        v <- current[[headers[c]]][order[i]]
        if (is.list(v)) v <- v[[1L]]
        if (length(v) > 1L || (!is.null(v) && !is.atomic(v)) ||
            (is.numeric(v) && any(is.infinite(v))))
          .lt_marine_abort("Reporting payload must be a finite scalar or an explicit unavailable/blank cell.",
            "reporting_payload")
        # Original TSV -> XLSX assembly writes empty strings/missing scalars as
        # empty worksheet cells. Preserve whitespace and real numeric zero.
        if (length(v) == 1L && (is.na(v) || identical(v, ""))) v <- NULL
        x$values[r, c] <- list(v); x$assigned[r, c] <- TRUE
      }
    }
  } else if (!is.null(current)) {
    .lt_marine_abort("A reference-only sheet cannot accept caller result values.", "reporting_fields")
  }
  # Schema metadata is applied only AFTER exact fresh key coverage. It never
  # defines a gene set, a model, a fold, or a random world.
  for (r in seq.int(2L, nrow(x$roles))) {
    for (c in which(x$roles[r, ] == "SCHEMA_ORDER_METADATA")) {
      k <- match(headers[c], keys)
      if (!is.na(k)) v <- x$schema$ordered_keys[[r - 1L]][[k]] else {
        if (source_id != "OUT_TABLE_S6" || sheet != "Fig5A_fold_summary" || headers[c] != "row_order")
          .lt_marine_abort("Unknown schema-order field; no inferred value.", "reporting_schema")
        v <- r - 1L
      }
      x$values[r, c] <- list(v); x$assigned[r, c] <- TRUE
    }
  }
  if (!all(x$assigned))
    .lt_marine_abort("Reporting sheet has unfilled/unclassified cells; no silent blank fallback.", "reporting_fields")
  # Verify imported/schema keys too, including the six historical reference rows.
  final_keys <- .lt_marine_report_keys(lapply(seq.int(2L, nrow(x$values)), function(r)
    as.list(x$values[r, key_cols])))
  if (!identical(final_keys, expected_keys))
    .lt_marine_abort("Final reporting key identity disagrees with the admitted sheet.", "reporting_keys")
  x$provenance <- list(authority_id = x$authority_id, authority_sha256 = x$authority_sha256,
    projection_source_sha256 = x$source_sha256, original_workbook_sha256 = x$schema$source_sha256,
    roles = as.list(table(x$roles)), imported_reference_qualifier = "imported_reporting_reference",
    historical_reference_qualifier = if (any(x$roles == "HISTORICAL_REFERENCE_IMPORTED_NOT_RECOMPUTED"))
      "historical_reference_imported_not_recomputed" else NULL,
    key_coverage_before_archived_order = TRUE, original_workbook_read = FALSE,
    no_computed_value_import = TRUE, scientific_or_target_pass_claimed = FALSE)
  x
}
