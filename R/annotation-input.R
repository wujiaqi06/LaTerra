# Structural keyed-input reader, reusable beyond the frozen Marine fixture.
# Selecting declared columns does not admit a table as scientific authority.
.lt_read_keyed_character_fields <- function(path, columns, key) {
  .lt_assert_unique_ids(columns, "declared annotation columns")
  if (!is.character(key) || length(key) != 1L || is.na(key) || !key %in% columns) {
    stop("Annotation key must name one declared column.", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path) || dir.exists(path)) {
    stop("Supply an existing annotation table file.", call. = FALSE)
  }
  x <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character(), quote = "\"",
    comment.char = "", strip.white = FALSE)
  .lt_assert_unique_ids(names(x), "annotation column names")
  if (!all(columns %in% names(x))) {
    stop(paste0("Missing annotation column(s): ",
      paste(setdiff(columns, names(x)), collapse = ", ")), call. = FALSE)
  }
  .lt_assert_unique_ids(x[[key]], "annotation gene keys")
  # No sorting, normalization, coercion of annotation labels, or path resolution.
  x[, columns, drop = FALSE]
}

# This profile-specific input boundary is separate from structural validation.
# Its specification must come from the explicit frozen Marine profile, never
# from a guessed replacement table or an archived predictor/coefficient result.
.lt_marine_turnover_annotation_input <- function(path, specification) {
  fields <- c("gene_id", "candidate_module", "annotation_status",
              "module_annotation_source")
  required <- c("delivery_sha256", "full_source_sha256", "authority_id",
                "authority_sha256")
  if (!is.list(specification) || !all(required %in% names(specification))) {
    .lt_marine_abort("Turnover lookup requires its explicit source/view authority specification.",
      "annotation_authority")
  }
  for (n in required) {
    z <- specification[[n]]
    if (!is.character(z) || length(z) != 1L || is.na(z) || !nzchar(z)) {
      .lt_marine_abort(paste0("Invalid turnover annotation authority field: ", n),
        "annotation_authority")
    }
  }
  for (n in c("delivery_sha256", "full_source_sha256", "authority_sha256")) {
    if (!grepl("^[0-9a-f]{64}$", specification[[n]])) {
      .lt_marine_abort(paste0("Invalid SHA-256 in turnover annotation specification: ", n),
        "annotation_authority")
    }
  }
  if (!file.exists(path) || dir.exists(path)) {
    .lt_marine_abort(paste0("Missing turnover annotation view: ", path,
      ". Supply the admitted four-column view; historical paths are never opened."),
      "missing_input")
  }
  actual <- .lt_marine_sha(path)
  if (!identical(actual, specification$delivery_sha256)) {
    .lt_marine_abort("Frozen turnover annotation view does not match its declared snapshot.",
      "authority_mismatch")
  }
  view <- .lt_read_keyed_character_fields(path, fields, "gene_id")
  header <- names(utils::read.delim(path, nrows = 0L, check.names = FALSE,
    colClasses = "character", na.strings = character(), comment.char = ""))
  if (!identical(header, fields)) {
    .lt_marine_abort("The Marine turnover delivery must contain only the four admitted fields in order.",
      "annotation_fields")
  }
  list(annotation = view, provenance = list(
    path = normalizePath(path, mustWork = TRUE), delivery_sha256 = actual,
    full_source_sha256 = specification$full_source_sha256,
    authority_id = specification$authority_id,
    authority_sha256 = specification$authority_sha256,
    derivation = "named-column projection preserving source strings and row order",
    allowed_fields = fields, row_count = nrow(view),
    source_paths_are_inert_strings = TRUE,
    feature_selection_authority = FALSE,
    original_full_source_required_at_runtime = FALSE))
}
