# INTEROP001 consumer boundary. No coordinate construction or scientific fitting.
.lt_sx_need <- function(ok, message, code = "structure") {
  if (!isTRUE(ok)) stop(structure(list(message = paste0("SplitAligner exchange: ", message),
    call = NULL), class = c(paste0("lt_error_exchange_", code), "error", "condition")))
}

.lt_sx_fields <- function(x, fields, where) {
  .lt_sx_need(is.list(x) && !is.null(names(x)) && !anyDuplicated(names(x)) &&
    setequal(names(x), fields), paste(where, "has missing, duplicate or unknown fields."))
}

.lt_sx_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
.lt_sx_keys <- function(x, empty = FALSE) is.character(x) && (empty || length(x) > 0L) &&
  !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
.lt_sx_equal <- function(x, y, where) .lt_sx_need(identical(x, y), paste(where, "disagrees with its authoritative source."))
.lt_sx_columns <- function(x, columns, where) {
  .lt_sx_need(is.data.frame(x) && !anyDuplicated(names(x)) && all(columns %in% names(x)),
    paste(where, "is missing required data-frame columns."))
}
.lt_sx_components <- function() c("schema", "producer", "contract", "axes", "identity",
  "tree", "coordinates", "labels", "data", "preserved", "annotation")
.lt_sx_arrays <- function() c("values", "coord_state", "value_reason", "source_state", "source_reason", "numeric_status")

.lt_sx_header <- function(x) {
  .lt_sx_need(inherits(x, "splitaligner_exchange"), "Expected a versioned splitaligner_exchange RDS, not a raw matrix.")
  .lt_sx_fields(x, c(.lt_sx_components(), "manifest"), "envelope")
  .lt_sx_need(identical(x$schema, list(id = "SplitAlignerR.exchange",
    version = "0.2.0-development", profile = "single-fixed-primitive-branch-length-v1")),
    "Unsupported schema/version/profile. Old 0.1.0-development has H01 and is refused without resealing; this reader requires explicit 0.2.0-development.", "version")
  .lt_sx_fields(x$contract, c("mode", "payload", "state_schema", "source_state_schema", "split_key_schema", "numeric_policy"), "contract")
  .lt_sx_need(identical(x$contract$mode, "fixed"), "Only single fixed mode is supported; free/paired are not admitted.", "mode")
  expected <- list(state_schema = "LT-compatible-observed-plus-structural-v1",
    source_state_schema = "single-tree-state-v1", split_key_schema = "SplitAligner-canonical-split-key-v2",
    numeric_policy = "finite-double-v1; negative finite values preserved; LT S2 admission is separate")
  for (n in names(expected)) .lt_sx_equal(x$contract[[n]], expected[[n]], paste("contract", n))
  p <- x$contract$payload
  .lt_sx_fields(p, c("type", "name", "units", "scale", "origin", "spec_id"), "payload")
  .lt_sx_need(all(vapply(p, .lt_sx_string, logical(1))) && identical(p$type, "branch_length") &&
    identical(p$scale, "identity"), "Only explicit identity-scale branch-length payloads are supported.", "payload")
  .lt_sx_fields(x$producer, c("package", "version", "core_version", "source"), "producer")
  .lt_sx_need(identical(x$producer$package, "SplitAlignerR") && .lt_sx_string(x$producer$version) &&
    .lt_sx_string(x$producer$core_version), "Missing or unsupported producer identity.")
  p <- x$producer$source
  .lt_sx_fields(p, c("source_commit", "source_tree", "source_parent", "dirty_scope", "source_files", "provenance_status"), "producer source")
  .lt_sx_need(all(vapply(p[c("source_commit", "source_tree", "source_parent")], function(v)
    .lt_sx_string(v) && grepl("^[0-9a-f]{40}$", v), logical(1))), "Producer Git identities must be full hashes.")
  .lt_sx_need(.lt_sx_string(p$provenance_status) && p$provenance_status %in% c("development_unreleased", "synthetic_fixture") &&
    is.character(p$dirty_scope) && !anyNA(p$dirty_scope) && !anyDuplicated(p$dirty_scope) &&
    is.character(p$source_files) && .lt_sx_keys(names(p$source_files)) && !anyNA(p$source_files) &&
    all(grepl("^[0-9a-f]{64}$", p$source_files)), "Unknown or incomplete producer provenance.")
  invisible(x)
}

# Validate ALL ORIGINAL arrays before lt_matrix can replace any dimnames.
.lt_sx_data <- function(x) {
  a <- x$axes
  .lt_sx_fields(a, c("orientation", "branch_margin", "gene_ids", "coordinate_ids"), "axes")
  .lt_sx_need(identical(a$branch_margin, 1L) || identical(a$branch_margin, 2L), "An explicit integer branch_margin is required.", "axis")
  .lt_sx_equal(a$orientation, if (a$branch_margin == 2L) "gene_by_branch" else "branch_by_gene", "orientation")
  .lt_sx_need(.lt_sx_keys(a$gene_ids) && .lt_sx_keys(a$coordinate_ids), "Gene and scientific coordinate keys must be unique and nonempty.", "axis")
  .lt_sx_fields(x$data, .lt_sx_arrays(), "data")
  axes <- if (a$branch_margin == 2L) list(a$gene_ids, a$coordinate_ids) else list(a$coordinate_ids, a$gene_ids)
  for (n in .lt_sx_arrays()) {
    v <- x$data[[n]]
    .lt_sx_need(is.matrix(v) && if (n == "values") is.double(v) else is.character(v),
      paste(n, "has an unsupported matrix storage type."), "axis")
    .lt_sx_need(identical(dim(v), as.integer(lengths(axes))) && identical(dimnames(v), axes),
      paste(n, "must have the complete original ordered axes; no per-array repair is allowed."), "axis")
  }
  d <- x$data[.lt_sx_arrays()]
  if (a$branch_margin == 1L) d <- lapply(d, t) # One explicit operation on all six.
  .lt_sx_need(!any(is.nan(d$values) | is.infinite(d$values)) && !any(d$values < 0, na.rm = TRUE),
    "S2 admission rejects negative, NaN or infinite branch lengths; no clipping or replacement.", "payload")
  .lt_sx_need(!anyNA(d$source_state) && all(d$source_state %in% c("mapped", "NA_struct", "NA_fuse", "NA_topo")),
    "Unknown or missing source state.", "state")
  state <- d$source_state; state[state == "mapped"] <- "observed"
  .lt_sx_equal(d$coord_state, state, "coordinate state mapping")
  observed <- state == "observed"; available <- is.finite(d$values)
  .lt_sx_need(!any(!observed & available), "Coordinate absence cannot contain a numeric payload.", "state")
  .lt_sx_need(!anyNA(d$source_reason) && all(nzchar(d$source_reason)) &&
    !anyNA(d$numeric_status) && all(nzchar(d$numeric_status)), "Source reason and numeric status must be preserved explicitly.", "state")
  missing <- observed & !available
  .lt_sx_need(all(d$numeric_status[missing] %in% c("missing_branch_length", "software_failure_marker")),
    "Mapped-unavailable payload needs a supported explicit numeric reason.", "state")
  reason <- matrix(NA_character_, nrow(d$values), ncol(d$values), dimnames = dimnames(d$values))
  reason[missing] <- d$numeric_status[missing]
  .lt_sx_equal(d$value_reason, reason, "current-payload value reason")
  d
}

# The producer owns canonical encoding, native B/reference agreement and graph
# semantics. No local reimplementation or caller-supplied validation callback.
.lt_sx_backend <- function() {
  .lt_sx_need(requireNamespace("SplitAlignerR", quietly = TRUE),
    "Install the supported SplitAlignerR 0.1.0.9002 development backend in your chosen library.", "dependency")
  .lt_sx_need(identical(as.character(utils::packageVersion("SplitAlignerR")), "0.1.0.9002"),
    "Unsupported SplitAlignerR backend version; no silent version fallback.", "dependency")
  schema <- SplitAlignerR::splitaligner_exchange_schema()
  .lt_sx_need(identical(schema$schema_version, "0.2.0-development") &&
    identical(schema$transport$content_hash, "SAR-C14N-1-SHA256"), "Backend schema/hash contract is unsupported.", "dependency")
  root <- find.package("SplitAlignerR")
  files <- c("DESCRIPTION", "NAMESPACE", "R/SplitAlignerR.rdb", "R/SplitAlignerR.rdx",
    paste0("libs/", list.files(file.path(root, "libs"), recursive = TRUE)))
  .lt_sx_need(all(file.exists(file.path(root, files))), "Installed backend runtime files are incomplete.", "dependency")
  list(package = "SplitAlignerR", version = "0.1.0.9002", path = root,
    schema_version = "0.2.0-development", hash_schema = "SAR-C14N-1-SHA256",
    public_validator = "SplitAlignerR::validate_splitaligner_exchange",
    runtime_files = data.frame(path = files, sha256 = vapply(file.path(root, files), function(p)
      digest::digest(file = p, algo = "sha256", serialize = FALSE), character(1)), row.names = NULL))
}

.lt_sx_validate_upstream <- function(x) {
  backend <- .lt_sx_backend()
  # Mandatory even when the consumer's structural relationships already agree.
  tryCatch(SplitAlignerR::validate_splitaligner_exchange(x), error = function(e)
    .lt_sx_need(FALSE, paste("Upstream public validation rejected the envelope:", conditionMessage(e)), "upstream"))
  backend
}

#' Read a versioned SplitAlignerR exchange into explicit S2 inputs
#'
#' Experimental development interface, not Marine/release certification.
#' Supports only the declared single-fixed primitive branch-length profile.
#' All six original arrays are checked before construction; a declared transpose
#' transforms all six together. Coordinates are scientific keys, not B aliases.
#' Complete original exchange, result, tree and ledgers are retained unchanged.
#'
#' Schema 0.2.0-development requires the supported SplitAlignerR 0.1.0.9002
#' public validator for canonical hashes and native label/reference agreement.
#' Old 0.1.0-development (INTEROP001-H01) and unknown versions are refused,
#' without resealing, automatic conversion, locale changes or validator bypass.
#' Scoped H01 evidence covers R4.4.2 C/UTF-8 only, not other R versions or release.
#' This import boundary does not add golden SHA requirements to ordinary matrices.
#'
#' @param file Existing exchange RDS. Historical source paths in it are never
#'   reopened. Read only files from a trusted producer.
#' @param allow_synthetic_migration Explicitly permit producer-declared synthetic
#'   label-migration fixtures. Never authorizes real historical migration.
#' @return A list with `matrix` (`lt_matrix`), full `tree`, S2 `coordinates`,
#'   unchanged `exchange`, and import `provenance`. Save this complete list to
#'   preserve all upstream evidence. Traits, groups, folds and recipes must still
#'   be supplied explicitly to [lt_run_s2()]; reading performs no analysis.
#' @export
lt_read_splitaligner_exchange <- function(file, allow_synthetic_migration = FALSE) {
  .lt_sx_need(.lt_sx_string(file) && file.exists(file) && !dir.exists(file), "Supply an existing exchange RDS file.", "file")
  .lt_sx_need(identical(allow_synthetic_migration, FALSE) || identical(allow_synthetic_migration, TRUE),
    "allow_synthetic_migration must be TRUE or FALSE.")
  file <- normalizePath(file, mustWork = TRUE)
  before <- digest::digest(file = file, algo = "sha256", serialize = FALSE)
  x <- readRDS(file)
  .lt_sx_need(identical(before, digest::digest(file = file, algo = "sha256", serialize = FALSE)),
    "File changed during import.", "file")
  .lt_sx_header(x)
  data <- .lt_sx_data(x)
  out <- .lt_sx_adapt(x, data, allow_synthetic_migration)
  # No public option skips this independent producer/native/canonical boundary.
  out$provenance$validation_backend <- .lt_sx_validate_upstream(x)
  out$provenance$input_file <- file
  out$provenance$input_file_sha256 <- before
  out$provenance$synthetic_migration_allowed <- allow_synthetic_migration
  out$matrix$metadata$splitaligner_import <- out$provenance
  out
}
