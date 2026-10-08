#' Allowed matrix state labels
#'
#' Returns the coordinate-state vocabulary frozen for the LT001-A1A matrix
#' contract. Downstream value-reason codes are deliberately not global.
#'
#' @return A character vector.
#' @export
lt_allowed_states <- function() {
  c("observed", "NA_struct", "NA_fuse", "NA_topo", "residual_NA")
}

.lt_runtime_coordinate_provenance <- function(provenance, branch_ids) {
  required <- c(
    "coordinate_system", "reference_tree_sha256",
    "ordered_branch_ledger_sha256", "ordered_taxon_ledger_sha256",
    "branch_label_contract", "source_method"
  )
  if (is.null(provenance)) provenance <- list()
  .lt_assert_named_list(provenance, "coordinate_provenance")
  if (!length(setdiff(required, names(provenance)))) return(provenance)
  defaults <- list(
    coordinate_system = "user_supplied_keyed_matrix",
    reference_tree_sha256 = NA_character_,
    ordered_branch_ledger_sha256 = .lt_hash(branch_ids),
    ordered_taxon_ledger_sha256 = NA_character_,
    branch_label_contract = "user_supplied_scientific_branch_keys",
    source_method = "user_supplied_matrix",
    audit_identity_status = "UNFROZEN_RUNTIME_OBJECT"
  )
  utils::modifyList(defaults, provenance)
}

.lt_gene_display_labels <- function(x) {
  x$gene_display_labels %||% x$gene_ids
}

.lt_branch_display_labels <- function(x) {
  x$branch_display_labels %||% x$branch_ids
}

#' Construct a coordinate-valid La Terra matrix
#'
#' `lt_matrix()` stores numeric values separately from coordinate truth and
#' current-payload absence reasons. It stores coordinate provenance but never
#' recomputes branch identity.
#'
#' @param values Dense base-R numeric gene-by-branch matrix.
#' @param coord_state Character coordinate-state matrix with the same dimensions
#'   as `values`.
#' @param payload Named typed-payload list containing `type`, `name`, `units`,
#'   `scale`, `origin`, and `spec_id`.
#' @param value_reason Character matrix describing why an observed coordinate
#'   lacks a value for the current payload. Defaults to all `NA_character_`.
#' @param gene_ids Ordered unique gene identifiers.
#' @param branch_ids Ordered unique branch-coordinate identifiers.
#' @param coordinate_provenance Optional named list containing
#'   `coordinate_system`, reference-tree and ordered-ledger SHA-256 identities,
#'   `branch_label_contract`, and `source_method`. Ordinary user matrices may
#'   omit this argument; La Terra creates an explicitly unfrozen runtime
#'   provenance record and never requires a pre-existing golden SHA-256.
#' @param metadata Optional named metadata list.
#' @param gene_display_labels Optional mutable display labels. Scientific gene
#'   identity remains `gene_ids`.
#' @param branch_display_labels Optional mutable display labels. Scientific
#'   coordinate identity remains `branch_ids`.
#'
#' @return An object of class `lt_matrix`.
#' @export
lt_matrix <- function(values,
                      coord_state,
                      payload,
                      value_reason = matrix(
                        NA_character_,
                        nrow = nrow(values),
                        ncol = ncol(values),
                        dimnames = dimnames(values)
                      ),
                      gene_ids = rownames(values),
                      branch_ids = colnames(values),
                      coordinate_provenance = NULL,
                      metadata = list(),
                      gene_display_labels = gene_ids,
                      branch_display_labels = branch_ids) {
  if (!is.matrix(values) || !is.numeric(values)) {
    .lt_abort("`values` must be a numeric matrix.")
  }
  if (!is.matrix(coord_state) || !is.character(coord_state)) {
    .lt_abort("`coord_state` must be a character matrix.")
  }
  if (!is.matrix(value_reason) || !is.character(value_reason)) {
    .lt_abort("`value_reason` must be a character matrix.")
  }
  if (!identical(dim(values), dim(coord_state)) ||
      !identical(dim(values), dim(value_reason))) {
    .lt_abort("Matrix/coord_state/value_reason dimension mismatch.")
  }
  if (nrow(values) < 1L || ncol(values) < 1L) {
    .lt_abort("Matrix fields must have at least one gene and one branch.")
  }
  if (is.null(gene_ids) || is.null(branch_ids)) {
    .lt_abort("`gene_ids` and `branch_ids` must be supplied or present as dimnames.")
  }

  dimnames(values) <- list(gene_ids, branch_ids)
  dimnames(coord_state) <- list(gene_ids, branch_ids)
  dimnames(value_reason) <- list(gene_ids, branch_ids)
  coordinate_provenance <- .lt_runtime_coordinate_provenance(
    coordinate_provenance, branch_ids
  )

  x <- structure(
    list(
      values = values,
      coord_state = coord_state,
      value_reason = value_reason,
      payload = payload,
      gene_ids = gene_ids,
      branch_ids = branch_ids,
      gene_display_labels = gene_display_labels,
      branch_display_labels = branch_display_labels,
      coordinate_provenance = coordinate_provenance,
      metadata = metadata
    ),
    class = "lt_matrix"
  )
  validate_lt_matrix(x)
  x
}

#' Validate a La Terra matrix
#'
#' @param x Object to validate.
#' @return `x`, invisibly, or an error describing the violated invariant.
#' @export
validate_lt_matrix <- function(x) {
  if (!inherits(x, "lt_matrix")) {
    .lt_abort("`x` must inherit from `lt_matrix`.")
  }
  if (!is.matrix(x$values) || !is.numeric(x$values)) {
    .lt_abort("`x$values` must be a numeric matrix.")
  }
  if (!is.matrix(x$coord_state) || !is.character(x$coord_state)) {
    .lt_abort("`x$coord_state` must be a character matrix.")
  }
  if (!is.matrix(x$value_reason) || !is.character(x$value_reason)) {
    .lt_abort("`x$value_reason` must be a character matrix.")
  }
  if (!identical(dim(x$values), dim(x$coord_state)) ||
      !identical(dim(x$values), dim(x$value_reason))) {
    .lt_abort("Matrix/coord_state/value_reason dimension mismatch.")
  }
  if (nrow(x$values) < 1L || ncol(x$values) < 1L) {
    .lt_abort("Matrix fields must have at least one gene and one branch.")
  }
  if (length(x$gene_ids) != nrow(x$values)) {
    .lt_abort("`gene_ids` length must equal the number of matrix rows.")
  }
  if (length(x$branch_ids) != ncol(x$values)) {
    .lt_abort("`branch_ids` length must equal the number of matrix columns.")
  }
  .lt_assert_unique_ids(x$gene_ids, "gene_ids")
  .lt_assert_unique_ids(x$branch_ids, "branch_ids")

  gene_display_labels <- .lt_gene_display_labels(x)
  branch_display_labels <- .lt_branch_display_labels(x)
  if (!is.character(gene_display_labels) ||
      length(gene_display_labels) != length(x$gene_ids) ||
      anyNA(gene_display_labels) || any(!nzchar(gene_display_labels))) {
    .lt_abort("`gene_display_labels` must be non-empty labels aligned to gene scientific keys.")
  }
  if (!is.character(branch_display_labels) ||
      length(branch_display_labels) != length(x$branch_ids) ||
      anyNA(branch_display_labels) || any(!nzchar(branch_display_labels))) {
    .lt_abort("`branch_display_labels` must be non-empty labels aligned to branch scientific keys.")
  }

  if (!identical(rownames(x$values), x$gene_ids) ||
      !identical(colnames(x$values), x$branch_ids) ||
      !identical(dimnames(x$coord_state), dimnames(x$values)) ||
      !identical(dimnames(x$value_reason), dimnames(x$values))) {
    .lt_abort("Matrix fields and ordered ledgers must have identical axes.")
  }

  if (anyNA(x$coord_state)) {
    .lt_abort("`coord_state` must use an explicit coordinate-state label.")
  }
  unknown <- setdiff(unique(as.vector(x$coord_state)), lt_allowed_states())
  if (length(unknown) > 0L) {
    .lt_abort(sprintf("Unknown coordinate-state labels: %s.", paste(unknown, collapse = ", ")))
  }

  if (any(is.nan(x$values)) || any(is.infinite(x$values), na.rm = TRUE)) {
    .lt_abort("`values` must not contain NaN or infinite values.")
  }
  observed <- x$coord_state == "observed"
  finite_value <- !is.na(x$values)
  explicit_reason <- !is.na(x$value_reason) & nzchar(trimws(x$value_reason))

  if (any(!observed & finite_value)) {
    .lt_abort("Every non-observed coordinate state must have an NA value.")
  }
  if (any(!observed & !is.na(x$value_reason))) {
    .lt_abort("Coordinate-absent cells must not carry a current-payload value reason.")
  }
  if (any(observed & finite_value & !is.na(x$value_reason))) {
    .lt_abort("Observed finite values must have an NA `value_reason`.")
  }
  if (any(observed & !finite_value & !explicit_reason)) {
    .lt_abort("Observed coordinates with NA values require a non-empty `value_reason`.")
  }

  .lt_validate_payload(x$payload)

  .lt_assert_named_list(x$coordinate_provenance, "coordinate_provenance")
  required_provenance <- c(
    "coordinate_system", "reference_tree_sha256",
    "ordered_branch_ledger_sha256", "ordered_taxon_ledger_sha256",
    "branch_label_contract", "source_method"
  )
  missing_provenance <- setdiff(required_provenance, names(x$coordinate_provenance))
  if (length(missing_provenance) > 0L) {
    .lt_abort(sprintf(
      "`coordinate_provenance` is missing: %s.",
      paste(missing_provenance, collapse = ", ")
    ))
  }
  hash_fields <- c(
    "reference_tree_sha256", "ordered_branch_ledger_sha256",
    "ordered_taxon_ledger_sha256"
  )
  for (field in setdiff(required_provenance, hash_fields)) {
    .lt_assert_scalar_character(
      x$coordinate_provenance[[field]],
      paste0("coordinate_provenance$", field)
    )
  }
  for (field in hash_fields) {
    value <- x$coordinate_provenance[[field]]
    if (!is.character(value) || length(value) != 1L ||
        (!is.na(value) && !grepl("^[0-9a-f]{64}$", value))) {
      .lt_abort(paste0("`coordinate_provenance$", field,
                       "` must be SHA-256 or NA for an unfrozen runtime object."))
    }
  }
  optional_identity <- intersect(
    c("source_software", "source_version", "source_commit"),
    names(x$coordinate_provenance)
  )
  for (field in optional_identity) {
    .lt_assert_scalar_character(
      x$coordinate_provenance[[field]],
      paste0("coordinate_provenance$", field)
    )
  }
  .lt_assert_named_list(x$metadata, "metadata")
  invisible(x)
}

#' Change display labels without changing scientific keys
#'
#' Display labels are deliberately outside scientific dependency identities and
#' audit-semantic hashes. This helper returns a structurally validated copy;
#' numerical fits remain compatible because `gene_ids` and `branch_ids` are
#' unchanged.
#'
#' @param x An `lt_matrix`, `lt_rate`, or `lt_rate_common_domain`.
#' @param gene_labels Optional labels aligned to scientific gene keys.
#' @param branch_labels Optional labels aligned to scientific branch keys.
#' @return A copy of `x` with updated display labels.
#' @export
lt_set_display_labels <- function(x, gene_labels = NULL, branch_labels = NULL) {
  if (!inherits(x, c("lt_matrix", "lt_rate", "lt_rate_common_domain"))) {
    .lt_abort(paste(
      "Display labels are supported for lt_matrix, lt_rate, and",
      "lt_rate_common_domain objects."
    ))
  }
  gene_keys <- if (inherits(x, "lt_matrix")) x$gene_ids else
    x$ordered_gene_ledger
  branch_keys <- if (inherits(x, "lt_matrix")) x$branch_ids else
    x$ordered_branch_ledger
  if (!is.null(gene_labels)) {
    if (!is.character(gene_labels) || length(gene_labels) != length(gene_keys) ||
        anyNA(gene_labels) || any(!nzchar(gene_labels))) {
      .lt_abort("`gene_labels` must be non-empty and aligned to gene scientific keys.")
    }
    x$gene_display_labels <- gene_labels
  }
  if (!is.null(branch_labels)) {
    if (!is.character(branch_labels) ||
        length(branch_labels) != length(branch_keys) ||
        anyNA(branch_labels) || any(!nzchar(branch_labels))) {
      .lt_abort("`branch_labels` must be non-empty and aligned to branch scientific keys.")
    }
    x$branch_display_labels <- branch_labels
  }
  if (inherits(x, "lt_matrix")) {
    validate_lt_matrix(x)
  } else if (inherits(x, "lt_rate")) {
    validate_lt_rate(x)
  } else {
    validate_lt_rate_common_domain(x)
  }
  x
}

.lt_validate_payload <- function(payload) {
  .lt_assert_named_list(payload, "payload")
  required <- c("type", "name", "units", "scale", "origin", "spec_id")
  missing <- setdiff(required, names(payload))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("`payload` is missing: %s.", paste(missing, collapse = ", ")))
  }
  for (field in c("type", "name", "units", "scale", "origin")) {
    .lt_assert_scalar_character(payload[[field]], paste0("payload$", field))
  }
  if (!payload$origin %in% c("imported", "derived")) {
    .lt_abort("`payload$origin` must be `imported` or `derived`.")
  }
  .lt_assert_nullable_scalar_character(payload$spec_id, "payload$spec_id")
  if (identical(payload$origin, "derived") && is.na(payload$spec_id)) {
    .lt_abort("A derived payload requires a non-empty `payload$spec_id`.")
  }
  invisible(payload)
}
