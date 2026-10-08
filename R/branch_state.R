# First-class branch-associated state objects.

.lt_branch_state_availability_levels <- function() {
  c("available", "unavailable", "unresolved")
}

.lt_branch_state_semantic_payload <- function(x) {
  x[c(
    "state_id", "trait_id", "branch_ids", "values", "type",
    "availability", "reason", "coding", "provenance"
  )]
}

.lt_level_key <- function(x) {
  if (length(x) != 1L || is.na(x)) {
    .lt_abort("Binary coding levels must be non-missing scalar values.")
  }
  if (is.factor(x)) {
    return(paste0("character:", encodeString(as.character(x), quote = "\"")))
  }
  if (is.numeric(x)) {
    if (!is.finite(x)) {
      .lt_abort("Binary numeric coding levels must be finite.")
    }
    return(paste0("numeric:", sprintf("%a", as.double(x))))
  }
  paste0(typeof(x), ":", encodeString(as.character(x), quote = "\""))
}

#' Construct a branch-associated state object
#'
#' A branch state is not a taxon trait. Scientific branch identifiers are
#' exact keys and are never normalized, case-folded, or guessed. Availability
#' and its reason are explicit so numeric `NA` is not asked to carry scientific
#' meaning by itself.
#'
#' @param state_id Stable branch-state identifier.
#' @param trait_id Identifier of the source trait concept.
#' @param branch_ids Ordered, unique scientific branch keys.
#' @param values Atomic values aligned to `branch_ids`.
#' @param type One of `binary`, `continuous`, `ordinal`, or `categorical`.
#' @param availability Per-branch `available`, `unavailable`, or `unresolved`.
#' @param reason Per-branch reason; `NA` for available values and non-empty for
#'   unavailable or unresolved values.
#' @param coding Named coding declaration. Binary states require explicit
#'   `reference_level` and `focal_level` entries.
#' @param provenance Named source/provenance list. File hashes are optional for
#'   ordinary runtime use.
#' @param metadata Named non-authoritative metadata list.
#' @return An `lt_branch_state`.
#' @export
lt_branch_state <- function(state_id, trait_id, branch_ids, values,
                            type = c("binary", "continuous", "ordinal",
                                     "categorical"),
                            availability, reason, coding = list(),
                            provenance = list(), metadata = list()) {
  type <- match.arg(type)
  core <- list(
    state_id = state_id,
    trait_id = trait_id,
    branch_ids = branch_ids,
    values = values,
    type = type,
    availability = availability,
    reason = reason,
    coding = coding,
    provenance = provenance
  )
  x <- structure(c(core, list(
    branch_state_identity = .lt_hash(core),
    metadata = metadata
  )), class = "lt_branch_state")
  validate_lt_branch_state(x)
  x
}

#' Validate a branch-associated state object
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_branch_state <- function(x) {
  if (!inherits(x, "lt_branch_state")) {
    .lt_abort("`x` must inherit from `lt_branch_state`.")
  }
  .lt_assert_scalar_character(x$state_id, "state_id")
  .lt_assert_scalar_character(x$trait_id, "trait_id")
  .lt_assert_unique_ids(x$branch_ids, "branch_ids")
  n <- length(x$branch_ids)
  if (!is.atomic(x$values) || is.list(x$values) || length(x$values) != n) {
    .lt_abort("Branch-state `values` must be an atomic vector aligned to `branch_ids`.")
  }
  if (!is.character(x$availability) || length(x$availability) != n ||
      anyNA(x$availability) ||
      any(!x$availability %in% .lt_branch_state_availability_levels())) {
    .lt_abort("Branch-state availability uses an unknown state or is not aligned to branch IDs.")
  }
  if (!is.character(x$reason) || length(x$reason) != n) {
    .lt_abort("Branch-state `reason` must be a character vector aligned to branch IDs.")
  }
  if (!x$type %in% c("binary", "continuous", "ordinal", "categorical")) {
    .lt_abort("Unknown branch-state type.")
  }
  .lt_assert_named_list(x$coding, "coding")
  .lt_assert_named_list(x$provenance, "provenance")
  .lt_assert_named_list(x$metadata, "metadata")

  available <- x$availability == "available"
  absent <- !available
  if (any(is.na(x$values[available]))) {
    .lt_abort("Available branch states require a present value.")
  }
  if (any(!is.na(x$values[absent]))) {
    .lt_abort("Unavailable or unresolved branch states require an NA value.")
  }
  if (any(!is.na(x$reason[available]))) {
    .lt_abort("Available branch states require reason = NA.")
  }
  if (any(is.na(x$reason[absent]) | !nzchar(trimws(x$reason[absent])))) {
    .lt_abort("Unavailable or unresolved branch states require a non-empty reason.")
  }

  if (x$type == "continuous") {
    if (!is.numeric(x$values)) {
      .lt_abort("Continuous branch-state values must be numeric.")
    }
    if (any(!is.finite(x$values[available]))) {
      .lt_abort("Available continuous branch-state values must be finite.")
    }
  } else if (is.numeric(x$values) &&
             any(is.infinite(x$values[available]) | is.nan(x$values[available]))) {
    .lt_abort("Branch-state values cannot contain NaN or Inf.")
  }

  if (x$type == "binary") {
    required <- c("reference_level", "focal_level")
    if (!all(required %in% names(x$coding))) {
      .lt_abort("Binary coding requires explicit reference_level and focal_level semantics.")
    }
    reference_key <- .lt_level_key(x$coding$reference_level)
    focal_key <- .lt_level_key(x$coding$focal_level)
    if (identical(reference_key, focal_key)) {
      .lt_abort("Binary reference_level and focal_level must be distinct.")
    }
    value_keys <- vapply(x$values[available], .lt_level_key, character(1L))
    if (any(!value_keys %in% c(reference_key, focal_key))) {
      .lt_abort("Available binary values must match the explicitly declared reference/focal levels.")
    }
  }

  expected_identity <- .lt_hash(.lt_branch_state_semantic_payload(x))
  .lt_assert_sha256(x$branch_state_identity, "branch_state_identity")
  if (!identical(x$branch_state_identity, expected_identity)) {
    .lt_abort("Branch-state semantic identity is stale or inconsistent.")
  }
  invisible(x)
}
