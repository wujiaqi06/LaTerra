# First-class externally generated branch-state null ensembles.

.lt_null_availability_levels <- function() {
  c("available", "unavailable", "unresolved")
}

.lt_null_ensemble_payload <- function(x) {
  x[c(
    "ensemble_id", "trait_id", "branch_ids", "replicate_ids", "values",
    "type", "coding", "availability", "null_hypothesis",
    "generator_provenance", "observed_realization_included"
  )]
}

#' Construct an external branch-state null ensemble
#'
#' Each column is one complete null world that is applied to every gene. This
#' constructor records an externally supplied ensemble; it does not generate,
#' shuffle, or simulate branch states.
#'
#' @param ensemble_id Stable ensemble identifier.
#' @param trait_id Trait concept represented by the null worlds.
#' @param branch_ids Ordered unique scientific branch keys.
#' @param replicate_ids Ordered unique null-world identifiers.
#' @param values Branch by replicate atomic matrix.
#' @param type Either `binary` or `continuous`.
#' @param coding Named coding declaration. Binary ensembles require explicit
#'   `reference_level` and `focal_level` entries.
#' @param availability Branch by replicate matrix using `available`,
#'   `unavailable`, or `unresolved`. A branch-length vector is expanded across
#'   replicates as an explicit convenience.
#' @param null_hypothesis Non-empty scientific declaration supplied by the
#'   external generator authority.
#' @param generator_provenance Named generator provenance list.
#' @param observed_realization_included Must be `FALSE`. The observed branch
#'   state is not a null replicate; the empirical plus-one correction handles
#'   the observed realization mathematically.
#' @param metadata Named non-authoritative metadata list.
#' @return An `lt_branch_state_ensemble`.
#' @export
lt_branch_state_ensemble <- function(
    ensemble_id, trait_id, branch_ids, replicate_ids, values,
    type = c("binary", "continuous"), coding = list(), availability,
    null_hypothesis, generator_provenance = list(),
    observed_realization_included = FALSE, metadata = list()) {
  type <- match.arg(type)
  values <- as.matrix(values)
  if (is.null(dimnames(values))) {
    dimnames(values) <- list(branch_ids, replicate_ids)
  }
  if (is.atomic(availability) && is.null(dim(availability))) {
    if (length(availability) != length(branch_ids)) {
      .lt_abort("Vector `availability` must have exactly one entry per branch before replicate expansion.")
    }
    availability <- matrix(
      availability, nrow = length(branch_ids), ncol = length(replicate_ids),
      dimnames = list(branch_ids, replicate_ids)
    )
  } else {
    availability <- as.matrix(availability)
    if (is.null(dimnames(availability))) {
      dimnames(availability) <- list(branch_ids, replicate_ids)
    }
  }
  core <- list(
    ensemble_id = ensemble_id,
    trait_id = trait_id,
    branch_ids = branch_ids,
    replicate_ids = replicate_ids,
    values = values,
    type = type,
    coding = coding,
    availability = availability,
    null_hypothesis = null_hypothesis,
    generator_provenance = generator_provenance,
    observed_realization_included = observed_realization_included
  )
  out <- structure(c(core, list(
    null_ensemble_identity = .lt_hash(core),
    metadata = metadata
  )), class = "lt_branch_state_ensemble")
  validate_lt_branch_state_ensemble(out)
  out
}

#' Validate an external branch-state null ensemble
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_branch_state_ensemble <- function(x) {
  if (!inherits(x, "lt_branch_state_ensemble")) {
    .lt_abort("`x` must inherit from `lt_branch_state_ensemble`.")
  }
  .lt_assert_scalar_character(x$ensemble_id, "ensemble_id")
  .lt_assert_scalar_character(x$trait_id, "trait_id")
  .lt_assert_unique_ids(x$branch_ids, "branch_ids")
  .lt_assert_unique_ids(x$replicate_ids, "replicate_ids")
  if (!length(x$branch_ids) || !length(x$replicate_ids)) {
    .lt_abort("Null ensembles require at least one branch and one replicate.")
  }
  expected_dim <- c(length(x$branch_ids), length(x$replicate_ids))
  expected_dimnames <- list(x$branch_ids, x$replicate_ids)
  if (!is.matrix(x$values) || !is.atomic(x$values) || is.list(x$values) ||
      !identical(dim(x$values), expected_dim) ||
      !identical(dimnames(x$values), expected_dimnames)) {
    .lt_abort("Null ensemble values must be a branch x replicate matrix with exact keyed dimnames.")
  }
  if (!is.matrix(x$availability) || !is.character(x$availability) ||
      !identical(dim(x$availability), expected_dim) ||
      !identical(dimnames(x$availability), expected_dimnames) ||
      anyNA(x$availability) ||
      any(!x$availability %in% .lt_null_availability_levels())) {
    .lt_abort("Null ensemble availability must be an exact keyed branch x replicate matrix of known states.")
  }
  if (!x$type %in% c("binary", "continuous")) {
    .lt_abort("Null ensemble type must be binary or continuous.")
  }
  .lt_assert_named_list(x$coding, "coding")
  .lt_assert_scalar_character(x$null_hypothesis, "null_hypothesis")
  .lt_assert_named_list(x$generator_provenance, "generator_provenance")
  if (!is.logical(x$observed_realization_included) ||
      length(x$observed_realization_included) != 1L ||
      is.na(x$observed_realization_included)) {
    .lt_abort("`observed_realization_included` must be one non-missing logical value.")
  }
  if (isTRUE(x$observed_realization_included)) {
    .lt_abort("The observed branch state must not be inserted into the null ensemble.")
  }
  .lt_assert_named_list(x$metadata, "metadata")
  available <- x$availability == "available"
  if (any(is.na(x$values[available]))) {
    .lt_abort("Available null-ensemble cells require a present value.")
  }
  if (any(!is.na(x$values[!available]))) {
    .lt_abort("Unavailable or unresolved null-ensemble cells require NA values.")
  }
  if (x$type == "continuous") {
    if (!is.numeric(x$values) || any(!is.finite(x$values[available]))) {
      .lt_abort("Available continuous null values must be finite numeric values.")
    }
  } else {
    required <- c("reference_level", "focal_level")
    if (!all(required %in% names(x$coding))) {
      .lt_abort("Binary null coding requires explicit reference_level and focal_level semantics.")
    }
    ref <- .lt_level_key(x$coding$reference_level)
    focal <- .lt_level_key(x$coding$focal_level)
    if (identical(ref, focal)) {
      .lt_abort("Binary null reference and focal levels must be distinct.")
    }
    keys <- vapply(x$values[available], .lt_level_key, character(1L))
    if (any(!keys %in% c(ref, focal))) {
      .lt_abort("Available binary null values must match the declared coding.")
    }
  }
  .lt_assert_sha256(x$null_ensemble_identity, "null_ensemble_identity")
  if (!identical(x$null_ensemble_identity,
                 .lt_hash(.lt_null_ensemble_payload(x)))) {
    .lt_abort("Null-ensemble semantic identity is stale or inconsistent.")
  }
  invisible(x)
}

print.lt_branch_state_ensemble <- function(x, ...) {
  cat("<lt_branch_state_ensemble>", x$ensemble_id, "\n")
  cat("  trait:", x$trait_id, " type:", x$type, "\n")
  cat("  branches:", length(x$branch_ids),
      " replicates:", length(x$replicate_ids), "\n")
  cat("  observed realization included: no\n")
  cat("  generator scientific validity: external / not adjudicated by T2B\n")
  invisible(x)
}
