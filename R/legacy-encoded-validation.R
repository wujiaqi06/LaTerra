# Source-bound validation and replay-only handling for historical finite logGBI.

.lt_legacy_encoded_inference_message <- function() {
  paste(
    "Finite encoded logGBI is retained for historical/diagnostic replay only;",
    "it is not a La Terra V1 ordinary inference representation.",
    "Boundary-aware inference belongs to a future project; use GBI for",
    "La Terra V1 ordinary inference."
  )
}

.lt_refuse_legacy_encoded_inference <- function(x, route) {
  if (inherits(x, "lt_logGBI_encoded")) {
    .lt_abort_class(
      .lt_legacy_encoded_inference_message(),
      "lt_error_legacy_encoded_inference_not_supported",
      list(route = route, representation = "legacy_global_k10")
    )
  }
  invisible(x)
}

.lt_abort_legacy_source_mismatch <- function(message, field = NULL) {
  .lt_abort_class(
    message,
    "lt_error_legacy_encoded_source_mismatch",
    list(field = field)
  )
}

.lt_source_bound_encoded_fields <- function(source_log_relative) {
  list(
    source_semantic_identity = source_log_relative$semantic_identity,
    source_audit_identity = source_log_relative$audit_identity,
    source_representation_id = source_log_relative$representation_id,
    source_representation_contract_version =
      source_log_relative$representation_contract_version,
    source_canonicalization_version =
      source_log_relative$canonicalization_version,
    source_domain_identity = source_log_relative$domain_identity,
    source_ordered_axes_identity = source_log_relative$ordered_axes_identity,
    source_positive_log_values_sha256 =
      .lt_hash(source_log_relative$positive_log_values),
    source_ratio_state_sha256 = .lt_hash(source_log_relative$ratio_state),
    source_positive_domain_identity =
      source_log_relative$representation_provenance$positive_domain_identity,
    source_exact_zero_identity =
      source_log_relative$representation_provenance$exact_zero_identity
  )
}

#' Validate a finite encoded view against its supplied boundary source
#'
#' Structural validation of an `lt_logGBI_encoded` object checks only its local
#' consistency. This source-bound validator additionally derives every
#' scientific identity, state, positive value, and global-k10 coordinate from a
#' separately supplied `lt_log_relative` source. Local hashes are corruption
#' checks, not provenance authentication.
#'
#' @param encoded An `lt_logGBI_encoded` diagnostic/replay view.
#' @param source_log_relative Its separately supplied authoritative
#'   `lt_log_relative` dependency.
#' @return `encoded`, invisibly.
#' @export
validate_lt_logGBI_encoded_source <- function(encoded, source_log_relative) {
  validate_lt_logGBI_encoded(encoded)
  validate_lt_log_relative(source_log_relative)

  gene_index <- match(
    encoded$ordered_gene_ledger, source_log_relative$ordered_gene_ledger
  )
  branch_index <- match(
    encoded$ordered_branch_ledger, source_log_relative$ordered_branch_ledger
  )
  if (anyNA(gene_index) || anyNA(branch_index)) {
    .lt_abort_legacy_source_mismatch(
      "Encoded scientific axes are not an exact keyed subset of the supplied source.",
      "ordered_scientific_axes"
    )
  }

  source_fields <- .lt_source_bound_encoded_fields(source_log_relative)
  source_fields$source_object_identity <-
    .lt_legacy_source_object_identity(source_fields)
  for (field in names(source_fields)) {
    if (!identical(encoded[[field]], source_fields[[field]])) {
      .lt_abort_legacy_source_mismatch(
        paste0("Encoded source-bound field `", field,
               "` differs from the supplied source."),
        field
      )
    }
  }

  expected_state <- source_log_relative$ratio_state[
    gene_index, branch_index, drop = FALSE
  ]
  if (!identical(encoded$ratio_state, expected_state)) {
    .lt_abort_legacy_source_mismatch(
      "Encoded ratio states differ from the supplied source.",
      "ratio_state"
    )
  }

  m_ref <- source_log_relative$representation_provenance$original_source_m_ref
  L_min_plus <- source_log_relative$representation_provenance$L_min_plus
  sentinel <-
    source_log_relative$representation_provenance$legacy_global_k10_sentinel
  if (!identical(encoded$m_ref, m_ref) ||
      !identical(encoded$L_min_plus, L_min_plus) ||
      !identical(encoded$sentinel, sentinel) ||
      !identical(L_min_plus, log(m_ref)) ||
      !identical(sentinel, L_min_plus - 10)) {
    .lt_abort_legacy_source_mismatch(
      "Encoded global-k10 m_ref/L_min_plus/sentinel differs from source provenance.",
      "m_ref"
    )
  }

  expected_values <- source_log_relative$positive_log_values[
    gene_index, branch_index, drop = FALSE
  ]
  state <- as.integer(expected_state)
  expected_values[state == 2L] <- sentinel
  if (!identical(encoded$values, expected_values)) {
    .lt_abort_legacy_source_mismatch(
      "Encoded positive logs or exact-zero coordinates differ from the supplied source.",
      "values"
    )
  }

  expected_m_ref_fields <- source_fields
  expected_m_ref_fields$m_ref <- m_ref
  expected_m_ref_fields$L_min_plus <- L_min_plus
  expected_m_ref_fields$sentinel <- sentinel
  expected_m_ref_provenance <-
    .lt_legacy_m_ref_provenance(expected_m_ref_fields)
  if (!identical(encoded$m_ref_provenance, expected_m_ref_provenance) ||
      !identical(encoded$m_ref_provenance_identity,
                 .lt_hash(expected_m_ref_provenance))) {
    .lt_abort_legacy_source_mismatch(
      "Encoded m_ref provenance is not bound to the supplied source.",
      "m_ref_provenance"
    )
  }

  invisible(encoded)
}

.lt_legacy_replay_authority_fields <- function() {
  c(
    "authority_id", "purpose", "source_semantic_identity",
    "source_audit_identity", "source_object_identity", "encoded_identity",
    "mode", "m_ref", "L_min_plus", "sentinel",
    "positive_log_values_identity", "positive_domain_identity",
    "exact_zero_identity", "view_domain_identity",
    "view_ordered_axes_identity"
  )
}

#' Validate a frozen historical finite-logGBI replay snapshot
#'
#' This function performs identity validation only. It cannot compute an
#' association, calibration, p-value, or adjusted p-value. The replay authority
#' must be supplied independently; La Terra does not authenticate it with a
#' hidden token or package-local secret.
#'
#' @param encoded Candidate `lt_logGBI_encoded` snapshot.
#' @param source_log_relative Separately supplied authoritative boundary source.
#' @param frozen_replay_authority Externally frozen named list containing the
#'   fields returned in `authority` by a previous successful call.
#' @return A list with `status`, `purpose`, and an `authority` snapshot suitable
#'   for external freezing. This is validation evidence, not inference.
#' @export
validate_legacy_logGBI_replay <- function(
    encoded, source_log_relative, frozen_replay_authority) {
  validate_lt_logGBI_encoded_source(encoded, source_log_relative)
  if (!is.list(frozen_replay_authority) ||
      !identical(names(frozen_replay_authority),
                 .lt_legacy_replay_authority_fields())) {
    .lt_abort_class(
      "A separately frozen legacy replay authority with exact named fields is required.",
      "lt_error_legacy_replay_authority_invalid"
    )
  }
  if (!identical(frozen_replay_authority$purpose,
                 "HISTORICAL_REPLAY_VALIDATION") ||
      !identical(frozen_replay_authority$mode, "legacy_global_k10")) {
    .lt_abort_class(
      "Legacy replay authority purpose or mode is invalid.",
      "lt_error_legacy_replay_authority_invalid"
    )
  }
  .lt_assert_scalar_character(
    frozen_replay_authority$authority_id, "authority_id"
  )
  expected <- list(
    authority_id = frozen_replay_authority$authority_id,
    purpose = "HISTORICAL_REPLAY_VALIDATION",
    source_semantic_identity = encoded$source_semantic_identity,
    source_audit_identity = encoded$source_audit_identity,
    source_object_identity = encoded$source_object_identity,
    encoded_identity = encoded$encoded_identity,
    mode = encoded$mode,
    m_ref = encoded$m_ref,
    L_min_plus = encoded$L_min_plus,
    sentinel = encoded$sentinel,
    positive_log_values_identity = encoded$positive_log_values_identity,
    positive_domain_identity = encoded$positive_domain_identity,
    exact_zero_identity = encoded$exact_zero_identity,
    view_domain_identity = encoded$view_domain_identity,
    view_ordered_axes_identity = encoded$view_ordered_axes_identity
  )
  if (!identical(frozen_replay_authority, expected)) {
    .lt_abort_class(
      "Candidate encoded snapshot differs from the externally frozen replay authority.",
      "lt_error_legacy_replay_authority_mismatch"
    )
  }
  list(
    status = "EXACT_FROZEN_REPLAY_MATCH",
    purpose = "HISTORICAL_REPLAY_VALIDATION",
    inference_authorized = FALSE,
    authority = expected
  )
}
