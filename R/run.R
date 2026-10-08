#' Empty La Terra provenance container
#'
#' Creates the complete LT001-A1A provenance schema without computing hashes,
#' ledgers, or environment records.
#'
#' @return A named list with all required future provenance fields.
#' @export
lt_empty_provenance <- function() {
  list(
    input_sha256 = data.frame(
      role = character(), path = character(), sha256 = character(),
      stringsAsFactors = FALSE
    ),
    configuration_sha256 = NA_character_,
    software_version = list(package = "LaTerra", version = NA_character_),
    environment = list(),
    ordered_gene_ledger = data.frame(
      order = integer(), gene_id = character(), stringsAsFactors = FALSE
    ),
    ordered_branch_ledger = data.frame(
      order = integer(), branch_id = character(), stringsAsFactors = FALSE
    ),
    taxon_domain_ledger = data.frame(
      order = integer(), taxon_id = character(), eligible = logical(),
      exclusion_reason = character(), stringsAsFactors = FALSE
    ),
    eligibility_masks = list(),
    specification_ledger = data.frame(
      specification_id = character(), kind = character(), version = character(),
      definition = character(), sha256 = character(), stringsAsFactors = FALSE
    ),
    filter_ledger = data.frame(
      order = integer(), stage = character(), decision_id = character(),
      target_type = character(), target_id = character(), decision = character(),
      reason = character(), stringsAsFactors = FALSE
    ),
    seed_ledger = data.frame(
      scope = character(), seed = integer(), rng_kind = character(),
      normal_kind = character(), sample_kind = character(), r_version = character(),
      stringsAsFactors = FALSE
    ),
    fold_ledger = data.frame(
      fold_id = character(), taxon_id = character(), role = character(),
      group = character(),
      stringsAsFactors = FALSE
    ),
    commands = character(),
    output_hashes = data.frame(
      role = character(), path = character(), sha256 = character(),
      stringsAsFactors = FALSE
    )
  )
}

#' Validate a La Terra provenance container
#'
#' This validator freezes audit-container structure only. Populating a valid
#' container or completing a run does not constitute scientific certification.
#'
#' @param provenance Named provenance list.
#' @return `provenance`, invisibly.
#' @export
validate_lt_provenance <- function(provenance) {
  .lt_assert_named_list(provenance, "provenance")
  template <- lt_empty_provenance()
  missing <- setdiff(names(template), names(provenance))
  unknown <- setdiff(names(provenance), names(template))
  if (length(missing) > 0L || length(unknown) > 0L) {
    detail <- c(
      if (length(missing) > 0L) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(unknown) > 0L) paste0("unknown: ", paste(unknown, collapse = ", "))
    )
    .lt_abort(paste0("Invalid provenance schema (", paste(detail, collapse = "; "), ")."))
  }

  .lt_assert_data_frame_columns(
    provenance$input_sha256, c("role", "path", "sha256"), "provenance$input_sha256"
  )
  .lt_assert_data_frame_columns(
    provenance$ordered_gene_ledger, c("order", "gene_id"),
    "provenance$ordered_gene_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$ordered_branch_ledger, c("order", "branch_id"),
    "provenance$ordered_branch_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$taxon_domain_ledger,
    c("order", "taxon_id", "eligible", "exclusion_reason"),
    "provenance$taxon_domain_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$specification_ledger,
    c("specification_id", "kind", "version", "definition", "sha256"),
    "provenance$specification_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$filter_ledger,
    c("order", "stage", "decision_id", "target_type", "target_id", "decision", "reason"),
    "provenance$filter_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$seed_ledger,
    c("scope", "seed", "rng_kind", "normal_kind", "sample_kind", "r_version"),
    "provenance$seed_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$fold_ledger, c("fold_id", "taxon_id", "role", "group"),
    "provenance$fold_ledger"
  )
  .lt_assert_data_frame_columns(
    provenance$output_hashes, c("role", "path", "sha256"),
    "provenance$output_hashes"
  )

  .lt_assert_nullable_scalar_character(
    provenance$configuration_sha256, "provenance$configuration_sha256"
  )
  if (!is.na(provenance$configuration_sha256)) {
    .lt_assert_sha256(
      provenance$configuration_sha256, "provenance$configuration_sha256"
    )
  }
  .lt_assert_named_list(provenance$software_version, "provenance$software_version")
  .lt_assert_named_list(provenance$environment, "provenance$environment")
  .lt_assert_named_list(provenance$eligibility_masks, "provenance$eligibility_masks")
  if (!is.character(provenance$commands)) {
    .lt_abort("`provenance$commands` must be a character vector.")
  }
  if (nrow(provenance$fold_ledger) > 0L &&
      any(!provenance$fold_ledger$role %in% c("train", "test", "excluded"))) {
    .lt_abort("`provenance$fold_ledger$role` must be train, test, or excluded.")
  }
  invisible(provenance)
}

#' Construct an audited La Terra run container
#'
#' @param run_id Stable run identifier.
#' @param matrix Optional `lt_matrix`.
#' @param trait Optional `lt_trait`.
#' @param domain Optional `lt_domain`.
#' @param analysis_spec Optional `lt_analysis_spec`.
#' @param config Parsed configuration list.
#'   An empty list means that no configuration has been attached to a draft
#'   run. Every non-empty configuration is passed through
#'   [validate_lt_config()] before the run exists.
#' @param provenance Provenance container, normally from
#'   `lt_empty_provenance()`.
#' @param status One of `draft`, `configured`, `running`, `completed`, or
#'   `failed`.
#' @param results Named result container; no scientific result is constructed.
#'
#' @return An object of class `lt_run`.
#' @export
lt_run <- function(run_id,
                   matrix = NULL,
                   trait = NULL,
                   domain = NULL,
                   analysis_spec = NULL,
                   config = list(),
                   provenance = lt_empty_provenance(),
                   status = c("draft", "configured", "running", "completed", "failed"),
                   results = list()) {
  status <- match.arg(status)
  .lt_validate_run_config(config)
  x <- structure(
    list(
      run_id = run_id,
      matrix = matrix,
      trait = trait,
      domain = domain,
      analysis_spec = analysis_spec,
      config = config,
      provenance = provenance,
      status = status,
      results = results
    ),
    class = "lt_run"
  )
  validate_lt_run(x)
  x
}

.lt_validate_run_config <- function(config) {
  .lt_assert_named_list(config, "config")
  if (!length(config)) {
    return(invisible(config))
  }
  validate_lt_config(config)
  invisible(config)
}

#' Validate an audited La Terra run container
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_run <- function(x) {
  if (!inherits(x, "lt_run")) {
    .lt_abort("`x` must inherit from `lt_run`.")
  }
  .lt_assert_scalar_character(x$run_id, "run_id")
  if (!is.null(x$matrix)) validate_lt_matrix(x$matrix)
  if (!is.null(x$trait)) validate_lt_trait(x$trait)
  if (!is.null(x$domain)) validate_lt_domain(x$domain)
  if (!is.null(x$trait) && !is.null(x$domain)) {
    validate_lt_trait_domain(x$trait, x$domain)
  }
  if (!is.null(x$analysis_spec)) validate_lt_analysis_spec(x$analysis_spec)
  .lt_validate_run_config(x$config)
  validate_lt_provenance(x$provenance)
  if (!x$status %in% c("draft", "configured", "running", "completed", "failed")) {
    .lt_abort("Unknown run status.")
  }
  .lt_assert_named_list(x$results, "results")
  invisible(x)
}
