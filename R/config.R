.lt_config_required_sections <- function() {
  c(
    "schema_version", "recipe_id", "execution", "matrix", "trait", "domain",
    "rate", "branch_state", "screen", "preprocess", "validation", "model",
    "outputs"
  )
}

#' Read or write a La Terra YAML configuration
#'
#' These helpers parse and serialize declarative configuration only. They do
#' not execute a recipe.
#'
#' @param path YAML file path.
#' @param validate Whether to apply the authoritative R-native structural
#'   validator.
#' @param config Named configuration list.
#'
#' @return `lt_read_config()` returns a named list. `lt_write_config()` returns
#'   the normalized path invisibly.
#' @export
lt_read_config <- function(path, validate = TRUE) {
  .lt_assert_scalar_character(path, "path")
  if (!file.exists(path)) {
    .lt_abort(sprintf("Configuration file does not exist: %s", path))
  }
  config <- yaml::read_yaml(path)
  if (validate) validate_lt_config(config)
  config
}

#' @rdname lt_read_config
#' @export
lt_write_config <- function(config, path, validate = TRUE) {
  .lt_assert_scalar_character(path, "path")
  if (validate) validate_lt_config(config)
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  yaml::write_yaml(config, path)
  invisible(normalizePath(path, mustWork = TRUE))
}

#' Validate a La Terra configuration container
#'
#' This R-native semantic validator is authoritative for V1. The YAML/JSON
#' Schema under `inst/schema/` is the machine-readable interoperability
#' contract and does not add a runtime JSON-Schema dependency.
#'
#' @param config Named configuration list.
#' @return `config`, invisibly.
#' @export
validate_lt_config <- function(config) {
  .lt_assert_named_list(config, "config")
  missing <- setdiff(.lt_config_required_sections(), names(config))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("Configuration is missing sections: %s.", paste(missing, collapse = ", ")))
  }
  .lt_assert_scalar_character(as.character(config$schema_version), "schema_version")
  if (!identical(as.character(config$schema_version), "1.1.0")) {
    .lt_abort("`schema_version` must be `1.1.0` for the LT001-A1A contract.")
  }
  .lt_assert_scalar_character(config$recipe_id, "recipe_id")

  list_sections <- setdiff(
    .lt_config_required_sections(),
    c("schema_version", "recipe_id")
  )
  for (name in list_sections) {
    .lt_assert_named_list(config[[name]], paste0("config$", name))
  }
  null_modules <- config$null_model
  if (is.null(null_modules)) {
    null_modules <- list()
  }
  if (!is.list(null_modules)) {
    .lt_abort("`config$null_model` must be a list of method modules.")
  }

  if (!is.logical(config$execution$enabled) ||
      length(config$execution$enabled) != 1L || is.na(config$execution$enabled)) {
    .lt_abort("`config$execution$enabled` must be TRUE or FALSE.")
  }
  if (!is.logical(config$execution$scientific_algorithms_implemented) ||
      length(config$execution$scientific_algorithms_implemented) != 1L ||
      is.na(config$execution$scientific_algorithms_implemented)) {
    .lt_abort("`config$execution$scientific_algorithms_implemented` must be TRUE or FALSE.")
  }

  matrix_required <- c(
    "input", "payload", "coord_state", "value_reason", "coordinate_provenance"
  )
  matrix_missing <- setdiff(matrix_required, names(config$matrix))
  if (length(matrix_missing) > 0L) {
    .lt_abort(sprintf("Matrix configuration is missing: %s.", paste(matrix_missing, collapse = ", ")))
  }
  .lt_validate_matrix_input(config$matrix$input, "config$matrix$input")
  .lt_validate_config_payload(config$matrix$payload)

  coord_state <- config$matrix$coord_state
  .lt_assert_named_list(coord_state, "config$matrix$coord_state")
  .lt_assert_scalar_character(coord_state$mode, "config$matrix$coord_state$mode")
  if (!coord_state$mode %in% c("all_observed", "separate_input")) {
    .lt_abort("Coordinate-state mode must be `all_observed` or `separate_input`.")
  }
  if (identical(coord_state$mode, "separate_input")) {
    if (is.null(coord_state$input)) {
      .lt_abort("Missing-coordinate input requires an explicit coordinate-state source.")
    }
    .lt_validate_file_input(coord_state$input, "config$matrix$coord_state$input")
  } else if (!is.null(coord_state$input)) {
    .lt_abort("All-observed coordinate-state mode must not declare a separate source.")
  }
  allowed <- coord_state$allowed_states
  if (!is.character(allowed) || length(allowed) != length(lt_allowed_states()) ||
      anyDuplicated(allowed) || !setequal(allowed, lt_allowed_states())) {
    .lt_abort("Matrix configuration must declare the complete frozen state vocabulary.")
  }
  if (!is.logical(coord_state$numeric_zero_is_observed) ||
      length(coord_state$numeric_zero_is_observed) != 1L ||
      !isTRUE(coord_state$numeric_zero_is_observed)) {
    .lt_abort("Matrix configuration must declare numeric zero as observed.")
  }

  value_reason <- config$matrix$value_reason
  .lt_assert_named_list(value_reason, "config$matrix$value_reason")
  .lt_assert_scalar_character(value_reason$mode, "config$matrix$value_reason$mode")
  if (!value_reason$mode %in% c("none", "separate_input")) {
    .lt_abort("Value-reason mode must be `none` or `separate_input`.")
  }
  if (identical(value_reason$mode, "separate_input")) {
    if (is.null(value_reason$input)) {
      .lt_abort("Separate value reasons require an explicit input source.")
    }
    .lt_validate_file_input(value_reason$input, "config$matrix$value_reason$input")
  } else if (!is.null(value_reason$input)) {
    .lt_abort("Value-reason mode `none` must not declare a separate source.")
  }

  .lt_validate_coordinate_provenance(config$matrix$coordinate_provenance)

  module_names <- c("rate", "branch_state", "screen", "preprocess", "validation", "model")
  for (name in module_names) {
    .lt_validate_method_module(config[[name]], paste0("config$", name))
  }
  if (identical(config$rate$method, "logGBI")) {
    .lt_abort_class(
      paste(
        "Configuration method `logGBI` is a historical finite encoding.",
        "Configure GBI and construct `lt_log_relative(gbi)` explicitly;",
        "use `lt_migrate_logGBI(...)` for historical objects, or",
        "legacy_global_k10 only for diagnostic or source-bound replay validation."
      ),
      "lt_error_logGBI_constructor_moved",
      list(requested_method = "logGBI", source = "configuration")
    )
  }
  for (i in seq_along(null_modules)) {
    .lt_validate_method_module(
      null_modules[[i]],
      paste0("config$null_model[[", i, "]]"),
      allow_disabled = TRUE
    )
  }

  outputs_required <- c("root", "formats", "required_provenance")
  outputs_missing <- setdiff(outputs_required, names(config$outputs))
  if (length(outputs_missing) > 0L) {
    .lt_abort(sprintf("Outputs configuration is missing: %s.", paste(outputs_missing, collapse = ", ")))
  }
  if (!setequal(config$outputs$required_provenance, names(lt_empty_provenance()))) {
    .lt_abort("Outputs must declare the complete frozen provenance contract.")
  }

  if (!is.null(config$authority)) {
    .lt_validate_authority_set(config$authority)
  }

  invisible(config)
}

.lt_validate_file_input <- function(input, name) {
  .lt_assert_named_list(input, name)
  missing <- setdiff(c("path", "sha256", "format"), names(input))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("`%s` is missing: %s.", name, paste(missing, collapse = ", ")))
  }
  .lt_assert_scalar_character(input$path, paste0(name, "$path"))
  .lt_assert_scalar_character(input$format, paste0(name, "$format"))
  if (!input$format %in% c("tsv", "csv", "rds", "archive_member")) {
    .lt_abort(sprintf("`%s$format` must be tsv, csv, rds, or archive_member.", name))
  }
  if (!is.null(input$sha256)) {
    .lt_assert_sha256(input$sha256, paste0(name, "$sha256"))
  }
  if (identical(input$format, "archive_member")) {
    .lt_assert_scalar_character(input$member, paste0(name, "$member"))
    if (!is.null(input$member_sha256)) {
      .lt_assert_sha256(input$member_sha256, paste0(name, "$member_sha256"))
    }
  } else if (!is.null(input$member) || !is.null(input$member_sha256)) {
    .lt_abort(sprintf("`%s` may declare archive-member fields only for archive_member format.", name))
  }
  invisible(input)
}

.lt_validate_matrix_input <- function(input, name) {
  .lt_validate_file_input(input, name)
  missing <- setdiff(c("gene_axis", "branch_axis", "gene_id_column"), names(input))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("`%s` is missing: %s.", name, paste(missing, collapse = ", ")))
  }
  if (!identical(input$gene_axis, "rows") || !identical(input$branch_axis, "columns")) {
    .lt_abort("V1 matrix input must declare gene rows and branch columns.")
  }
  .lt_assert_scalar_character(input$gene_id_column, paste0(name, "$gene_id_column"))
  invisible(input)
}

.lt_validate_config_payload <- function(payload) {
  .lt_assert_named_list(payload, "config$matrix$payload")
  required <- c("type", "name", "units", "scale", "origin", "spec_id")
  missing <- setdiff(required, names(payload))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("Matrix payload is missing: %s.", paste(missing, collapse = ", ")))
  }
  for (field in c("type", "name", "units", "scale", "origin")) {
    .lt_assert_scalar_character(payload[[field]], paste0("config$matrix$payload$", field))
  }
  if (!payload$origin %in% c("imported", "derived")) {
    .lt_abort("Matrix payload origin must be `imported` or `derived`.")
  }
  if (!is.null(payload$spec_id)) {
    .lt_assert_scalar_character(payload$spec_id, "config$matrix$payload$spec_id")
  }
  if (identical(payload$origin, "derived") && is.null(payload$spec_id)) {
    .lt_abort("A derived matrix payload requires `spec_id`.")
  }
  invisible(payload)
}

.lt_validate_coordinate_provenance <- function(provenance) {
  .lt_assert_named_list(provenance, "config$matrix$coordinate_provenance")
  required <- c(
    "coordinate_system", "reference_tree_sha256",
    "ordered_branch_ledger_sha256", "ordered_taxon_ledger_sha256",
    "branch_label_contract", "source_method"
  )
  missing <- setdiff(required, names(provenance))
  if (length(missing) > 0L) {
    .lt_abort(sprintf("Coordinate provenance is missing: %s.", paste(missing, collapse = ", ")))
  }
  for (field in c("coordinate_system", "branch_label_contract", "source_method")) {
    .lt_assert_scalar_character(
      provenance[[field]], paste0("config$matrix$coordinate_provenance$", field)
    )
  }
  for (field in c(
    "reference_tree_sha256", "ordered_branch_ledger_sha256",
    "ordered_taxon_ledger_sha256"
  )) {
    .lt_assert_sha256(
      provenance[[field]], paste0("config$matrix$coordinate_provenance$", field)
    )
  }
  for (field in intersect(
    c(
      "reference_tree_path", "ledger_hash_canonicalization",
      "source_software", "source_version", "source_commit"
    ),
    names(provenance)
  )) {
    .lt_assert_scalar_character(
      provenance[[field]], paste0("config$matrix$coordinate_provenance$", field)
    )
  }
  if (!is.null(provenance$mapping_artifact)) {
    .lt_validate_file_input(
      provenance$mapping_artifact,
      "config$matrix$coordinate_provenance$mapping_artifact"
    )
  }
  for (field in intersect(c("already_mapped", "recompute_coordinates"), names(provenance))) {
    value <- provenance[[field]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      .lt_abort(sprintf(
        "`config$matrix$coordinate_provenance$%s` must be TRUE or FALSE.", field
      ))
    }
  }
  if (!is.null(provenance$recompute_coordinates) && provenance$recompute_coordinates) {
    .lt_abort("La Terra configuration must not request coordinate recomputation.")
  }
  invisible(provenance)
}

.lt_validate_authority_set <- function(authority) {
  .lt_assert_named_list(authority, "config$authority")
  .lt_assert_scalar_character(
    authority$contract_version, "config$authority$contract_version"
  )
  if (!is.list(authority$records) || length(authority$records) < 1L) {
    .lt_abort("`config$authority$records` must contain at least one authority record.")
  }
  allowed_roles <- c(
    "scientific_specification", "executable_reference", "data_result",
    "source_code_identity"
  )
  roles <- character(length(authority$records))
  for (i in seq_along(authority$records)) {
    record <- authority$records[[i]]
    name <- paste0("config$authority$records[[", i, "]]")
    .lt_assert_named_list(record, name)
    missing <- setdiff(c("role", "identity", "governs"), names(record))
    if (length(missing) > 0L) {
      .lt_abort(sprintf("`%s` is missing: %s.", name, paste(missing, collapse = ", ")))
    }
    for (field in c("role", "identity", "governs")) {
      .lt_assert_scalar_character(record[[field]], paste0(name, "$", field))
    }
    if (!record$role %in% allowed_roles) {
      .lt_abort(sprintf("`%s$role` is not a recognized authority role.", name))
    }
    roles[[i]] <- record$role
    if (!is.null(record$sha256)) {
      .lt_assert_sha256(record$sha256, paste0(name, "$sha256"))
    }
    for (field in intersect(c("tag", "tag_object", "peeled_commit"), names(record))) {
      .lt_assert_scalar_character(record[[field]], paste0(name, "$", field))
    }
    for (field in intersect(c("tag_object", "peeled_commit"), names(record))) {
      if (!grepl("^[0-9a-f]{40}$", record[[field]])) {
        .lt_abort(sprintf("`%s$%s` must be a 40-character Git object ID.", name, field))
      }
    }
  }
  if (anyDuplicated(roles)) {
    .lt_abort("Authority roles must be unique within one authority set.")
  }
  invisible(authority)
}

.lt_validate_method_module <- function(module, name, allow_disabled = TRUE) {
  .lt_assert_named_list(module, name)
  required <- c("enabled", "method", "contract_version", "parameters")
  missing <- setdiff(required, names(module))
  unknown <- setdiff(names(module), required)
  if (length(missing) > 0L || length(unknown) > 0L) {
    .lt_abort(sprintf("`%s` must contain exactly: %s.", name, paste(required, collapse = ", ")))
  }
  if (!is.logical(module$enabled) || length(module$enabled) != 1L || is.na(module$enabled)) {
    .lt_abort(sprintf("`%s$enabled` must be TRUE or FALSE.", name))
  }
  if (!allow_disabled && !module$enabled) {
    .lt_abort(sprintf("`%s` must be enabled.", name))
  }
  .lt_assert_scalar_character(module$method, paste0(name, "$method"))
  .lt_assert_scalar_character(module$contract_version, paste0(name, "$contract_version"))
  .lt_assert_named_list(module$parameters, paste0(name, "$parameters"))
  invisible(module)
}
