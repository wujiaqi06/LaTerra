#' @export
print.lt_matrix <- function(x, ...) {
  validate_lt_matrix(x)
  counts <- table(factor(x$coord_state, levels = lt_allowed_states()))
  cat("<lt_matrix>\n")
  cat("  dimensions: ", nrow(x$values), " genes x ", ncol(x$values), " branches\n", sep = "")
  cat("  payload: ", x$payload$type, "/", x$payload$name, " (", x$payload$origin, ")\n", sep = "")
  cat("  coordinate states: ", paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", "), "\n", sep = "")
  cat("  value reasons: ", sum(!is.na(x$value_reason)), " cells\n", sep = "")
  cat("  coordinate system: ", x$coordinate_provenance$coordinate_system, "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_trait <- function(x, ...) {
  validate_lt_trait(x)
  cat("<lt_trait> ", x$trait_id, "\n", sep = "")
  cat("  type: ", x$type, "\n", sep = "")
  cat("  taxa: ", length(x$taxon_ids), "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_domain <- function(x, ...) {
  validate_lt_domain(x)
  cat("<lt_domain> ", x$domain_id, "\n", sep = "")
  cat("  taxa: ", length(x$taxon_ids), "\n", sep = "")
  cat("  eligible: ", sum(x$eligible), "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_rate_spec <- function(x, ...) {
  validate_lt_rate_spec(x)
  cat("<lt_rate_spec> ", x$rate_id, "\n", sep = "")
  cat("  method: ", x$method, "\n", sep = "")
  cat("  trait-independent: ", x$trait_independent, "\n", sep = "")
  invisible(x)
}

#' @export
print.lt_analysis_spec <- function(x, ...) {
  validate_lt_analysis_spec(x)
  cat("<lt_analysis_spec> ", x$analysis_id, "\n", sep = "")
  cat("  rate: ", x$rate$rate_id, "\n", sep = "")
  cat("  scientific execution: not implemented\n")
  invisible(x)
}

#' @export
print.lt_run <- function(x, ...) {
  validate_lt_run(x)
  cat("<lt_run> ", x$run_id, "\n", sep = "")
  cat("  status: ", x$status, "\n", sep = "")
  cat("  provenance fields: ", length(x$provenance), "\n", sep = "")
  invisible(x)
}

#' Plot a La Terra matrix
#'
#' Plot semantics are intentionally not frozen in LT001-A1A.
#'
#' @param x An `lt_matrix`.
#' @param ... Reserved for a future plot specification.
#' @return No value; always errors in this skeleton release.
#' @export
plot.lt_matrix <- function(x, ...) {
  validate_lt_matrix(x)
  .lt_not_implemented("lt_matrix plot")
}
#' @export
print.lt_branch_state <- function(x, ...) {
  validate_lt_branch_state(x)
  counts <- table(factor(
    x$availability, levels = .lt_branch_state_availability_levels()
  ))
  cat("<lt_branch_state>", x$state_id, "\n")
  cat("  trait:", x$trait_id, "\n")
  cat("  type:", x$type, "\n")
  cat("  branches:", length(x$branch_ids), "\n")
  cat("  availability:", paste(names(counts), counts, sep = "=",
                                collapse = ", "), "\n")
  invisible(x)
}

#' @export
print.lt_association_domain <- function(x, ...) {
  validate_lt_association_domain(x)
  cat("<lt_association_domain>", x$domain_id, "\n")
  cat("  mode:", x$domain_mode, "\n")
  cat("  representations:",
      paste(x$rate_identities$representation, collapse = ", "), "\n")
  cat("  genes x branches:", length(x$gene_ids), "x",
      length(x$branch_ids), "\n")
  cat("  eligible cells:", sum(x$eligible), "\n")
  invisible(x)
}

#' @export
print.lt_association <- function(x, ...) {
  validate_lt_association(x)
  cat("<lt_association>", x$association_id, "\n")
  cat("  representation:", x$rate_representation, "\n")
  cat("  branch state:", x$branch_state_id, "(", x$trait_type, ")\n")
  cat("  estimable genes:", sum(x$status == "ESTIMABLE"), "/",
      length(x$gene_ids), "\n")
  cat("  inference:", x$inference_status, "\n")
  invisible(x)
}
