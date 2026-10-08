#' Construct a La Terra analysis specification
#'
#' `lt_analysis_spec()` is a declarative container. It does not execute any
#' scientific operation.
#'
#' @param analysis_id Stable analysis identifier.
#' @param rate An `lt_rate_spec`.
#' @param branch_state,screen,preprocess,validation,model,outputs Named
#'   configuration lists for the corresponding workflow layers.
#' @param null_model A list of zero or more named null-model module lists.
#' @param metadata Optional named metadata list.
#'
#' @return An object of class `lt_analysis_spec`.
#' @export
lt_analysis_spec <- function(analysis_id,
                             rate,
                             branch_state = list(),
                             screen = list(),
                             preprocess = list(),
                             validation = list(),
                             model = list(),
                             null_model = list(),
                             outputs = list(),
                             metadata = list()) {
  x <- structure(
    list(
      analysis_id = analysis_id,
      rate = rate,
      branch_state = branch_state,
      screen = screen,
      preprocess = preprocess,
      validation = validation,
      model = model,
      null_model = null_model,
      outputs = outputs,
      metadata = metadata
    ),
    class = "lt_analysis_spec"
  )
  validate_lt_analysis_spec(x)
  x
}

#' Validate a La Terra analysis specification
#'
#' @param x Object to validate.
#' @return `x`, invisibly.
#' @export
validate_lt_analysis_spec <- function(x) {
  if (!inherits(x, "lt_analysis_spec")) {
    .lt_abort("`x` must inherit from `lt_analysis_spec`.")
  }
  .lt_assert_scalar_character(x$analysis_id, "analysis_id")
  validate_lt_rate_spec(x$rate)
  component_names <- c(
    "branch_state", "screen", "preprocess", "validation", "model",
    "outputs", "metadata"
  )
  for (name in component_names) {
    .lt_validate_component_list(x[[name]], name)
  }
  if (!is.list(x$null_model)) {
    .lt_abort("`null_model` must be a list of module configurations.")
  }
  for (i in seq_along(x$null_model)) {
    .lt_validate_component_list(
      x$null_model[[i]], paste0("null_model[[", i, "]]"))
  }
  invisible(x)
}
