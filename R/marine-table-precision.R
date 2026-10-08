.lt_marine_reporting_precision_clarification <- function() {
  list(id = "TARGET001-REPORTING-CLARIFICATION002",
    sha256 = "15a4eca9ee406c56a96f16b326aa83548ea665623182f6deace2fec742aff20e",
    qualifier = "presentation_rule_frozen_by_main_control_original_generator_unlocated",
    fields = "OUT_TABLE_S5/predictor_annotation!E2:G149",
    rule = 'as.numeric(sprintf("%.6g", fresh_value))', locale = "C",
    role = "final_numeric_worksheet_view_not_scientific_coefficient_source")
}

.lt_marine_table_six_significant <- function(x) {
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
    .lt_marine_abort("S5 numeric report view requires finite full-precision inputs.", "reporting_payload")
  if (Sys.getlocale("LC_NUMERIC") != "C")
    .lt_marine_abort("Set LC_NUMERIC=C for the frozen numeric report view.", "reporting_locale")
  value <- as.numeric(sprintf("%.6g", x))
  if (any(!is.finite(value)) || any(sign(value) != sign(x)) || any((value == 0) != (x == 0)))
    .lt_marine_abort("S5 numeric presentation changed finiteness, sign or zero state.", "reporting_payload")
  value
}

.lt_marine_s5_worksheet_view <- function(full_precision_annotation) {
  fields <- c("marine_coef", "aquatic_coef", "max_abs_coef")
  if (!is.data.frame(full_precision_annotation) || anyDuplicated(names(full_precision_annotation)) ||
      !all(c("gene", fields) %in% names(full_precision_annotation)))
    .lt_marine_abort("S5 view requires the complete fresh fitted-annotation table.", "reporting_payload")
  .lt_assert_unique_ids(full_precision_annotation$gene, "S5 report-view gene keys")
  if (!identical(full_precision_annotation$max_abs_coef,
      pmax(abs(full_precision_annotation$marine_coef), abs(full_precision_annotation$aquatic_coef))))
    .lt_marine_abort("S5 max_abs_coef is not the full-precision derived maximum.", "reporting_payload")
  before <- serialize(full_precision_annotation, NULL, version = 3)
  view <- full_precision_annotation
  for (field in fields) view[[field]] <- .lt_marine_table_six_significant(full_precision_annotation[[field]])
  other <- setdiff(names(view), fields)
  if (!identical(serialize(full_precision_annotation, NULL, version = 3), before) ||
      !identical(view[other], full_precision_annotation[other]))
    .lt_marine_abort("Final worksheet view changed non-presentation fields.", "reporting_payload")
  list(values = view, provenance = list(authority = .lt_marine_reporting_precision_clarification(),
    full_precision_derived_annotation_sha256 = digest::digest(before, serialize = FALSE, algo = "sha256"),
    original_annotation_unchanged = TRUE, all_other_derived_fields_unchanged = TRUE,
    downstream_derivation_from_formatted_columns = FALSE))
}
