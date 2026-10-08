.lt_marine_reporting_run_reader <- function(run_dir) {
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  manifest_path <- file.path(run_dir, "SHA256SUMS")
  if (!file.exists(manifest_path))
    .lt_marine_abort("Complete the Marine M2 run before exporting its tables.", "reporting_dependency")
  lines <- readLines(manifest_path)
  if (!length(lines) || any(!grepl("^[0-9a-f]{64}  .+", lines)))
    .lt_marine_abort("Malformed run output manifest.", "reporting_dependency")
  paths <- substring(lines, 67L)
  if (anyDuplicated(paths) || any(startsWith(paths, "/")) || any(grepl("(^|/)\\.\\.(/|$)", paths)))
    .lt_marine_abort("Ambiguous or escaping run output paths.", "reporting_dependency")
  if (!setequal(list.files(run_dir, recursive = TRUE), c(paths, "SHA256SUMS")))
    .lt_marine_abort("Scientific run output inventory differs from its frozen manifest.", "reporting_dependency")
  read <- function(name, kind = c("rds", "tsv")) {
    kind <- match.arg(kind)
    i <- match(name, paths)
    if (is.na(i)) .lt_marine_abort(paste0("Missing declared this-run output: ", name), "reporting_dependency")
    p <- normalizePath(file.path(run_dir, name), mustWork = TRUE)
    if (!startsWith(p, paste0(run_dir, "/")))
      .lt_marine_abort("This-run reporting dependency escapes its output directory.", "reporting_dependency")
    .lt_marine_verify(p, substr(lines[i], 1L, 64L), paste0("this_run_output_", name))
    if (kind == "rds") readRDS(p) else
      utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  }
  run <- read("run.rds"); validate_lt_run(run)
  if (!identical(run$status, "completed") || !identical(run$results$profile, "marine_NC_submitted") ||
      !identical(run$results$milestone, "M2"))
    .lt_marine_abort("Tables require a completed explicit Marine M2 run, not a substituted generic object.", "reporting_dependency")
  authority <- read("M2_authority.rds")
  if (!identical(authority, .lt_marine_m2_authority()))
    .lt_marine_abort("Run scientific-input authority differs from the current frozen Marine contract.", "reporting_dependency")
  list(read = read, run = run, manifest_sha256 = .lt_marine_sha(manifest_path))
}

.lt_marine_reporting_write_xlsx <- function(sheets, path) {
  .lt_marine_table_file_dependencies()
  wb <- openxlsx::createWorkbook(creator = "LaTerra")
  header <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#F2F2F2", wrapText = TRUE,
    valign = "top")
  body <- openxlsx::createStyle(wrapText = TRUE, valign = "top")
  for (x in sheets) {
    sheet <- x$schema$sheet; openxlsx::addWorksheet(wb, sheet)
    v <- x$values; headers <- unlist(x$schema$headers, use.names = FALSE)
    openxlsx::writeData(wb, sheet, t(headers), colNames = FALSE, rowNames = FALSE)
    for (j in seq_len(ncol(v))) {
      column <- v[-1L, j]
      live <- column[!vapply(column, is.null, TRUE)]
      types <- unique(vapply(live, function(x) if (is.numeric(x)) "numeric" else typeof(x), ""))
      if (length(types) > 1L) {
        # Preserve rare genuinely mixed typed columns, not formatted strings.
        for (r in seq_along(column)) if (!is.null(column[[r]]))
          openxlsx::writeData(wb, sheet, column[[r]], startRow = r + 1L, startCol = j, colNames = FALSE)
      } else {
        blank <- if (!length(types) || types == "numeric") NA_real_ else
          if (types == "logical") NA else NA_character_
        col <- unlist(lapply(column, function(x) if (is.null(x)) blank else x), use.names = FALSE)
        openxlsx::writeData(wb, sheet, col, startRow = 2L, startCol = j, colNames = FALSE, keepNA = FALSE)
      }
    }
    openxlsx::addStyle(wb, sheet, header, rows = 1L, cols = seq_len(ncol(v)), gridExpand = TRUE)
    openxlsx::addStyle(wb, sheet, body, rows = seq.int(2L, nrow(v)), cols = seq_len(ncol(v)), gridExpand = TRUE)
    widths <- vapply(seq_len(ncol(v)), function(j) {
      text <- vapply(v[, j], function(x) if (is.null(x)) "" else as.character(x), "")
      min(64, max(14, max(nchar(text), na.rm = TRUE) + 2))
    }, 0)
    openxlsx::setColWidths(wb, sheet, seq_len(ncol(v)), widths)
    heights <- .lt_marine_table_row_heights(v, widths)
    heights[1L] <- max(42, heights[1L])
    openxlsx::setRowHeights(wb, sheet, seq_len(nrow(v)), heights)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)
  }
  tmp <- tempfile("lt_marine_table_unfinalized_", fileext = ".xlsx")
  on.exit(unlink(tmp), add = TRUE)
  openxlsx::saveWorkbook(wb, tmp, overwrite = FALSE)
  .lt_marine_finalize_table_file(tmp, path)
}

#' Export the frozen Marine reporting tables from a completed run
#'
#' A terminal reporting step: uses this-run fitted/screened outputs and the
#' explicitly admitted metadata projection. It does not refit models, recalculate
#' GBI or traits, generate random worlds, or read archived result workbooks.
#' The five original worksheet schemas are preserved. Historical reporting
#' references and later-frozen presentation rules are disclosed in companions.
#' Completion is not independent TARGET001 or scientific certification.
#'
#' @param data_root Root of the same frozen input bundle used by the M2 run.
#' @param run_dir Directory produced by [lt_run_marine_m2()].
#' @param output_dir New directory outside the input and scientific-run roots.
#' @return Invisibly a list of workbook paths and reporting provenance.
#' @export
lt_export_marine_tables <- function(data_root, run_dir, output_dir) {
  .lt_assert_scalar_character(data_root, "data_root")
  .lt_assert_scalar_character(run_dir, "run_dir")
  .lt_assert_scalar_character(output_dir, "output_dir")
  .lt_marine_table_file_dependencies()
  data_root <- normalizePath(data_root, mustWork = TRUE)
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  parent <- normalizePath(dirname(output_dir), mustWork = TRUE)
  if (file.exists(output_dir)) .lt_marine_abort("Choose a new reporting directory; existing reports are not overwritten.", "output_exists")
  if (any(vapply(c(data_root, run_dir), function(p) parent == p || startsWith(parent, paste0(p, "/")), TRUE)))
    .lt_marine_abort("Reporting outputs must be outside the input and scientific-run roots.", "output_location")
  reader <- .lt_marine_reporting_run_reader(run_dir)
  read <- reader$read
  # Only the already consumed input views needed by terminal reporting are read.
  a <- .lt_marine_m2_authority()$inputs
  read_input <- function(spec) {
    path <- file.path(data_root, spec$path)
    .lt_marine_verify(path, spec$sha256, "already_admitted_reporting_input")
    if (!spec$sha256 %in% reader$run$provenance$input_sha256$sha256)
      .lt_marine_abort("Reporting input was not consumed by this scientific run.", "reporting_dependency")
    utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  }
  traits <- read_input(lt_marine_profile()$inputs$traits)
  annotation <- read_input(a$final_modules); worlds <- read_input(a$fig5b_worlds)
  output_dir <- file.path(parent, basename(output_dir))
  if (!dir.create(output_dir)) .lt_marine_abort("Cannot reserve reporting directory.", "io")
  tryCatch({
    models <- read("computed_model_inputs.rds")
    s1 <- .lt_marine_table_s1(traits)
    screen <- function(trait, kind) read(paste0("M1/", trait, ".", kind, ".tsv"), "tsv")
    mt <- screen("marine_binary", "tested"); ms <- screen("marine_binary", "significant")
    at <- screen("aquatic_v2", "tested"); az <- screen("aquatic_v2", "significant")
    s2m <- .lt_marine_table_s2(mt, ms, "marine_binary")
    s2a <- .lt_marine_table_s2(at, az, "aquatic_v2")
    s4 <- .lt_marine_table_s4_sensitivity(models)
    overlap <- .lt_marine_table_overlap(mt, at, ms, az)
    s5 <- .lt_marine_table_s5(read("full_data.rds"), annotation)
    # All scientific joins, markers, ordering and lists are complete before
    # constructing the separately authorized final numeric worksheet view.
    s5_view <- .lt_marine_s5_worksheet_view(s5$predictor_annotation)
    s5_view$provenance$upstream_full_precision_object <- list(path = "full_data.rds",
      source_root = "scientific", sha256 = .lt_marine_sha(file.path(run_dir, "full_data.rds")))
    ids <- .lt_marine_sensitivity_ids()
    sensitivity <- stats::setNames(lapply(ids, function(id) read(paste0("nested_phase13_", id, ".rds"))), ids)
    s6 <- .lt_marine_table_s6(sensitivity, read("Fig5B.rds"), worlds, read("turnover.rds"))
    current <- list(OUT_TABLE_S1 = list(SupplementaryTableS1 = s1),
      OUT_TABLE_S2 = list(marine_significant_genes = s2m, aquatic_significant_genes = s2a),
      OUT_TABLE_S4 = c(s4, list(Sheet1 = overlap)),
      OUT_TABLE_S5 = list(predictor_annotation = s5_view$values,
        display_module_summary = s5$display_module_summary),
      OUT_TABLE_S6 = s6[setdiff(names(s6), "provenance")])
    authority <- .lt_marine_reporting_authority(); sheets <- list()
    for (schema in authority$schemas) {
      id <- schema$source_id; sheet <- schema$sheet
      d <- current[[id]][[sheet]]
      if (!is.null(d)) d <- .lt_marine_reporting_fields(id, sheet, d)
      sheets[[paste(id, sheet, sep = "/")]] <- .lt_marine_assemble_reporting_sheet(id, sheet, d)
    }
    saveRDS(sheets, file.path(output_dir, "logical_sheets.rds"), version = 3)
    output_names <- c(OUT_TABLE_S1 = "Supplementary_Table_S1_aquatic_marine_traits.xlsx",
      OUT_TABLE_S2 = "Supplementary_Table_S2_ttest_significant_genes.xlsx",
      OUT_TABLE_S4 = "Supplementary_Table_S4_ttest_drop_sensitivity.xlsx",
      OUT_TABLE_S5 = "Supplementary_Table_S5_Figure4C_predictor_annotation.xlsx",
      OUT_TABLE_S6 = "Supplementary_Table_S6_Figure5_sensitivity_permutation_turnover.xlsx")
    files <- stats::setNames(file.path(output_dir, output_names), names(output_names))
    file_finalization <- stats::setNames(lapply(names(files), function(id)
      .lt_marine_reporting_write_xlsx(sheets[startsWith(names(sheets), paste0(id, "/"))], files[id])), names(files))
    installed <- system.file(package = "LaTerra")
    runtime_files <- list.files(installed, recursive = TRUE, full.names = TRUE)
    runtime_hashes <- data.frame(path = substring(runtime_files, nchar(installed) + 2L),
      sha256 = vapply(runtime_files, .lt_marine_sha, ""), stringsAsFactors = FALSE)
    prov <- list(profile = "marine_NC_submitted", stage = "terminal_reporting_only",
      scientific_run_id = reader$run$run_id, scientific_run_manifest_sha256 = reader$manifest_sha256,
      scientific_run_input_sha256 = reader$run$provenance$input_sha256,
      scientific_source_software = reader$run$provenance$software_version,
      reporting_package_version = as.character(utils::packageVersion("LaTerra")),
      reporting_installed_runtime_files = runtime_hashes,
      reporting_environment = utils::capture.output(utils::sessionInfo()),
      command = paste0("lt_export_marine_tables(data_root=", encodeString(data_root, quote = '"'),
        ", run_dir=", encodeString(run_dir, quote = '"'), ", output_dir=", encodeString(output_dir, quote = '"'), ")"),
      reporting_authority = authority[c("authority_id", "authority_sha256", "source_sha256")],
      clarification = .lt_marine_reporting_clarification(),
      precision_clarification = .lt_marine_reporting_precision_clarification(),
      S5_final_numeric_view = s5_view$provenance,
      sheets = lapply(sheets, `[[`, "provenance"), S5 = s5$provenance, S6 = s6$provenance,
      workbook_file_finalization = file_finalization,
      S4_historical_context = attr(overlap, "historical_context_check"),
      existing_scientific_outputs_modified = FALSE, original_workbook_consumed = FALSE,
      new_scientific_calculation = FALSE, TARGET001_acceptance = "NOT_CLAIMED")
    saveRDS(prov, file.path(output_dir, "reporting_provenance.rds"), version = 3)
    writeLines(c("Marine reporting tables, assembled from the identified scientific run.",
      "No original result workbook was consumed by this reporting operation.",
      "S4 Sheet1 rows 24-29: historical_reference_imported_not_recomputed; discussion-only reference.",
      "S6 overview's old no-analysis-rerun wording describes original table compilation, not this fresh scientific replay.",
      "S5 coefficient phrases and S6 percentage display: presentation_rule_frozen_by_main_control_original_generator_unlocated.",
      "S5 E:G only: final numeric worksheet view at six significant digits; full-precision scientific coefficients/derived fields unchanged. Original generator unlocated.",
      "Historical Fig5B direction and Phase13 screening-map qualifiers remain mandatory.",
      "Main predictor-versus-predictor evidence and cross-layer descriptive context remain separate.",
      "These are computed reports, not an independent TGT021, M2/M3, scientific or release PASS."),
      file.path(output_dir, "REPORTING_QUALIFIERS.txt"))
    yaml::write_yaml(list(status = "TABLES_ASSEMBLED_NOT_REPLAY_CERTIFIED", workbook_count = 5L,
      sheet_count = length(sheets), scientific_run_manifest_sha256 = reader$manifest_sha256),
      file.path(output_dir, "REPORTING_STATUS.yaml"))
    completed_run <- .lt_marine_complete_reporting_run(reader, prov, run_dir, output_dir, output_names)
    paths <- list.files(output_dir, full.names = TRUE)
    writeLines(paste(vapply(paths, .lt_marine_sha, ""), basename(paths), sep = "  "), file.path(output_dir, "SHA256SUMS"))
    invisible(list(files = files, provenance = prov, run = completed_run))
  }, error = function(e) {
    yaml::write_yaml(list(status = "FAILED_REPORTING_PARTIAL_RETAINED", message = conditionMessage(e),
      condition_class = class(e), scientific_outputs_modified = FALSE), file.path(output_dir, "REPORTING_STATUS.yaml"))
    stop(e)
  })
}
