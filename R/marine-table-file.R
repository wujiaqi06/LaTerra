.lt_marine_table_file_dependencies <- function() {
  missing <- c("openxlsx", "xml2", "zip")[!vapply(c("openxlsx", "xml2", "zip"),
    requireNamespace, TRUE, quietly = TRUE)]
  if (length(missing)) .lt_marine_abort(paste0("Install these optional R packages to export .xlsx tables: ",
    paste(missing, collapse = ", "), ". No alternate format is silently substituted."), "dependency")
}

.lt_marine_table_row_heights <- function(values, widths) {
  if (!is.matrix(values) || length(widths) != ncol(values) || any(!is.finite(widths) | widths <= 2))
    .lt_marine_abort("Invalid worksheet layout dimensions.", "reporting_file")
  lines <- matrix(1, nrow(values), ncol(values))
  for (j in seq_len(ncol(values))) for (r in seq_len(nrow(values))) {
    x <- values[[r, j]]
    if (is.null(x)) next
    pieces <- strsplit(as.character(x), "\n", fixed = TRUE)[[1L]]
    # Conservative wrapping allowance for glyphs/spaces, without touching text.
    lines[r, j] <- max(1, sum(pmax(1, ceiling(nchar(pieces, type = "width") / ((widths[j] - 2) * .85)))))
  }
  pmin(409.5, pmax(18, 15 * apply(lines, 1L, max) + 4))
}

.lt_marine_table_part_path <- function(base, target) {
  if (!is.character(target) || length(target) != 1L || is.na(target) || !nzchar(target) ||
      grepl("[\\\\?#]", target) || grepl("^[[:alpha:]][[:alnum:]+.-]*:", target))
    .lt_marine_abort("Invalid local workbook part target.", "reporting_file")
  parts <- strsplit(if (startsWith(target, "/")) substring(target, 2L) else
    paste(base, target, sep = "/"), "/", fixed = TRUE)[[1L]]
  stack <- character()
  for (p in parts) {
    if (p %in% c("", ".")) next
    if (p == "..") {
      if (!length(stack)) .lt_marine_abort("Workbook part escapes its package.", "reporting_file")
      stack <- utils::head(stack, -1L)
    } else stack <- c(stack, p)
  }
  paste(stack, collapse = "/")
}

# Some openxlsx versions emit dangling, unused drawing/VML relationships even
# for value-only sheets. Remove only those unused links, never cell/worksheet
# content, and require every other internal package relationship to resolve.
.lt_marine_finalize_table_file <- function(source, destination) {
  .lt_marine_table_file_dependencies()
  if (file.exists(destination)) .lt_marine_abort("Existing workbook is not overwritten.", "output_exists")
  tmp <- tempfile("lt_marine_xlsx_parts_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  entries <- utils::unzip(source, list = TRUE)$Name
  if (!length(entries) || anyDuplicated(entries) || any(startsWith(entries, "/")) ||
      any(grepl("(^|/)\\.\\.(/|$)|\\\\", entries)))
    .lt_marine_abort("Unsafe or ambiguous workbook archive.", "reporting_file")
  utils::unzip(source, exdir = tmp)
  files <- list.files(tmp, recursive = TRUE, all.files = TRUE)
  hashes <- vapply(file.path(tmp, files), .lt_marine_sha, "")
  names(hashes) <- files
  removed <- list(); changed <- character()
  ns <- "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
  for (relpath in files[grepl("(^|/)_rels/.*[.]rels$", files)]) {
    owner <- if (relpath == "_rels/.rels") "" else
      file.path(dirname(dirname(relpath)), sub("[.]rels$", "", basename(relpath)))
    doc <- xml2::read_xml(file.path(tmp, relpath))
    relationships <- xml2::xml_find_all(doc, "/*[local-name()='Relationships']/*[local-name()='Relationship']")
    ids <- xml2::xml_attr(relationships, "Id")
    if (anyNA(ids) || anyDuplicated(ids)) .lt_marine_abort("Ambiguous workbook relationship IDs.", "reporting_file")
    for (i in seq_along(relationships)) {
      node <- relationships[[i]]
      if (identical(xml2::xml_attr(node, "TargetMode"), "External")) next
      target <- xml2::xml_attr(node, "Target")
      resolved <- .lt_marine_table_part_path(if (nzchar(owner)) dirname(owner) else ".", target)
      if (resolved %in% files) next
      type <- xml2::xml_attr(node, "Type")
      allowed <- type %in% paste0(ns, c("/drawing", "/vmlDrawing")) &&
        grepl("^xl/worksheets/[^/]+[.]xml$", owner) && owner %in% files
      if (!allowed) .lt_marine_abort(paste0("Missing required workbook part: ", resolved), "reporting_file")
      sheet <- xml2::read_xml(file.path(tmp, owner))
      used <- xml2::xml_text(xml2::xml_find_all(sheet,
        paste0("//@*[namespace-uri()='", ns, "']")))
      if (ids[i] %in% used)
        .lt_marine_abort("A referenced worksheet drawing is missing; it cannot be silently removed.", "reporting_file")
      removed[[length(removed) + 1L]] <- list(relationship_part = relpath, id = ids[i], target = target,
        reason = "unused_dangling_drawing_relationship_only")
      xml2::xml_remove(node); changed <- union(changed, relpath)
    }
    if (relpath %in% changed) xml2::write_xml(doc, file.path(tmp, relpath))
  }
  # Explicit byte-identity proof for all non-relationship parts, including values,
  # shared strings, styles and source worksheet XML; no cell edit is possible.
  unchanged <- setdiff(files, changed)
  if (!identical(unname(hashes[unchanged]), unname(vapply(file.path(tmp, unchanged), .lt_marine_sha, ""))))
    .lt_marine_abort("Workbook finalization changed a non-relationship part.", "reporting_file")
  zip::zipr(destination, files = files, root = tmp, mode = "mirror")
  list(removed_unused_relationships = removed, non_relationship_parts_byte_identical = TRUE,
    cell_values_modified = FALSE)
}
