# Negative conformance scanner for current logGBI documentation.

.lt_loggbi_forbidden_claim_patterns <- function() {
  c(
    default_production_finite_sentinel = paste0(
      "default\\s+(production\\s+)?loggbi.{0,180}",
      "(uses?|assigns?|represents?|encodes?).{0,120}",
      "finite\\s+([a-z-]+\\s+){0,3}sentinel"
    ),
    default_min_loggbi_zero = paste0(
      "zero\\s+(is\\s+)?(represented|encoded).{0,160}",
      "min\\s*\\(\\s*loggbi\\s*\\).{0,100}(by\\s+default|default)"
    ),
    generic_calibration_encoded_production = paste0(
      "generic\\s+calibration.{0,160}(consumes?|uses?).{0,160}",
      "production\\s+zero[- ]encoded\\s+loggbi"
    ),
    finite_loggbi_current_primary = paste0(
      "finite\\s+zero[- ]state\\s+loggbi.{0,160}",
      "(is|as).{0,100}(the\\s+)?current\\s+primary\\s+representation"
    ),
    current_primary_finite_loggbi = paste0(
      "current\\s+primary\\s+representation.{0,160}",
      "(is|uses?).{0,100}finite\\s+([a-z-]+\\s+){0,3}loggbi"
    ),
    production_global_k10_default = paste0(
      "production\\s+loggbi.{0,160}",
      "(min\\s*\\(\\s*log\\s*\\(\\s*gbi|global[- ]k10).{0,140}",
      "(zero|default)"
    ),
    acknowledgement_unlocks_encoded_inference = paste0(
      "acknowledg[a-z]*.{0,120}",
      "(allows?|authorizes?|enables?|unlocks?|permits?).{0,120}",
      "(encoded|finite[- ]sentinel).{0,100}inference"
    ),
    encoded_inference_requires_acknowledgement = paste0(
      "(encoded|finite[- ]sentinel).{0,100}inference.{0,120}",
      "(requires?|available.{0,30}through).{0,80}acknowledg"
    )
  )
}

.lt_loggbi_doc_files <- function(root) {
  root <- normalizePath(root, mustWork = TRUE)
  direct <- unlist(lapply(c("DESCRIPTION", "README*", "NEWS*"), function(x) {
    Sys.glob(file.path(root, x))
  }), use.names = FALSE)
  recursive_roots <- c(
    "man", "inst/doc", "vignettes", "inst/schema", "inst/recipes",
    "inst/examples", "examples", "demo"
  )
  recursive <- unlist(lapply(recursive_roots, function(x) {
    path <- file.path(root, x)
    if (!dir.exists(path)) return(character())
    list.files(path, recursive = TRUE, full.names = TRUE,
               all.files = FALSE, no.. = TRUE)
  }), use.names = FALSE)
  paths <- unique(c(direct, recursive))
  # Installed runnable examples may ship serialized input data. These are not
  # UTF-8 documentation and must not be readLines()'d as prose.
  paths <- paths[!grepl("\\.(rds|rda|rdata)$", paths, ignore.case = TRUE)]
  paths[file.exists(paths) & !file.info(paths)$isdir]
}

.lt_loggbi_doc_sections <- function(lines, path, current_version) {
  if (!grepl("^NEWS", basename(path), ignore.case = TRUE)) {
    return(list(list(
      text = paste(lines, collapse = " "),
      scope = "current",
      section = basename(path)
    )))
  }
  headers <- grep(
    "^#{1,6}\\s+LaTerra\\s+[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+",
    lines
  )
  if (!length(headers)) {
    return(list(list(
      text = paste(lines, collapse = " "),
      scope = "current",
      section = "unversioned_NEWS"
    )))
  }
  starts <- c(1L, headers)
  starts <- unique(starts)
  ends <- c(starts[-1L] - 1L, length(lines))
  Map(function(lo, hi) {
    block <- lines[lo:hi]
    version_line <- grep(
      "^#{1,6}\\s+LaTerra\\s+[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+",
      block, value = TRUE
    )
    version <- if (length(version_line)) sub(
      ".*LaTerra\\s+([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+).*",
      "\\1", version_line[[1L]]
    ) else NA_character_
    historical <- !is.na(version) &&
      utils::compareVersion(version, current_version) < 0L
    list(
      text = paste(block, collapse = " "),
      scope = if (historical) "historical_versioned_NEWS" else "current",
      section = if (is.na(version)) "NEWS_preamble" else version
    )
  }, starts, ends)
}

.lt_scan_loggbi_documentation <- function(
    root, current_version = as.character(utils::packageVersion("LaTerra")),
    paths = NULL) {
  .lt_assert_scalar_character(root, "root")
  .lt_assert_scalar_character(current_version, "current_version")
  if (is.null(paths)) paths <- .lt_loggbi_doc_files(root)
  if (!is.character(paths) || anyNA(paths)) {
    .lt_abort("`paths` must be a character vector of documentation files.")
  }
  patterns <- .lt_loggbi_forbidden_claim_patterns()
  findings <- list()
  k <- 0L
  for (path in paths) {
    if (!file.exists(path) || file.info(path)$isdir) next
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    sections <- .lt_loggbi_doc_sections(lines, path, current_version)
    for (section in sections) {
      if (!identical(section$scope, "current")) next
      normalized <- tolower(gsub("[[:space:]]+", " ", section$text))
      for (pattern_id in names(patterns)) {
        hit <- regexpr(patterns[[pattern_id]], normalized, perl = TRUE)
        if (hit[[1L]] < 0L) next
        k <- k + 1L
        start <- max(1L, hit[[1L]] - 80L)
        stop <- min(nchar(normalized),
                    hit[[1L]] + attr(hit, "match.length") + 80L)
        findings[[k]] <- data.frame(
          path = normalizePath(path, mustWork = TRUE),
          section = section$section,
          pattern_id = pattern_id,
          excerpt = substr(normalized, start, stop),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(findings)) {
    return(data.frame(
      path = character(), section = character(), pattern_id = character(),
      excerpt = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, findings)
}

.lt_assert_loggbi_documentation_conforms <- function(
    root, current_version = as.character(utils::packageVersion("LaTerra")),
    paths = NULL) {
  findings <- .lt_scan_loggbi_documentation(
    root, current_version = current_version, paths = paths
  )
  if (nrow(findings)) {
    .lt_abort_class(
      paste0(
        "Current documentation reintroduces a forbidden finite-logGBI claim: ",
        findings$pattern_id[[1L]], " in ", findings$path[[1L]], "."
      ),
      "lt_error_logGBI_documentation_nonconformant",
      list(findings = findings)
    )
  }
  invisible(findings)
}
