# Frozen Marine/S2 input adapter. This is not an upstream coordinate caller.

.lt_marine_abort <- function(message, code = "input") {
  stop(structure(list(message = message, call = NULL),
                 class = c(paste0("lt_error_marine_", code), "error", "condition")))
}

.lt_marine_sha <- function(path) digest::digest(file = path, algo = "sha256")

.lt_marine_verify <- function(path, expected, role) {
  if (!file.exists(path) || dir.exists(path)) {
    .lt_marine_abort(paste0("Missing ", role, ": ", path,
      ". Supply the frozen Dryad inputs and TARGET001-AUTHORITY-ADDENDUM001 supplement; no historical-path fallback is used."),
      "missing_input")
  }
  actual <- .lt_marine_sha(path)
  if (!identical(actual, expected)) {
    .lt_marine_abort(paste0("Frozen snapshot mismatch for ", role, ": ", path,
      ". Expected ", expected, "; found ", actual,
      ". This profile requires the frozen S2 input, not a replacement dataset."),
      "authority_mismatch")
  }
  data.frame(role = role, path = path, sha256 = actual, stringsAsFactors = FALSE)
}

.lt_marine_read_tsv <- function(path, character_only = FALSE) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                   na.strings = if (character_only) character(0) else "NA",
                   colClasses = if (character_only) "character" else NA,
                   quote = "\"", comment.char = "")
}

.lt_marine_token_levels <- function() {
  c("finite_numeric", "NA", "NA_fuse", "NA_struct", "NA_topo", "residual_NA",
    "NaN", "Inf", "-Inf", "empty")
}

.lt_marine_parse_matrix <- function(path) {
  dat <- .lt_marine_read_tsv(path, TRUE)
  if (ncol(dat) < 2L || nrow(dat) < 1L) {
    .lt_marine_abort("Matrix must contain a gene column and branch columns.", "matrix")
  }
  genes <- dat[[1L]]; branches <- names(dat)[-1L]
  .lt_assert_unique_ids(genes, "Marine gene IDs")
  .lt_assert_unique_ids(branches, "Marine branch IDs")
  tokens <- as.matrix(dat[-1L]); dimnames(tokens) <- list(genes, branches)
  levels <- .lt_marine_token_levels()
  missing_labels <- levels[-1L]
  token_names <- tokens
  token_names[token_names == ""] <- "empty"
  codes <- match(token_names, missing_labels)
  missing <- !is.na(codes)
  num <- suppressWarnings(as.numeric(tokens))
  bad <- !missing & !is.finite(num)
  if (any(bad)) {
    .lt_marine_abort(paste0("Unexpected matrix token(s): ",
      paste(utils::head(unique(tokens[bad]), 5L), collapse = ", ")), "matrix_token")
  }
  num[missing] <- NA_real_
  values <- matrix(num, nrow(tokens), ncol(tokens), dimnames = dimnames(tokens))
  codes[!missing] <- 0L
  state <- matrix(as.raw(codes), nrow(tokens), ncol(tokens),
                  dimnames = dimnames(tokens))
  list(values = values, token_code = state, token_vocabulary = levels,
       gene_ids = genes, branch_ids = branches)
}

.lt_marine_parse_complete_state <- function(path) {
  dat <- .lt_marine_read_tsv(path, TRUE)
  if (ncol(dat) < 2L || nrow(dat) < 1L) {
    .lt_marine_abort("Complete state ledger needs gene and branch axes.", "state")
  }
  genes <- dat[[1L]]; branches <- names(dat)[-1L]
  .lt_assert_unique_ids(genes, "complete-state gene IDs")
  .lt_assert_unique_ids(branches, "complete-state branch IDs")
  # This vocabulary is the admitted fixed-coordinate artifact's vocabulary,
  # not an attempted redefinition of all La Terra missingness states.
  vocabulary <- c("observed", "NA_struct", "NA_fuse")
  code <- match(as.matrix(dat[-1L]), vocabulary) - 1L
  if (anyNA(code)) {
    .lt_marine_abort("Unexpected label in the frozen complete-state ledger.", "state_token")
  }
  list(gene_ids = genes, branch_ids = branches, vocabulary = vocabulary,
    state_code = matrix(as.raw(code), length(genes), length(branches),
      dimnames = list(genes, branches)))
}

.lt_marine_state_layers <- function(raw, complete, historical) {
  if (!identical(dimnames(raw$values), dimnames(complete$state_code))) {
    .lt_marine_abort("Raw and complete supplemental state have different ordered axes.", "axis")
  }
  if (!identical(historical$branch_ids, raw$branch_ids) ||
      any(!historical$gene_ids %in% raw$gene_ids)) {
    .lt_marine_abort("Historical partial provenance is not a keyed subset of the full axis.", "axis")
  }
  if (any((complete$state_code == as.raw(0L)) != is.finite(raw$values))) {
    .lt_marine_abort("Complete upstream states conflict with raw payload availability.", "state_payload")
  }
  # No numerical values are filled, masked or replaced here. The common-cell
  # historical reconciliation is frozen addendum evidence, not re-inferred.
  list(gene_ids = raw$gene_ids, branch_ids = raw$branch_ids,
    vocabulary = raw$token_vocabulary, raw_token_code = raw$token_code,
    complete_state = complete,
    historical_partial = list(gene_ids = historical$gene_ids,
      branch_ids = historical$branch_ids, token_code = historical$token_code,
      vocabulary = historical$token_vocabulary),
    interpretation = "Raw numerical authority; historical partial provenance; supplemental complete upstream state provenance. No state regeneration or numerical modification.")
}

.lt_marine_environment <- function(profile, mode) {
  packages <- c("ape", "castor")
  for (p in packages) {
    if (!requireNamespace(p, quietly = TRUE)) {
      .lt_marine_abort(paste0("Marine annotation requires package '", p,
        "'. Install it before running; no alternate ASR is used."), "dependency")
    }
  }
  current <- list(R = as.character(getRversion()), platform = R.version$platform)
  for (p in packages) current[[p]] <- utils::packageDescription(p)$Version
  matches <- vapply(names(current), function(k)
    identical(current[[k]], profile$reference_environment[[k]]), logical(1L))
  if (mode == "frozen" && !all(matches)) {
    .lt_marine_abort(paste0("Frozen reference environment differs at: ",
      paste(names(matches)[!matches], collapse = ", "),
      ". Restore the recorded R/package stack or explicitly request ",
      "environment='compatibility' (not official certification)."), "environment")
  }
  if (mode == "compatibility") {
    warning("Explicit compatibility execution: this run does not certify the official frozen environment.",
            call. = FALSE)
  }
  list(mode = mode, reference_match = all(matches), versions = current,
       session_info = utils::capture.output(utils::sessionInfo()),
       platform_and_package_match = matches)
}

.lt_marine_inputs <- function(data_root, profile, stage) {
  rows <- list(); input <- list()
  for (role in names(profile$inputs)) {
    spec <- profile$inputs[[role]]
    path <- file.path(data_root, spec$path)
    rows[[role]] <- .lt_marine_verify(path, spec$sha256, role)
    input[[role]] <- path
  }
  members <- utils::untar(input$branch_archive, list = TRUE)
  wanted <- vapply(profile$archive_members, `[[`, character(1L), "path")
  if (anyDuplicated(members) || any(!wanted %in% members) ||
      any(grepl("(^/|^[A-Za-z]:|\\\\|(^|/)\\.\\.(/|$))", members))) {
    .lt_marine_abort("Invalid or ambiguous branch archive member layout.", "archive")
  }
  # The complete archive SHA was verified before extraction. Only named members
  # are extracted; neither GBI oracles nor output tables are execution inputs.
  utils::untar(input$branch_archive, files = wanted, exdir = stage)
  for (role in names(profile$archive_members)) {
    spec <- profile$archive_members[[role]]
    path <- file.path(stage, spec$path)
    rows[[paste0("member_", role)]] <- .lt_marine_verify(path, spec$sha256, role)
    input[[role]] <- path
  }
  input$hash_ledger <- do.call(rbind, rows)
  input
}
