# File convenience layer. Scientific operators remain in lt_run_s2().
.lt_binary_tsv <- function(file, required) {
  .lt_assert_scalar_character(file, "file")
  if (!file.exists(file) || dir.exists(file)) stop("Supply an existing UTF-8 TSV file: ", file, call. = FALSE)
  file <- normalizePath(file, mustWork = TRUE)
  hash <- .lt_marine_sha(file)
  tab <- utils::read.delim(file, colClasses = "character", check.names = FALSE,
    quote = "", comment.char = "", na.strings = character(), fileEncoding = "UTF-8")
  if (anyDuplicated(names(tab)) || !all(required %in% names(tab)) || !nrow(tab))
    stop("TSV requires nonempty rows and unique columns: ", paste(required, collapse = ", "), call. = FALSE)
  if (!identical(hash, .lt_marine_sha(file))) stop("Input file changed while reading.", call. = FALSE)
  attr(tab, "input_source") <- list(path = file, sha256 = hash)
  tab
}

.lt_binary_number <- function(x, label) {
  # No factor coding, taxonomy inference, missing-value repair or recoding.
  value <- suppressWarnings(as.numeric(x))
  if (anyNA(value) || any(!is.finite(value)) || any(!nzchar(x)))
    stop(label, " must contain explicit finite numeric values, not blanks/NA/text.", call. = FALSE)
  value
}

lt_import <- function(file) {
  if (!is.character(file) || length(file) != 1L || !grepl("\\.rds$", file, ignore.case = TRUE))
    stop("lt_import requires a supported SplitAlignerR exchange .rds. Bare matrices and Perl files require a producer-supported versioned exchange; do not rename them.", call. = FALSE)
  lt_read_splitaligner_exchange(file)
}

lt_read_trait <- function(file, trait_name, taxon_column = "species") {
  .lt_assert_scalar_character(trait_name, "trait_name")
  .lt_assert_scalar_character(taxon_column, "taxon_column")
  if (identical(taxon_column, trait_name)) stop("Trait and taxon columns must differ.", call. = FALSE)
  tab <- .lt_binary_tsv(file, c(taxon_column, trait_name))
  .lt_assert_unique_ids(tab[[taxon_column]], "trait species")
  y <- .lt_binary_number(tab[[trait_name]], "Binary trait")
  if (!all(y %in% c(0, 0.5, 1)) || !all(c(0, 1) %in% y))
    stop("Binary trait requires numeric 0/1, optional explicit 0.5 exclusion, and both endpoint classes; recode deliberately outside LaTerra.", call. = FALSE)
  lt_trait(trait_name, tab[[taxon_column]], y, type = "binary",
    coding = list(negative = 0, positive = 1, excluded = 0.5),
    metadata = list(input_source = attr(tab, "input_source"), selected_column = trait_name))
}

lt_read_groups <- function(file) {
  tab <- .lt_binary_tsv(file, c("species", "genus"))
  .lt_assert_unique_ids(tab$species, "group species")
  if (anyNA(tab$genus) || any(!nzchar(tab$genus)))
    stop("Every species requires an explicit nonempty genus/group key; taxonomy is not inferred.", call. = FALSE)
  tab
}

lt_read_folds <- function(file) {
  tab <- .lt_binary_tsv(file, c("fold_id", "genus", "seed"))
  .lt_assert_unique_ids(tab$fold_id, "fold identifiers")
  .lt_assert_unique_ids(tab$genus, "fold genus keys")
  tab$seed <- .lt_binary_number(tab$seed, "Fold seed")
  if (any(tab$seed != trunc(tab$seed)) || any(abs(tab$seed) > .Machine$integer.max))
    stop("Fold seeds must be integer-valued R seeds.", call. = FALSE)
  tab
}

lt_binary_recipes <- function() {
  data.frame(recipe = "submitted_S2", version = "S2",
    validation_recipe = c("nested_phase12B", "nested_phase13", "global_screen_foldwise"),
    selection = c("fold-held branches excluded; nested_multiply arithmetic",
      "fold-held branches excluded; phase13_power arithmetic",
      "global supervised screen reused; not nested feature selection"),
    footprint = "inclusive gene q0.975/type7; arithmetic GE/BE; divide BE then GE; historical castor ASR; reference-minus-focal Welch; strict-prefix FDR0.01",
    stringsAsFactors = FALSE)
}

lt_run_binary <- function(input, trait, terminal_groups, folds, recipe,
                          validation_recipe, output_dir, cores = 1L, trait_name = NULL) {
  if (missing(recipe) || missing(validation_recipe))
    stop("Choose recipe and validation_recipe explicitly; inspect lt_binary_recipes().", call. = FALSE)
  .lt_s2_installed_runtime()
  if (is.character(input)) input <- lt_import(input)
  if (!is.list(input) || !all(c("matrix", "tree", "coordinates", "exchange", "provenance") %in% names(input)))
    stop("input must be the complete lt_import()/lt_read_splitaligner_exchange() result or an exchange RDS path; use lt_run_s2 for explicitly constructed inputs.", call. = FALSE)
  if (is.character(trait)) {
    if (is.null(trait_name)) stop("Choose trait_name explicitly when supplying a trait TSV path.", call. = FALSE)
    trait <- lt_read_trait(trait, trait_name)
  }
  if (is.character(terminal_groups)) terminal_groups <- lt_read_groups(terminal_groups)
  if (is.character(folds)) folds <- lt_read_folds(folds)
  matrix <- input$matrix
  # Preserve all six source arrays and producer evidence in supplied_inputs.rds;
  # a local metadata addition only, never a new numeric/coordinate authority.
  matrix$metadata$splitaligner_exchange <- input$exchange
  lt_run_s2(matrix, input$tree, input$coordinates, trait, terminal_groups,
    folds, recipe, validation_recipe, output_dir, cores)
}
