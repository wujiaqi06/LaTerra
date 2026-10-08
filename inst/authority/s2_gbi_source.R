args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript build_endpointfix_gbi_matrix.R <matrix_type> <classified_matrix> <gbi_output> <gene_effect_output> <branch_effect_output> <baseline_summary_output>")
}

matrix_type <- args[1]
classified_matrix_file <- args[2]
gbi_output_file <- args[3]
gene_effect_output_file <- args[4]
branch_effect_output_file <- args[5]
baseline_summary_output_file <- args[6]

for (output_dir in unique(dirname(c(
  gbi_output_file,
  gene_effect_output_file,
  branch_effect_output_file,
  baseline_summary_output_file
)))) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
}

mean_no_na <- function(x) {
  x0 <- x[!is.na(x)]
  if (length(x0) == 0) {
    return(NA_real_)
  }
  mean(x0)
}

median_no_na <- function(x) {
  x0 <- x[!is.na(x)]
  if (length(x0) == 0) {
    return(NA_real_)
  }
  median(x0)
}

trim_upper_975 <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  cutoff <- as.numeric(stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE))
  x[x >= cutoff] <- NA
  x
}

apply_rows_preserve_matrix <- function(m, fun) {
  out <- t(vapply(seq_len(nrow(m)), function(i) fun(m[i, ]), FUN.VALUE = numeric(ncol(m))))
  rownames(out) <- rownames(m)
  colnames(out) <- colnames(m)
  out
}

classify_token <- function(x) {
  state_tokens <- c("NA", "NA_fuse", "NA_struct", "NA_topo", "residual_NA", "NaN", "Inf", "-Inf", "")
  out <- rep("finite_numeric", length(x))
  out[x %in% state_tokens] <- x[x %in% state_tokens]
  out[out == ""] <- "empty"
  suppressWarnings(num <- as.numeric(x))
  out[is.na(num) & !(x %in% state_tokens)] <- "unexpected_state"
  out[!is.na(num) & !is.finite(num)] <- x[!is.na(num) & !is.finite(num)]
  out
}

classified_raw <- utils::read.delim(
  classified_matrix_file,
  na.strings = character(0),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(classified_raw) < 2) {
  stop("Input classified matrix must contain one gene column and at least one branch column.")
}

gene_ids <- classified_raw[[1]]
branch_ids <- colnames(classified_raw)[-1]

if (any(is.na(gene_ids)) || any(gene_ids == "")) {
  stop("Input matrix has empty gene IDs.")
}
if (any(duplicated(gene_ids))) {
  stop("Input matrix has duplicated gene IDs.")
}
if (any(duplicated(branch_ids))) {
  stop("Input matrix has duplicated branch IDs.")
}

value_chr <- as.matrix(classified_raw[, -1, drop = FALSE])
state_matrix <- matrix(classify_token(as.vector(value_chr)), nrow = nrow(value_chr), ncol = ncol(value_chr))
unexpected_n <- sum(state_matrix == "unexpected_state")
if (unexpected_n > 0) {
  stop(sprintf("Input matrix contains %d unexpected nonnumeric states.", unexpected_n))
}

non_numeric_tokens <- c("NA", "NA_fuse", "NA_struct", "NA_topo", "residual_NA", "NaN", "Inf", "-Inf", "empty")
value_chr[value_chr %in% non_numeric_tokens] <- NA
branch_length <- suppressWarnings(matrix(as.numeric(value_chr), nrow = nrow(value_chr), ncol = ncol(value_chr)))
rownames(branch_length) <- gene_ids
colnames(branch_length) <- branch_ids

if (any(!is.na(branch_length) & !is.finite(branch_length))) {
  stop("Input matrix contains non-finite numeric values after conversion.")
}

branch_length_trimmed <- apply_rows_preserve_matrix(branch_length, trim_upper_975)

gene_effect <- data.frame(
  gene = rownames(branch_length_trimmed),
  GE_median = apply(branch_length_trimmed, 1, median_no_na),
  GE_mean = apply(branch_length_trimmed, 1, mean_no_na),
  n_present_before_trim = apply(branch_length, 1, function(x) sum(!is.na(x))),
  n_present_after_trim = apply(branch_length_trimmed, 1, function(x) sum(!is.na(x))),
  n_trimmed = apply(!is.na(branch_length) & is.na(branch_length_trimmed), 1, sum),
  zero_count_before_trim = apply(branch_length, 1, function(x) sum(x == 0, na.rm = TRUE)),
  zero_count_after_trim = apply(branch_length_trimmed, 1, function(x) sum(x == 0, na.rm = TRUE)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
gene_effect$ratio_median_over_mean <- gene_effect$GE_median / gene_effect$GE_mean

branch_effect <- data.frame(
  branch = colnames(branch_length_trimmed),
  BE_median = apply(branch_length_trimmed, 2, median_no_na),
  BE_mean = apply(branch_length_trimmed, 2, mean_no_na),
  n_present_before_trim = apply(branch_length, 2, function(x) sum(!is.na(x))),
  n_present_after_trim = apply(branch_length_trimmed, 2, function(x) sum(!is.na(x))),
  n_trimmed = apply(!is.na(branch_length) & is.na(branch_length_trimmed), 2, sum),
  zero_count_before_trim = apply(branch_length, 2, function(x) sum(x == 0, na.rm = TRUE)),
  zero_count_after_trim = apply(branch_length_trimmed, 2, function(x) sum(x == 0, na.rm = TRUE)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
branch_effect$ratio_median_over_mean <- branch_effect$BE_median / branch_effect$BE_mean

bl_div_be <- sweep(branch_length_trimmed, 2, branch_effect$BE_mean, "/")
gbi <- sweep(bl_div_be, 1, gene_effect$GE_mean, "/")
gbi[!is.finite(gbi)] <- NA_real_

gbi_df <- data.frame(gene = rownames(gbi), gbi, check.names = FALSE)

mapped_before <- sum(!is.na(branch_length))
mapped_after <- sum(!is.na(branch_length_trimmed))
trimmed_n <- mapped_before - mapped_after
baseline_summary <- data.frame(
  matrix_type = matrix_type,
  metric = c(
    "n_genes",
    "n_branches",
    "n_cells_total",
    "n_mapped_before_trim",
    "n_mapped_after_trim",
    "n_trimmed_cells",
    "global_trim_fraction_of_mapped_before",
    "n_zero_before_trim",
    "n_zero_after_trim",
    "n_genes_with_nonfinite_GE_mean",
    "n_branches_with_nonfinite_BE_mean",
    "n_gene_effect_mean_zero",
    "n_branch_effect_mean_zero",
    "gbi_na_fraction"
  ),
  value = c(
    nrow(branch_length),
    ncol(branch_length),
    length(branch_length),
    mapped_before,
    mapped_after,
    trimmed_n,
    if (mapped_before == 0) NA_real_ else trimmed_n / mapped_before,
    sum(branch_length == 0, na.rm = TRUE),
    sum(branch_length_trimmed == 0, na.rm = TRUE),
    sum(!is.finite(gene_effect$GE_mean) | is.na(gene_effect$GE_mean)),
    sum(!is.finite(branch_effect$BE_mean) | is.na(branch_effect$BE_mean)),
    sum(gene_effect$GE_mean == 0, na.rm = TRUE),
    sum(branch_effect$BE_mean == 0, na.rm = TRUE),
    mean(is.na(gbi))
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

utils::write.table(gbi_df, file = gbi_output_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
utils::write.table(gene_effect, file = gene_effect_output_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
utils::write.table(branch_effect, file = branch_effect_output_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")
utils::write.table(baseline_summary, file = baseline_summary_output_file, quote = FALSE, sep = "\t", row.names = FALSE, na = "NA")

run_metadata <- data.frame(
  field = c(
    "matrix_type",
    "input_path",
    "gbi_output_path",
    "gene_effect_output_path",
    "branch_effect_output_path",
    "baseline_summary_output_path",
    "timestamp",
    "formula",
    "trimming_rule",
    "non_numeric_state_handling",
    "n_genes",
    "n_branches",
    "total_missing_count",
    "total_zero_count"
  ),
  value = c(
    matrix_type,
    normalizePath(classified_matrix_file, mustWork = FALSE),
    normalizePath(gbi_output_file, mustWork = FALSE),
    normalizePath(gene_effect_output_file, mustWork = FALSE),
    normalizePath(branch_effect_output_file, mustWork = FALSE),
    normalizePath(baseline_summary_output_file, mustWork = FALSE),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    "GBI = branch_length_trimmed / BE_mean / GE_mean",
    "within_gene_upper_97.5_percent_to_NA",
    "NA, NA_fuse, NA_struct, NA_topo, residual_NA, NaN, Inf, -Inf, empty treated as missing",
    nrow(branch_length),
    ncol(branch_length),
    sum(is.na(branch_length)),
    sum(branch_length == 0, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
metadata_file <- sub("\\.tsv$", ".run_metadata.tsv", gbi_output_file)
utils::write.table(run_metadata, file = metadata_file, quote = FALSE, sep = "\t", row.names = FALSE)

message("GBI build complete: ", gbi_output_file)
