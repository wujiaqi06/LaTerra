binary_example_files <- function() {
  dir <- system.file("examples", "binary_small", package = "LaTerra", mustWork = TRUE)
  stats::setNames(as.list(file.path(dir, c("exchange.rds", "trait.tsv", "groups.tsv", "folds.tsv"))),
    c("input", "trait", "terminal_groups", "folds"))
}

binary_need_stack <- function() {
  for (p in c("SplitAlignerR", "ape", "castor", "glmnet")) skip_if_not_installed(p)
  if (as.character(utils::packageVersion("SplitAlignerR")) != "0.1.0.9002") skip("Pinned 9002 producer required")
}

# Count points that reached the graphics tree, not rows merely kept by build().
binary_drawn_points <- function(grob) {
  if (inherits(grob, "points")) return(length(grob$x))
  if (inherits(grob, "gtable")) {
    panels <- grepl("^panel", grob$layout$name)
    children <- if (any(panels)) grob$grobs[panels] else grob$grobs
  } else children <- grob$children
  if (is.null(children) || !length(children)) return(0L)
  sum(vapply(children, binary_drawn_points, integer(1)))
}

binary_expect_all_oof_points <- function(plot, expected) {
  device_file <- tempfile("binary-point-render-", fileext = ".pdf")
  grDevices::pdf(device_file)
  device <- grDevices::dev.cur()
  on.exit({ grDevices::dev.off(device); unlink(device_file) }, add = TRUE)
  expect_warning(built <- ggplot2::ggplot_build(plot), NA)
  points <- built$data[[1]]
  expect_equal(nrow(points), expected)
  # Default alpha/fill NA means unspecified; required point aesthetics must exist.
  expect_false(anyNA(points[, c("x", "y", "shape", "colour", "size", "stroke")]))
  expect_warning(drawn <- ggplot2::ggplot_gtable(built), NA)
  expect_equal(binary_drawn_points(drawn), expected)
}

test_that("file readers preserve explicit keys, selection and complete ordered folds", {
  f <- binary_example_files()
  y <- lt_read_trait(f$trait, "binary_trait")
  expect_s3_class(y, "lt_trait")
  expect_equal(as.integer(table(y$values)), c(12L, 12L))
  expect_match(y$metadata$input_source$sha256, "^[0-9a-f]{64}$")
  groups <- lt_read_groups(f$terminal_groups)
  folds <- lt_read_folds(f$folds)
  expect_identical(folds$fold_id, c("first", "second", "third"))
  expect_equal(folds$seed, c(91, 92, 93))
  expect_setequal(folds$genus, groups$genus)
  expect_error(lt_read_trait(f$trait, "missing_column"), "requires")
  bad <- tempfile(fileext = ".tsv"); on.exit(unlink(bad))
  writeLines(c("species\tbinary_trait", "one\t0", "two\t2"), bad)
  expect_error(lt_read_trait(bad, "binary_trait"), "numeric 0/1")
  writeLines(c("species\tbinary_trait", "one\t0", "one\t1"), bad)
  expect_error(lt_read_trait(bad, "binary_trait"), "duplicate")
  writeLines(c("species\tbinary_trait", "one\t0", "two\tNA"), bad)
  expect_error(lt_read_trait(bad, "binary_trait"), "explicit finite")
  writeLines(c("fold_id\tgenus\tseed", "a\tone\t3.5"), bad)
  expect_error(lt_read_folds(bad), "integer-valued")
  expect_error(lt_import("bare_B123.tsv"), "Bare matrices and Perl")
  expect_error(lt_run_binary(), "Choose recipe")
  expect_identical(lt_binary_recipes()$validation_recipe,
    c("nested_phase12B", "nested_phase13", "global_screen_foldwise"))
  root <- normalizePath(test_path("..", ".."))
  expect_false(any(grepl("\\.rds$", LaTerra:::.lt_loggbi_doc_files(root))))
})

test_that("saved tables retain exact identity keys and literal quoted labels", {
  file <- tempfile(fileext = ".tsv"); on.exit(unlink(file))
  writeLines(c("species\tfold_id\tgenus\tresponse\tprobability\tstatus",
    "001\t01\t001\t0\t0.1\tfit", "NA\t02\t\"Group\"\t1\t0.9\tfit",
    "TRUE\t03\tTRUE\t1\tNA\tfit"), file)
  tab <- LaTerra:::.lt_binary_saved_table(file)
  expect_identical(tab$species, c("001", "NA", "TRUE"))
  expect_identical(tab$fold_id, c("01", "02", "03"))
  expect_identical(tab$genus, c("001", '"Group"', "TRUE"))
  expect_equal(tab$response, c(0, 1, 1))
  expect_equal(tab$probability, c(0.1, 0.9, NA_real_))
  keys <- c("gene", "gene_id", "branch", "branch_id", "species", "genus", "group", "fold_id", "run_id",
    "trait", "trait_id", "trait_name", "canonical_split_key", "terminal_taxon", "taxon_id", "scope", "role", "path", "sha256")
  ids <- c("001", "1", "NA", "TRUE", '"Key"')
  d <- stats::setNames(as.data.frame(rep(list(ids), length(keys)), stringsAsFactors = FALSE), keys)
  d$probability <- c(0.1, 0.9, NA_real_, 0.4, 0.5)
  LaTerra:::.lt_marine_write_tsv(d, file)
  restored <- LaTerra:::.lt_binary_saved_table(file)
  for (key in keys) expect_identical(restored[[key]], ids)
  expect_identical(restored$probability, d$probability)
  expect_equal(length(unique(restored$group)), 5L)
})

test_that("internal missing taxonomy is decoded without silently erasing unexpected identities", {
  file <- tempfile(fileext = ".tsv"); on.exit(unlink(file))
  writeLines(c("branch_id\tbranch_type\tterminal_taxon",
    "001\tterminal\tNA", "1\tinternal\tNA", "TRUE\tinternal\t"), file)
  tab <- LaTerra:::.lt_binary_saved_table(file)
  expect_identical(tab$branch_id, c("001", "1", "TRUE"))
  expect_identical(tab$terminal_taxon, c("NA", NA_character_, NA_character_))
  writeLines(c("branch_id\tbranch_type\tterminal_taxon",
    "1\tinternal\tUNEXPECTED_VALUE"), file)
  bytes <- readBin(file, "raw", n = file.info(file)$size)
  expect_error(LaTerra:::.lt_binary_saved_table(file), "unexpected terminal_taxon.*no identity is silently erased")
  expect_identical(readBin(file, "raw", n = file.info(file)$size), bytes)
})

test_that("installed runtime is namespace-bound and required before any S2 work or output creation", {
  expect_identical(LaTerra:::.lt_s2_installed_runtime(),
    file.path(getNamespaceInfo(asNamespace("LaTerra"), "path"), "R", "LaTerra.rdb"))
  local_mocked_bindings(.lt_s2_installed_runtime = function(...) stop("Install the source package; load_all unsupported"),
    .lt_s2_user_coordinates = function(...) stop("coordinate work forbidden"),
    .lt_marine_annotate = function(...) stop("ASR forbidden"),
    .lt_marine_baseline = function(...) stop("baseline forbidden"),
    .lt_marine_screen = function(...) stop("screen forbidden"),
    .lt_s2_fit_design = function(...) stop("fit forbidden"), .package = "LaTerra")
  output <- tempfile("preflight-run-")
  expect_error(lt_run_s2(recipe = "submitted_S2", validation_recipe = "nested_phase12B", output_dir = output),
    "Install the source package; load_all unsupported")
  expect_error(lt_run_binary(recipe = "submitted_S2", validation_recipe = "nested_phase12B", output_dir = output),
    "Install the source package; load_all unsupported")
  expect_false(file.exists(output))
})

test_that("public file workflow equals unchanged S2 and reports never refit or alter runs", {
  binary_need_stack()
  f <- binary_example_files()
  imported <- lt_import(f$input)
  expect_equal(dim(imported$matrix$values), c(32L, 27L))
  original <- serialize(imported, NULL, version = 3)
  parent <- tempfile("binary-test-"); dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  a <- file.path(parent, "public"); b <- file.path(parent, "lowlevel")
  result <- do.call(lt_run_binary, c(as.list(f), list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", trait_name = "binary_trait", output_dir = a)))
  expect_identical(serialize(imported, NULL, version = 3), original)
  low <- lt_run_s2(imported$matrix, imported$tree, imported$coordinates,
    lt_read_trait(f$trait, "binary_trait"), lt_read_groups(f$terminal_groups), lt_read_folds(f$folds),
    "submitted_S2", "nested_phase12B", b)
  for (file in c("S2_baseline_GBI.rds", "validation.rds")) expect_identical(readRDS(file.path(a, file)), readRDS(file.path(b, file)))
  for (file in c("branch_states.tsv", "screen_tested.tsv", "screen_significant.tsv", "OOF_predictions.tsv", "validation_summary.tsv"))
    expect_identical(readLines(file.path(a, file)), readLines(file.path(b, file)))
  expect_identical(readRDS(file.path(a, "supplied_inputs.rds"))$matrix$metadata$splitaligner_exchange, imported$exchange)
  files <- list.files(a, full.names = TRUE)
  before <- vapply(files, LaTerra:::.lt_marine_sha, character(1))
  # If saved views accidentally call any scientific operators this test fails.
  local_mocked_bindings(.lt_marine_baseline = function(...) stop("refit forbidden"),
    .lt_marine_screen = function(...) stop("rescreen forbidden"),
    .lt_s2_fit_design = function(...) stop("refit forbidden"), .package = "LaTerra")
  saved <- lt_read_binary_run(result)
  expect_s3_class(saved, "lt_binary_result")
  expect_identical(saved$integrity, "verified_file_hashes")
  expect_equal(saved$diagnostics$oof_coverage$n_unique_predicted_terminals, 24)
  expect_equal(saved$diagnostics$trim_genes$n_trimmed,
    unname(rowSums(readRDS(file.path(a, "S2_baseline_GBI.rds"))$trim_mask)))
  expect_true(is.data.frame(saved$diagnostics$input_cells))
  expect_true(is.data.frame(saved$diagnostics$input_exact_extrema))
  expect_identical(summary(saved)$source, normalizePath(a))
  report <- file.path(parent, "report")
  expect_identical(lt_report_binary(saved, report), normalizePath(report))
  expect_true(file.exists(file.path(report, "report.md")))
  expect_true(file.exists(file.path(report, "D4_validation_folds.tsv")))
  text <- readLines(file.path(report, "report.md"))
  expect_true(any(grepl("indexes the supplied_inputs.rds snapshot", text, fixed = TRUE)))
  expect_true(any(grepl("splitaligner_import", text, fixed = TRUE)))
  expect_true(any(grepl("attr(folds, 'input_source')", text, fixed = TRUE)))
  supplied <- readRDS(file.path(a, "supplied_inputs.rds"))
  expect_identical(supplied$trait$metadata$input_source, lt_read_trait(f$trait, "binary_trait")$metadata$input_source)
  expect_identical(attr(supplied$terminal_groups, "input_source"), attr(lt_read_groups(f$terminal_groups), "input_source"))
  expect_identical(attr(supplied$folds, "input_source"), attr(lt_read_folds(f$folds), "input_source"))
  expect_identical(readRDS(file.path(a, "run.rds"))$provenance$software_version$installed_runtime_sha256,
    LaTerra:::.lt_marine_sha(system.file("R/LaTerra.rdb", package = "LaTerra", mustWork = TRUE)))
  expect_error(lt_report_binary(saved, report), "exists")
  expect_error(lt_report_binary(saved, file.path(a, "figures")), "outside")
  if (requireNamespace("ggplot2", quietly = TRUE)) for (view in c("roc", "oof", "trim")) {
    plot <- lt_plot_binary(saved, view)
    expect_s3_class(plot, "ggplot")
    expect_false(grepl("PARTIAL / FAILED", paste(plot$labels$title, plot$labels$caption), fixed = TRUE))
    if (view == "oof") {
      d <- saved$tables$OOF_predictions
      binary_expect_all_oof_points(plot, sum(is.finite(d$probability) & d$response %in% c(0, 1)))
    } else expect_true(nrow(ggplot2::ggplot_build(plot)$data[[1]]) > 0)
  }
  expect_identical(vapply(list.files(a, full.names = TRUE), LaTerra:::.lt_marine_sha, character(1)), before)
  writeLines(c(readLines(file.path(a, "screen_tested.tsv")), ""), file.path(a, "screen_tested.tsv"))
  expect_error(lt_read_binary_run(a), "output hash mismatch")
  unlink(file.path(a, "gene_effect.tsv"))
  expect_error(lt_read_binary_run(a), "Completed-run integrity error: missing")
})

test_that("OOF rendering preserves all rows from an eight-group public run", {
  binary_need_stack()
  skip_if_not_installed("ggplot2")
  f <- binary_example_files()
  input <- lt_import(f$input)
  trait <- lt_read_trait(f$trait, "binary_trait")
  groups <- lt_read_groups(f$terminal_groups)
  keys <- c("001", "1", "003", "004", "005", "006", "007", "008")
  groups$genus <- rep(keys, 3)
  folds <- data.frame(fold_id = sprintf("fold%02d", 1:8),
    genus = keys, seed = 101:108)
  parent <- tempfile("binary-eight-groups-"); dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  run_dir <- file.path(parent, "run")
  run <- lt_run_binary(input, trait, groups, folds, "submitted_S2", "nested_phase12B", run_dir)
  files <- list.files(run_dir, full.names = TRUE)
  before <- vapply(files, LaTerra:::.lt_marine_sha, character(1))
  local_mocked_bindings(.lt_marine_baseline = function(...) stop("refit forbidden"),
    .lt_marine_screen = function(...) stop("rescreen forbidden"),
    .lt_s2_fit_design = function(...) stop("refit forbidden"), .package = "LaTerra")
  saved <- lt_read_binary_run(run)
  eligible <- saved$tables$OOF_predictions
  eligible <- eligible[is.finite(eligible$probability) & eligible$response %in% c(0, 1), , drop = FALSE]
  expect_equal(nrow(eligible), 24L)
  expect_equal(length(unique(eligible$genus)), 8L)
  expect_identical(unique(saved$tables$fold_ledger$group), keys)
  expect_equal(length(unique(saved$tables$fold_ledger$group)), 8L)
  report <- file.path(parent, "report")
  lt_report_binary(saved, report)
  expect_identical(readBin(file.path(run_dir, "fold_ledger.tsv"), "raw", n = file.info(file.path(run_dir, "fold_ledger.tsv"))$size),
    readBin(file.path(report, "fold_ledger.tsv"), "raw", n = file.info(file.path(report, "fold_ledger.tsv"))$size))
  exported <- LaTerra:::.lt_binary_saved_table(file.path(report, "fold_ledger.tsv"))
  expect_setequal(exported$group, eligible$genus)
  plot <- lt_plot_binary(saved, "oof")
  expect_identical(plot$data, eligible)
  binary_expect_all_oof_points(plot, nrow(eligible))
  expect_warning(ggplot2::ggsave(file.path(parent, "oof.png"), plot,
    width = 8, height = 6, dpi = 72), NA)
  expect_identical(vapply(files, LaTerra:::.lt_marine_sha, character(1)), before)
})

test_that("large-group OOF views have no aesthetic-induced omissions", {
  skip_if_not_installed("ggplot2")
  for (n_groups in c(24L, 128L, 512L)) {
    # These are rendering fixtures, not newly fitted scientific predictions.
    d <- data.frame(species = sprintf("tip%04d", seq_len(2L * n_groups)),
      genus = rep(sprintf("group%04d", seq_len(n_groups)), each = 2L),
      response = rep(c(0, 1), n_groups), probability = rep(c(0.2, 0.8), n_groups))
    before <- serialize(d, NULL, version = 3)
    saved <- list(run_dir = "render-only-fixture", completed = TRUE, tables = list(OOF_predictions = d))
    local_mocked_bindings(lt_read_binary_run = function(...) saved, .package = "LaTerra")
    plot <- lt_plot_binary("render-only-fixture", "oof")
    expect_identical(plot$data, d)
    expect_match(plot$labels$subtitle, paste(n_groups, "supplied groups"), fixed = TRUE)
    binary_expect_all_oof_points(plot, 2L * n_groups)
    expect_identical(serialize(d, NULL, version = 3), before)
  }
})

test_that("failed runs report partial dependencies and completed hash changes are errors", {
  binary_need_stack()
  parent <- tempfile("binary-failure-"); dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  f <- binary_example_files()
  local_mocked_bindings(.lt_marine_screen = function(...) stop("deliberate test failure"), .package = "LaTerra")
  # The directory name deliberately contains no failure words: labels must
  # reflect RUN_STATUS.yaml, not happen to match a user-selected pathname.
  run <- file.path(parent, "run")
  expect_error(do.call(lt_run_binary, c(as.list(f), list(recipe = "submitted_S2",
    validation_recipe = "nested_phase12B", trait_name = "binary_trait", output_dir = run))), "deliberate test failure")
  saved <- lt_read_binary_run(run)
  expect_false(saved$completed)
  expect_identical(saved$status$stage, "screen")
  expect_identical(saved$diagnostics$dependencies$status[3:4], c("incomplete_failed_stage", "not_run"))
  report <- file.path(parent, "partial_report")
  lt_report_binary(run, report)
  expect_true(any(grepl("PARTIAL failed run", readLines(file.path(report, "report.md")))))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    files <- list.files(run, full.names = TRUE)
    before <- vapply(files, LaTerra:::.lt_marine_sha, character(1))
    local_mocked_bindings(.lt_marine_baseline = function(...) stop("refit forbidden"),
      .lt_s2_fit_design = function(...) stop("refit forbidden"), .package = "LaTerra")
    plot <- lt_plot_binary(saved, "trim")
    expect_match(plot$labels$title, "PARTIAL / FAILED (stage: screen)", fixed = TRUE)
    expect_match(plot$labels$caption, "PARTIAL / FAILED (stage: screen)", fixed = TRUE)
    expect_warning(ggplot2::ggsave(file.path(parent, "partial_trim.png"), plot, width = 8, height = 6, dpi = 72), NA)
    expect_error(lt_plot_binary(saved, "oof"), "not_run_or_incomplete")
    expect_error(lt_plot_binary(saved, "roc"), "not_run_or_incomplete")
    expect_identical(vapply(files, LaTerra:::.lt_marine_sha, character(1)), before)
  }
})

test_that("provenance-stage failures retain OOF plots with explicit partial status and no refit", {
  binary_need_stack()
  skip_if_not_installed("ggplot2")
  f <- binary_example_files()
  input <- lt_import(f$input)
  trait <- lt_read_trait(f$trait, "binary_trait")
  groups <- lt_read_groups(f$terminal_groups); folds <- lt_read_folds(f$folds)
  parent <- tempfile("binary-late-stage-"); dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  run_dir <- file.path(parent, "run")
  local_mocked_bindings(lt_empty_provenance = function(...) stop("deliberate provenance failure"), .package = "LaTerra")
  expect_error(lt_run_binary(input, trait, groups, folds, "submitted_S2", "nested_phase12B", run_dir),
    "deliberate provenance failure")
  saved <- lt_read_binary_run(run_dir)
  expect_false(saved$completed)
  expect_identical(saved$status$stage, "provenance")
  expect_true(nrow(saved$tables$OOF_predictions) > 0L)
  files <- list.files(run_dir, full.names = TRUE)
  before <- vapply(files, LaTerra:::.lt_marine_sha, character(1))
  local_mocked_bindings(.lt_marine_baseline = function(...) stop("refit forbidden"),
    .lt_marine_screen = function(...) stop("rescreen forbidden"),
    .lt_s2_fit_design = function(...) stop("refit forbidden"), .package = "LaTerra")
  for (view in c("oof", "roc", "trim")) {
    plot <- lt_plot_binary(saved, view)
    expect_match(plot$labels$title, "PARTIAL / FAILED (stage: provenance)", fixed = TRUE)
    expect_match(plot$labels$caption, "PARTIAL / FAILED (stage: provenance)", fixed = TRUE)
    expect_warning(ggplot2::ggsave(file.path(parent, paste0(view, ".png")), plot, width = 8, height = 6, dpi = 72), NA)
  }
  expect_identical(vapply(files, LaTerra:::.lt_marine_sha, character(1)), before)
})
