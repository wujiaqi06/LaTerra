# Runnable installed example. All data are synthetic, not biological evidence.
library(LaTerra)
example_dir <- system.file("examples", "binary_small", package = "LaTerra", mustWork = TRUE)
run_dir <- tempfile("LaTerra-binary-run-")
report_dir <- tempfile("LaTerra-binary-report-")
imported <- lt_import(file.path(example_dir, "exchange.rds"))
result <- lt_run_binary(imported, trait = file.path(example_dir, "trait.tsv"),
  trait_name = "binary_trait", terminal_groups = file.path(example_dir, "groups.tsv"),
  folds = file.path(example_dir, "folds.tsv"), recipe = "submitted_S2",
  validation_recipe = "nested_phase12B", output_dir = run_dir)
saved <- lt_read_binary_run(result)
print(saved)
lt_report_binary(saved, report_dir)
cat("Scientific run:", run_dir, "\nExternal report:", report_dir, "\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  roc <- lt_plot_binary(saved, "roc")
  print(roc)
}
