test_that("S5 presentation uses significant digits with finite/sign/zero guards", {
  f <- LaTerra:::.lt_marine_table_six_significant
  expect_identical(f(c(0, -0, .123456789, -123.456789, 1.23456789e-5)),
    c(0, 0, .123457, -123.457, 1.23457e-5))
  x <- c(.Machine$double.xmin, .Machine$double.xmin * .Machine$double.eps, .Machine$double.xmax)
  expect_true(all(is.finite(f(x)) & f(x) > 0))
  for (x in list(NA_real_, Inf, -Inf, NaN, numeric(), "1", c(1, NA_real_)))
    expect_error(f(x), class = "lt_error_marine_reporting_payload")
})

test_that("S5 worksheet precision cannot flow back to markers, colors, lists or science", {
  # Standalone fresh fitted/annotation fixture, not historical values.
  headers <- unlist(LaTerra:::.lt_marine_reporting_layout("OUT_TABLE_S5", "predictor_annotation")$schema$headers)
  a <- as.data.frame(setNames(rep(list(rep("", 2)), length(headers)), headers),
    stringsAsFactors = FALSE, check.names = FALSE)
  a$gene <- c("g1", "g2"); a$display_module <- "fixture"; a$annotation_confidence <- "high"
  for (field in c("counted_in_Figure4C_circle_size", "displayed_in_approved_representative_gene_table",
    "recommended_for_Figure4C_original", "recommended_for_Figure4C_display", "reported_in_Supplementary_Table_S5"))
    a[[field]] <- c(TRUE, TRUE)
  a$recommended_for_TableS5_only <- a$keep_unassigned <- c(FALSE, FALSE)
  a$approved_representative_display_order_within_cell <- c(2, 1)
  fit <- function(b) list(beta = b, fit = list(beta = matrix(b, dimnames = list(names(b), "s0"))))
  full <- list(models = list(fix_marine_binary = fit(c(g1 = -.10000004, g2 = .10000003)),
    fix_aquatic_v2 = fit(c(g1 = .10000003, g2 = -.10000002))))
  full_before <- serialize(full, NULL, version = 3)
  fresh <- LaTerra:::.lt_marine_table_s5(full, a)
  before <- serialize(fresh, NULL, version = 3)
  view <- LaTerra:::.lt_marine_s5_worksheet_view(fresh$predictor_annotation)
  expect_identical(serialize(full, NULL, version = 3), full_before)
  expect_identical(serialize(fresh, NULL, version = 3), before)
  expect_identical(view$values$marine_coef, c(-.1, .1))
  expect_identical(view$values$max_abs_coef, c(.1, .1))
  expect_identical(view$values$coef_marker_marine_gt_0_1, c("#", "#"))
  expect_identical(view$values$label_color_direction, c("slow_negative_blue", "fast_positive_red"))
  expect_identical(view$values$Figure4C_cell_all_genes, rep("g1;g2", 2))
  expect_identical(view$values$Figure4C_cell_representative_genes, rep("g2;g1", 2))
  other <- setdiff(names(view$values), c("marine_coef", "aquatic_coef", "max_abs_coef"))
  expect_identical(view$values[other], fresh$predictor_annotation[other])
  expect_false(view$provenance$downstream_derivation_from_formatted_columns)
  expect_identical(view$provenance$authority$id, "TARGET001-REPORTING-CLARIFICATION002")
  bad <- fresh$predictor_annotation; bad$max_abs_coef <- .1
  expect_error(LaTerra:::.lt_marine_s5_worksheet_view(bad), class = "lt_error_marine_reporting_payload")
})
