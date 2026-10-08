test_that("S2 display BH does not replace selected-prefix membership", {
  d <- data.frame(gene = c("g1", "g2", "g3"), tvalue = c(2, -3, 4),
    pvalue = c(.01, .03, .04), mean1 = c(2, 3, 5), mean2 = c(1, 4, 4),
    marine = c(5, 6, 7), non_marine = c(8, 9, 10))
  s <- d[c(2, 1), ]
  before <- serialize(list(d, s), NULL, version = 3)
  x <- LaTerra:::.lt_marine_table_s2(d, s, "marine_binary")
  expect_identical(x$gene_id, c("g2", "g1"))
  expect_equal(x$FDR, c(.04, .03))
  expect_identical(x$direction, c("fast", "slow"))
  expect_identical(x$mean_trait_state_1, s$mean2)
  expect_identical(x$n_state_1, s$marine)
  expect_identical(serialize(list(d, s), NULL, version = 3), before)
  changed <- s; changed$tvalue[1] <- 30
  expect_error(LaTerra:::.lt_marine_table_s2(d, changed, "marine_binary"),
    class = "lt_error_marine_reporting_keys")
  expect_error(LaTerra:::.lt_marine_table_bh(c(.01, NA_real_)), class = "lt_error_marine_reporting_payload")
  expect_equal(LaTerra:::.lt_marine_table_bh(c(.01, .01, .2)), c(.015, .015, .2))
})

test_that("admitted percentage formatter derives values and fails on undefined inputs", {
  f <- LaTerra:::.lt_marine_table_percent
  expect_identical(f(876 / 894), "98.0%")
  expect_identical(f(.973), "97.3%")
  expect_identical(f(0), "0.0%")
  expect_identical(f(1), "100.0%")
  for (x in list(NA_real_, NaN, Inf, -Inf, -.1, 1.1, numeric(), c(.1, .2), "0.98"))
    expect_error(f(x), class = "lt_error_marine_reporting_payload")
})

test_that("coefficient wording preserves both axes and unrounded mixed signs", {
  f <- LaTerra:::.lt_marine_table_direction
  expect_identical(f(1, -2, TRUE, TRUE), "marine fast direction positive; aquatic slow direction negative")
  expect_identical(f(-1, 2, TRUE, TRUE), "marine slow direction negative; aquatic fast direction positive")
  expect_identical(f(1e-300, 0, TRUE, FALSE), "marine fast direction positive")
  expect_identical(f(0, -1e-300, FALSE, TRUE), "aquatic slow direction negative")
  for (v in list(0, -0, NA_real_, NaN, Inf, -Inf, numeric(), "1"))
    expect_error(f(v, 0, TRUE, FALSE), class = "lt_error_marine_reporting_payload")
  expect_error(f(1, 0, "True", FALSE), class = "lt_error_marine_reporting_payload")
  expect_error(f(0, 0, FALSE, FALSE), class = "lt_error_marine_reporting_payload")
  expect_identical(LaTerra:::.lt_marine_reporting_clarification()$qualifier,
    "presentation_rule_frozen_by_main_control_original_generator_unlocated")
})

s5_reporting_fixture <- function() {
  headers <- unlist(LaTerra:::.lt_marine_reporting_layout("OUT_TABLE_S5", "predictor_annotation")$schema$headers)
  a <- as.data.frame(setNames(rep(list(rep("", 3)), length(headers)), headers),
    stringsAsFactors = FALSE, check.names = FALSE)
  a$gene <- c("g1", "g2", "g3")
  a$display_module <- "test module"
  a$annotation_confidence <- "high"
  for (f in c("counted_in_Figure4C_circle_size", "displayed_in_approved_representative_gene_table",
    "recommended_for_Figure4C_original", "recommended_for_Figure4C_display", "reported_in_Supplementary_Table_S5"))
    a[[f]] <- rep("True", 3)
  a$recommended_for_TableS5_only <- a$keep_unassigned <- rep("False", 3)
  # Shared g1/g2 are deliberately in an approved representative order that
  # differs from new coefficient order. It must not be reselected/reordered.
  a$approved_representative_display_order_within_cell <- c(1, 2, 1)
  fit <- function(b) list(beta = b, fit = list(beta = matrix(b, dimnames = list(names(b), "s0"))))
  list(annotation = a, full = list(models = list(
    fix_marine_binary = fit(c(g1 = .3, g2 = -.6, g3 = 0)),
    fix_aquatic_v2 = fit(c(g1 = .1, g2 = .2, g3 = -.7)))))
}

test_that("S5 uses fitted support and only approved annotation fields", {
  f <- s5_reporting_fixture()
  x <- LaTerra:::.lt_marine_table_s5(f$full, f$annotation)
  expect_identical(x$predictor_annotation$gene, c("g1", "g2", "g3"))
  expect_equal(x$predictor_annotation$marine_coef, c(.3, -.6, 0))
  expect_identical(x$predictor_annotation$Figure4C_cell_representative_genes[1:2], rep("g1;g2", 2))
  expect_identical(x$predictor_annotation$Figure4C_cell_all_genes[1:2], rep("g2;g1", 2))
  expect_equal(x$display_module_summary$n_shared, 2)
  bad_source <- f$annotation
  for (c in c("marine_coef", "aquatic_coef", "max_abs_coef")) bad_source[[c]] <- 1e99
  bad_source$selected_in_marine <- "False"
  bad_source$coefficient_direction_summary <- "stale archived answer"
  same <- LaTerra:::.lt_marine_table_s5(f$full, bad_source)
  expect_identical(x, same)
  expect_false(x$provenance$representative_selection_recomputed)
  expect_false(x$provenance$original_numeric_annotation_columns_consumed)
  bad_fit <- f$full; bad_fit$models$fix_marine_binary$beta[1] <- 77
  expect_error(LaTerra:::.lt_marine_table_s5(bad_fit, f$annotation), class = "lt_error_marine_reporting_dependency")
  expect_error(LaTerra:::.lt_marine_table_s5(f$full, f$annotation[-1, ]), class = "lt_error_marine_reporting_keys")
})

test_that("S5 stops at a fresh priority tie and never falls back to gene order", {
  f <- s5_reporting_fixture()
  f$full$models$fix_marine_binary$beta[1] <- .6
  f$full$models$fix_marine_binary$fit$beta[1, 1] <- .6
  expect_error(LaTerra:::.lt_marine_table_s5(f$full, f$annotation), class = "lt_error_marine_reporting_priority_tie")
  f <- s5_reporting_fixture(); f$annotation$approved_representative_display_order_within_cell[1] <- 2
  expect_error(LaTerra:::.lt_marine_table_s5(f$full, f$annotation), class = "lt_error_marine_reporting_keys")
  f <- s5_reporting_fixture(); f$annotation$counted_in_Figure4C_circle_size[1] <- "yes"
  expect_error(LaTerra:::.lt_marine_table_s5(f$full, f$annotation), class = "lt_error_marine_reporting_payload")
})

test_that("S4 final display aliases preserve original metric keys and historical context gate", {
  d <- data.frame(gene = paste0("g", 1:4), tvalue = c(2, 3, -1, -2))
  # A tiny different context must reach the historical-context guard, rather
  # than failing on recovered snake-case source IDs versus final display IDs.
  expect_error(LaTerra:::.lt_marine_table_overlap(d, d, d[1:3, ], d[c(1, 2, 4), ]),
    class = "lt_error_marine_reporting_historical_context")
})

test_that("S6 retains the original final builder's True/False textual column", {
  expect_identical(LaTerra:::.lt_marine_table_boolean_text(c(TRUE, FALSE)), c("True", "False"))
  for (x in list(c(TRUE, NA), c(1, 0), c("True", "False")))
    expect_error(LaTerra:::.lt_marine_table_boolean_text(x), class = "lt_error_marine_reporting_payload")
})
