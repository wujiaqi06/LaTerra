sx_plot_bundle <- function() {
  skip_if_not_installed("SplitAlignerR")
  lt_read_splitaligner_exchange(test_path("fixtures", "splitaligner_exchange_v020.rds"))
}

# Independent tiny-tree descendant calculation, only for inspecting plots.
sx_test_edge_taxa <- function(tree, node) {
  if (node <= length(tree$tip.label)) return(tree$tip.label[node])
  sort(unlist(lapply(tree$edge[tree$edge[, 1] == node, 2], function(n)
    sx_test_edge_taxa(tree, n))), method = "radix")
}
sx_test_14_values <- function(a, gene) {
  tr <- attr(a, "display_tree")
  keys <- vapply(a$child_node, function(n) paste(sx_test_edge_taxa(tr, n), collapse = ","), character(1))
  expected <- if (gene == "gene01") c(A = 0, B = 2, C = 4, D = 5, E = 7, "A,B" = 3, "C,D" = 6) else
    c(A = 10, B = 20, C = 40, D = 50, E = 70, "A,B" = 30, "C,D" = 60)
  expect_setequal(keys, names(expected))
  expect_identical(a$gene_branch_value, unname(expected[keys]))
  expect_identical(a$tree_branch_length, rep(1, 7))
  expect_identical(a$node, as.integer(tr$edge[, 2]))
  expect_false(any(setdiff(tr$edge[, 1], tr$edge[, 2]) %in% a$node))
}

test_that("public adapter remaps values for actual changed node IDs and postorder", {
  x <- sx_plot_bundle(); original <- x
  trees <- list(x$tree, ape::reorder.phylo(x$tree, "postorder"),
    ape::read.tree(text = "((D:1,C:1):1,E:1,(B:1,A:1):1);"))
  for (tree in trees) for (gene in c("gene01", "gene02")) {
    a <- lt_annotate_splitaligner(x, gene, tree = tree)
    sx_test_14_values(a, gene)
    expect_identical(attr(a, "display_tree"), tree)
    expect_identical(a$display_label, a$display_alias)
    expect_identical(a$plot_value, a$gene_branch_value)
  }
  expect_identical(x, original)
  expect_false(identical(trees[[3]]$tip.label, x$tree$tip.label))
  s <- lt_annotate_splitaligner(x, "gene01", field = "coord_state")
  expect_true(all(s$plot_value == "observed"))
})

test_that("APE public call labels actual plotted edge rows, not traversal guesses", {
  x <- sx_plot_bundle()
  file <- tempfile(fileext = ".pdf"); grDevices::pdf(file, width = 10, height = 6)
  on.exit({ grDevices::dev.off(); unlink(file) }, add = TRUE)
  trees <- list(x$tree, ape::reorder.phylo(x$tree, "postorder"),
    ape::read.tree(text = "((D:1,C:1):1,E:1,(B:1,A:1):1);"))
  for (tree in trees) for (gene in c("gene01", "gene02")) {
    a <- lt_plot_splitaligner(x, gene, backend = "ape", tree = tree)
    sx_test_14_values(a, gene)
    actual <- get("last_plot.phylo", envir = getFromNamespace(".PlotPhyloEnv", "ape"))
    expect_identical(actual$edge[a$edge_row, 2], a$child_node)
    expect_identical(actual$edge[a$edge_row, 1], a$parent_node)
    expect_true(all(is.finite(actual$xx[a$node])))
    expect_true(all(is.finite(actual$yy[a$node])))
    expect_setequal(attr(a, "display_tree")$tip.label, LETTERS[1:5])
  }
})

test_that("ggtree public object retains all keyed values and draws exact label positions", {
  skip_if_not_installed("ggtree"); skip_if_not_installed("ggplot2")
  x <- sx_plot_bundle()
  tr <- ape::read.tree(text = "((D:1,C:1):1,E:1,(B:1,A:1):1);")
  for (gene in c("gene01", "gene02")) {
    warnings <- character()
    p <- withCallingHandlers(lt_plot_splitaligner(x, gene, backend = "ggtree", tree = tr),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    # Retained upstream ggtree3.14/ggplot4 compatibility warnings, not a blanket
    # production suppression. Any different warning fails this test.
    expect_true(all(grepl("deprecated|Arguments in .* must be used", warnings)))
    expect_true(inherits(p, "ggplot"))
    a <- attr(p, "lt_annotation"); sx_test_14_values(a, gene)
    expect_identical(nrow(p$data), 8L)
    label_layer <- Filter(function(layer) is.data.frame(layer$data) && "label" %in% names(layer$data) &&
      identical(layer$data$label, a$plot_label), p$layers)
    expect_length(label_layer, 1L)
    nodes <- as.data.frame(p$data); layer <- label_layer[[1]]$data
    expect_identical(layer$y_pos, nodes$y[match(a$child_node, nodes$node)])
    expect_identical(layer$x_pos, (nodes$x[match(a$child_node, nodes$node)] +
      nodes$x[match(a$parent_node, nodes$node)]) / 2)
  }
})

test_that("display edits work but stale science and ambiguous sources fail clearly", {
  x <- sx_plot_bundle(); x$matrix$branch_display_labels <- paste("edited", 1:7)
  x$matrix$metadata$note <- "display only"
  a <- lt_annotate_splitaligner(x, "gene01")
  expect_true(all(startsWith(a$display_label, "edited ")))
  sx_test_14_values(a, "gene01")
  y <- x; y$matrix$values[1, 1] <- 42
  expect_error(lt_annotate_splitaligner(y, "gene01"), "stale", class = "lt_error_exchange_stale")
  expect_silent(validate_lt_matrix(y$matrix))
  expect_error(lt_annotate_splitaligner(matrix(1, 1, 7), "gene01"), "Bare B", class = "lt_error_exchange_plot_input")
  expect_error(lt_annotate_splitaligner(x, "missing"), class = "lt_error_exchange_gene")
  expect_error(lt_annotate_splitaligner(x, "gene01", field = "made_up"), class = "lt_error_exchange_field")
  wrong <- ape::read.tree(text = "((D:1,A:1):1,E:1,(B:1,C:1):1);")
  expect_error(lt_annotate_splitaligner(x, "gene01", tree = wrong), class = "lt_error_exchange_annotation")
})

test_that("plot labels keep exact zero distinct from unavailable states and reasons", {
  x <- sx_local_test()
  x <- sx_cell_test(x, 1, 2, "mapped", NA_real_, "missing_branch_length", "TEST_SOURCE")
  x <- sx_cell_test(x, 1, 3, "NA_struct", NA_real_, "not_applicable", "TEST_ABSENCE")
  a <- lt_annotate_splitaligner(sx_import_test(x), "gene01")
  expect_true(any(grepl("NA [missing_branch_length]", a$plot_label, fixed = TRUE)))
  expect_true(any(grepl("NA [NA_struct]", a$plot_label, fixed = TRUE)))
  zero <- which(!is.na(a$gene_branch_value) & a$gene_branch_value == 0)
  expect_length(zero, 1L)
  expect_true(grepl("| 0", a$plot_label[zero], fixed = TRUE))
  expect_identical(a$coord_state[zero], "observed")
})
