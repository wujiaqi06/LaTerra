# Public display adapters. Upstream owns the graph/coordinate mapping; these
# functions only select fields and render its exact tree-instance annotation.
.lt_sx_plot_fields <- function() c("gene_branch_value", "coord_state",
  "value_reason", "source_state", "source_reason", "numeric_status", "tree_branch_length")

#' Automatically annotate a SplitAligner import for one gene
#' @param x Complete return value of [lt_read_splitaligner_exchange()].
#' @param gene One exact scientific gene key.
#' @param field One field from the upstream edge annotation.
#' @param tree Optional display phylo. Representation changes are mapped by
#'   SplitAlignerR, never by B numbers or user-supplied positional joins.
#' @return Complete edge annotation data frame, plus display labels and selected
#'   plot field. The display tree and backend provenance are attributes.
#' @export
lt_annotate_splitaligner <- function(x, gene, field = "gene_branch_value", tree = NULL) {
  .lt_sx_need(is.list(x) && inherits(x$matrix, "lt_matrix") &&
    inherits(x$exchange, "splitaligner_exchange"),
    "Use the complete lt_read_splitaligner_exchange() result. Bare B-labelled matrices lack a source tree and label provenance; re-export with SplitAlignerR 0.1.0.9002 export_splitaligner_exchange().",
    "plot_input")
  .lt_sx_need(.lt_sx_string(gene) && gene %in% x$matrix$gene_ids,
    "Select one exact gene key from x$matrix$gene_ids.", "gene")
  .lt_sx_need(.lt_sx_string(field) && field %in% .lt_sx_plot_fields(),
    paste("Unsupported field. Choose:", paste(.lt_sx_plot_fields(), collapse = ", ")), "field")
  ex <- x$exchange
  .lt_sx_header(ex); data <- .lt_sx_data(ex)
  .lt_sx_labels(ex, isTRUE(x$provenance$synthetic_migration_allowed))
  validate_lt_matrix(x$matrix)
  .lt_sx_need(identical(x$matrix$gene_ids, ex$axes$gene_ids) &&
    identical(x$matrix$branch_ids, ex$axes$coordinate_ids) &&
    all(vapply(c("values", "coord_state", "value_reason"), function(n)
      identical(x$matrix[[n]], data[[n]]), logical(1))),
    "Matrix science changed after import; this upstream annotation is stale. Re-export/re-import a coherent exchange. Display-label and non-authoritative metadata edits remain allowed.", "stale")
  backend <- .lt_sx_backend()
  if (is.null(tree)) tree <- x$tree
  # Public annotation includes the mandatory full upstream exchange validation.
  a <- tryCatch(SplitAlignerR::annotate_splitaligner_exchange(ex, gene, tree = tree),
    error = function(e) .lt_sx_need(FALSE, paste("Upstream public annotation rejected the source/display tree:",
      conditionMessage(e)), "annotation"))
  .lt_sx_need(is.data.frame(a) && nrow(a) == nrow(tree$edge) &&
    identical(a$edge_row, seq_len(nrow(tree$edge))) &&
    identical(a$child_node, as.integer(tree$edge[, 2])) &&
    identical(a$parent_node, as.integer(tree$edge[, 1])) &&
    !anyDuplicated(a$node) && !anyDuplicated(a$coordinate_id),
    "Upstream annotation is not a complete one-edge/one-coordinate mapping.", "annotation")
  j <- match(a$coordinate_id, x$matrix$branch_ids)
  .lt_sx_need(!anyNA(j) && identical(a$gene_branch_value, as.numeric(x$matrix$values[gene, j])),
    "Selected values disagree with the imported scientific keys.", "annotation")
  a$display_label <- x$matrix$branch_display_labels[j]
  a$plot_value <- a[[field]]
  value <- if (is.numeric(a$plot_value)) format(a$plot_value, digits = 6L, trim = TRUE) else a$plot_value
  value[is.na(a$plot_value)] <- "NA"
  if (identical(field, "gene_branch_value")) {
    missing <- is.na(a$gene_branch_value)
    reason <- ifelse(a$coord_state == "observed", a$numeric_status, a$coord_state)
    value[missing] <- paste0("NA [", reason[missing], "]")
  }
  a$plot_label <- paste(a$display_label, value, sep = " | ")
  attr(a, "display_tree") <- tree
  attr(a, "field") <- field
  attr(a, "payload") <- ex$contract$payload
  attr(a, "provenance") <- list(input_identity = ex$identity,
    display_tree_instance = unique(a$tree_instance_id), backend = backend,
    adapter = "SplitAlignerR::annotate_splitaligner_exchange", scientific_analysis_run = FALSE)
  a
}

.lt_sx_plot_colors <- function(a) {
  color <- rep("#2A5CAA", nrow(a))
  if (is.numeric(a$plot_value)) color[!is.na(a$plot_value) & a$plot_value == 0] <- "#A96500"
  color[is.na(a$plot_value)] <- "#666666"
  color
}

#' Plot a gene without manually matching branch labels or nodes
#' @inheritParams lt_annotate_splitaligner
#' @param backend APE draws to the current device; ggtree returns a ggplot.
#' @param label_size Positive label size (ggplot mm; APE equivalent cex).
#' @return APE invisibly returns the complete plotted annotation. ggtree returns
#'   a ggplot with the annotation in attribute `lt_annotation`.
#' @export
lt_plot_splitaligner <- function(x, gene, field = "gene_branch_value",
                                backend = c("ape", "ggtree"), tree = NULL, label_size = 3) {
  backend <- match.arg(backend)
  .lt_sx_need(is.numeric(label_size) && length(label_size) == 1L &&
    is.finite(label_size) && label_size > 0, "label_size must be finite and positive.", "plot")
  .lt_sx_need(is.list(x) && inherits(x$tree, "phylo"),
    "Supply the complete lt_read_splitaligner_exchange() result, not bare B labels.", "plot_input")
  if (is.null(tree)) tree <- x$tree
  .lt_sx_need(inherits(tree, "phylo"), "A display tree must be a phylo.", "plot")
  .lt_sx_need(requireNamespace("ape", quietly = TRUE), "Install ape to plot a tree.", "dependency")
  # plot.phylo uses a cladewise representation. Ask the upstream adapter for
  # THIS representation, rather than permuting annotation by B numbers.
  if (backend == "ape") tree <- ape::reorder.phylo(tree, "cladewise")
  a <- lt_annotate_splitaligner(x, gene, field, tree)
  lengths <- tree$edge.length
  .lt_sx_need(is.numeric(lengths) && length(lengths) == nrow(tree$edge) &&
    all(is.finite(lengths)) && all(lengths >= 0),
    "Plot requires finite nonnegative reference edge lengths; gene values are never substituted for missing geometry.", "plot_geometry")
  root <- setdiff(unique(tree$edge[, 1]), tree$edge[, 2])
  .lt_sx_need(length(root) == 1L && !root %in% a$node, "Root cannot carry an ordinary incoming matrix value.", "annotation")
  color <- .lt_sx_plot_colors(a)
  depth <- ape::node.depth.edgelength(tree); span <- max(depth)
  if (span == 0) span <- 1 # Display padding only; zero geometry is not modified.
  title <- paste(gene, "|", field)
  subtitle <- paste(length(tree$tip.label), "taxa;", nrow(a), "primitive branches; reference geometry retained")
  units <- attr(a, "payload")$units
  note <- paste("Edge label = display label |", field,
    if (field == "gene_branch_value") paste0(" (", units, ")") else "",
    "; 0 is not missing. Root has no matrix cell.")
  if (backend == "ape") {
    old <- graphics::par(mar = c(6, 2, 4.4, 2))
    on.exit(graphics::par(old), add = TRUE)
    ape::plot.phylo(tree, edge.color = "#777777", edge.width = 1.4,
      cex = label_size / 3.5, x.lim = c(-.15 * span, 1.5 * span),
      y.lim = c(.4, length(tree$tip.label) + .8), label.offset = .025 * span,
      no.margin = FALSE, show.node.label = FALSE, root.edge = FALSE)
    ape::edgelabels(a$plot_label, edge = a$edge_row, frame = "none", adj = c(.5, -.5),
      cex = label_size / 3.75, col = color)
    ape::nodelabels("root", node = root, frame = "none", adj = c(1.2, .5), cex = label_size / 3.75)
    graphics::title(main = title)
    graphics::mtext(subtitle, side = 3, line = .5, cex = .8)
    ape::axisPhylo(side = 1, backward = FALSE, cex.axis = .8)
    graphics::mtext("Reference-tree distance (stored tree units)", side = 1, line = 2.6, cex = .8)
    graphics::mtext(note, side = 1, line = 4.2, cex = .7)
    attr(a, "backend") <- "ape"
    return(invisible(a))
  }
  .lt_sx_need(requireNamespace("ggtree", quietly = TRUE) && requireNamespace("ggplot2", quietly = TRUE),
    "Install ggtree and ggplot2 for backend='ggtree'; APE is available separately.", "dependency")
  p <- ggtree::ggtree(tree, size = .5, color = "#777777")
  nodes <- as.data.frame(p$data)
  .lt_sx_need(!anyDuplicated(nodes$node), "ggtree emitted duplicate node keys.", "plot")
  child <- match(a$child_node, nodes$node); parent <- match(a$parent_node, nodes$node)
  .lt_sx_need(!anyNA(child) && !anyNA(parent), "ggtree lost a keyed display-tree node.", "plot")
  labels <- data.frame(x = (nodes$x[child] + nodes$x[parent]) / 2,
    y = nodes$y[child], label = a$plot_label, color = color)
  r <- nodes[nodes$node == root, , drop = FALSE]
  # aes uses explicit data columns; no scientific or positional B-label join.
  x_pos <- y_pos <- label <- colour <- NULL
  names(labels) <- c("x_pos", "y_pos", "label", "colour")
  root_label <- data.frame(x_pos = r$x, y_pos = r$y)
  p <- p + ggtree::geom_tiplab(size = label_size, offset = .025 * span, color = "#222222") +
    ggplot2::geom_text(data = labels, ggplot2::aes(x = x_pos, y = y_pos, label = label, color = colour),
      inherit.aes = FALSE, vjust = -.6, size = label_size) + ggplot2::scale_color_identity() +
    ggplot2::geom_text(data = root_label, ggplot2::aes(x = x_pos, y = y_pos), label = "root",
      inherit.aes = FALSE, hjust = 1.2, size = label_size, color = "#555555") +
    ggtree::theme_tree2() + ggplot2::xlim(-.15 * span, 1.5 * span) +
    ggplot2::labs(title = title, subtitle = subtitle, caption = note,
      x = "Reference-tree distance (stored tree units)") +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 18, 10, 18),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0))
  attr(p, "lt_annotation") <- a
  p
}
