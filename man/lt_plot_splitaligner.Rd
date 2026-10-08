\name{lt_plot_splitaligner}
\alias{lt_plot_splitaligner}
\alias{lt_annotate_splitaligner}
\title{Automatically annotate and plot a SplitAligner import}
\usage{
lt_annotate_splitaligner(x, gene, field = "gene_branch_value", tree = NULL)

lt_plot_splitaligner(x, gene, field = "gene_branch_value",
  backend = c("ape", "ggtree"), tree = NULL, label_size = 3)
}
\arguments{
\item{x}{Complete return value of \code{lt_read_splitaligner_exchange()}. Save this complete bundle, not just its numeric matrix.}
\item{gene}{One exact scientific key from \code{x$matrix$gene_ids}.}
\item{field}{One of \code{gene_branch_value}, \code{coord_state}, \code{value_reason}, \code{source_state}, \code{source_reason}, \code{numeric_status}, or \code{tree_branch_length}. These are imported fields, not new analysis results.}
\item{tree}{Optional display \code{phylo}; defaults to \code{x$tree}. Compatible rotation, tip/node renumbering, ladderization or edge reordering is mapped automatically by the supported upstream public adapter. No B-number or old-node-position assumption is made.}
\item{backend}{\code{ape} draws on the current graphics device. \code{ggtree} returns a ggplot object, printed normally or with \code{print()}.}
\item{label_size}{Finite positive label size in ggplot millimetres or equivalent APE text scaling. For dense trees, use a suitably large graphics device; labels are not silently dropped.}
}
\value{\code{lt_annotate_splitaligner()} returns a complete edge annotation data frame, retaining source/display aliases, scientific IDs, tree instance, edge/parent/child/node identities, geometry and all gene state/reason fields. Additional columns are \code{display_label}, \code{plot_value} and \code{plot_label}; attributes retain the display tree, selected field, payload and backend provenance. APE invisibly returns this annotation; ggtree returns a ggplot with the same annotation in \code{attr(p, "lt_annotation")}.}
\description{Plot directly from a validated imported object, without hand-written label/node joins. Uses \code{SplitAlignerR::annotate_splitaligner_exchange()} for the actual display tree and gene. No upstream graph algorithm is copied into La Terra and no scientific analysis runs.}
\details{Requires the supported SplitAlignerR 0.1.0.9002 schema-0.2 backend; ggtree and ggplot2 are optional and needed only for that backend. Original scientific values and reference-tree geometry remain separate. Zero is labelled zero; unavailable values carry the original reason/state, and the root never receives an ordinary matrix value. The plot requires finite nonnegative stored reference lengths; it never substitutes gene values for geometry. Species tip labels remain unchanged.

Editable branch display labels and non-authoritative metadata remain usable. Changing the matrix science without updating its upstream exchange makes this particular annotation dependency stale, not the entire matrix scientifically invalid. Bare ambiguous historical matrices, unsupported profiles, incompatible trees and missing source evidence fail with an explicit error. An already-migrated named-profile exchange is consumed as-is, without another conversion. This is a development interface and technical display, not a final publication design or release certification.}
\examples{
\dontrun{
x <- lt_read_splitaligner_exchange("exchange.rds")
lt_plot_splitaligner(x, gene = "gene01", backend = "ape")
p <- lt_plot_splitaligner(x, gene = "gene01", backend = "ggtree")
print(p)
lt_annotate_splitaligner(x, gene = "gene01", field = "coord_state")
}}
