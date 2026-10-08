# Run with your installed package/library selected by R_LIBS, for example:
# Rscript plot_splitaligner.R exchange.rds gene01 ape gene01.png
# Optional fifth argument: a compatible display-tree Newick file.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(4L, 5L))
library(LaTerra)
x <- lt_read_splitaligner_exchange(args[1])
display_tree <- if (length(args) == 5L) ape::read.tree(args[5]) else NULL
grDevices::png(args[4], width = 1600, height = 900, res = 150)
p <- lt_plot_splitaligner(x, gene = args[2], backend = args[3], tree = display_tree)
if (args[3] == "ggtree") print(p)
grDevices::dev.off()
