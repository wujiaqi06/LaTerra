# Plot an imported gene without matching labels yourself

This development interface supports SplitAlignerR 0.1.0.9002 exchange schema
0.2.0-development. Install the matching SplitAlignerR backend and LaTerra into
your chosen R library. APE is required; ggtree and ggplot2 are optional.

```r
library(LaTerra)
x <- lt_read_splitaligner_exchange("exchange.rds")
lt_plot_splitaligner(x, gene = "gene01", backend = "ape")

# Alternatively, return a plot for further ggplot styling:
p <- lt_plot_splitaligner(x, gene = "gene01", backend = "ggtree")
print(p)

# States/reasons are selectable, too:
lt_plot_splitaligner(x, gene = "gene01", field = "coord_state")
```

Replace `gene01` with a scientific gene key in `x$matrix$gene_ids`. No manual
crosswalk, `match`, join, node renumbering or tip-label replacement is needed.
The labels combine the current display label with the selected field; the
gene value does not replace the reference tree's branch length. Zero is not
missing. The root has no incoming matrix cell. Save the complete `x` bundle
to retain original arrays, tree, result, ledgers and migration evidence.

For a compatible alternate tree representation, supply `tree = display_tree`.
The package calls the upstream public adapter afresh for that tree instance.
For a data frame rather than a plot use `lt_annotate_splitaligner(x, "gene01")`.
Additional fields: `value_reason`, `source_state`, `source_reason`,
`numeric_status`, `tree_branch_length`. No statistical results are calculated.

The named `marine302-user-bound-20260908-v1` profile is handled automatically
when it is already embedded in a source-bound migration exchange. You need
not edit or splice historical crosswalk tables. The native or migrated label
scheme is retained; an already converted exchange is not converted again.
Arbitrary old matrices with only B labels are ambiguous. Obtain a supported
SplitAlignerR export containing the original tree, scientific coordinates and
label provenance; LaTerra does not guess an old label scheme. Schema 0.1
files need explicit producer re-export to 0.2, not a locally changed hash.

Branch display labels and non-authoritative metadata can be edited. Changing
the numeric matrix requires coherent upstream re-export/re-import before
using this source-bound display route, so a plot cannot silently show stale
values. This restriction is not a golden-SHA barrier on ordinary lt_matrix.

Use a taller/wider graphics device for hundreds of branches. No labels are
silently hidden, and final manuscript layout remains a separate design task.
