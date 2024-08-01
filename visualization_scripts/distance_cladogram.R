#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(treeio)
  library(ggpubr)
  library(ggtree)
  library(phytools)
  library(dplyr)})

# Generate a cladogram from ntSynt distance newick file

# Parse the input arguments
parser <- ArgumentParser(
  description = "Generate a cladodogram from ntSynt distance newick file,\noutput suggested file order for ribbon plots"
)
parser$add_argument("--nwk", help = "Newick file", required = TRUE)
parser$add_argument("-p", "--prefix",
                    help = "Output prefix for file ordering (default ntsynt_dists)", required = FALSE,
                    default = "ntsynt_dists")
parser$add_argument("--lim", help = "Ratio adjustment for xlimits (default 0.1)",
                    required = FALSE, default = 0.1, type = "double")
parser$add_argument("--png", help = "Output cladogram in PNG format", required = FALSE,
                    action = "store_true")

args <- parser$parse_args()
treefile <- treeio::read.newick(args$nwk)
treefile <- midpoint_root(treefile)

if (args$png) {
  tree_ggtree <- ggtree(treefile, branch.length = "none") +
    geom_tiplab() + hexpand(ratio = args$lim)

  ggsave(paste(args$prefix, ".cladogram.png", sep = ""), tree_ggtree, width = 25, height = 25, units = "cm",
         dpi = 300)
}

ordered_labels <- fortify(treefile) %>% filter(isTip) %>% arrange(desc(y)) %>% pull(label)
write.table(ordered_labels, paste(args$prefix, ".order.tsv", sep = ""),
            quote = FALSE, row.names = FALSE, col.names = FALSE)
