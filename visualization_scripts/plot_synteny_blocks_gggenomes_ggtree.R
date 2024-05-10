#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(argparse)
  library(gtools)
  library(scales)
  library(ggtree)
  library(treeio)
  library(phytools)
  library(dplyr)
  library(ggpubr)})
suppressWarnings(suppressPackageStartupMessages(library(gggenomes)))

# Example script for generating ntSynt synteny ribbon plots using gggenomes

# Parse the input arguments
parser <- ArgumentParser(description = "Plot the ntSynt synteny blocks and cladogram using gggenomes")
parser$add_argument("-s", "--sequences", help = "Input sequence lengths TSV", required = TRUE)
parser$add_argument("-l", "--links", help = "Synteny block links", required = TRUE)
parser$add_argument("-c", "--painting", help="File with chromosome painting information", required = FALSE)
parser$add_argument("--scale", help = "Length of scale bar in bases (default 1 Gbp)", default = 1e9,
                    required = FALSE, type = "double")
parser$add_argument("--tree", help = "Newick-formatted cladogram", required = TRUE)
parser$add_argument("--ratio", help = "Ratio adjustment for labels on left side of the ribbon plot. Increase if the labels are cut-off, decrease to decrease space between ribbon plot and cladogram",
                    default = 0.1, required = FALSE, type = "double")
parser$add_argument("-p", "--prefix",
                    help = "Output prefix for PNG image (default synteny_gggenomes_plot)", required = FALSE,
                    default = "synteny_gggenomes_plot")

args <- parser$parse_args()

# Read in and prepare sequences
sequences <- read.csv(args$sequences, sep = "\t", header = TRUE)

# https://stackoverflow.com/questions/32378108/using-gtoolsmixedsort-or-alternatives-with-dplyrarrange
input_order <- unique(sequences$bin_id)
input_chrom_order <- unique(sequences$seq_id)

mixedrank <- function(x) order(gtools::mixedorder(x))
sequences <- sequences %>%
  arrange(factor(bin_id, levels = input_order))

# Read in and prepare synteny links
links_ntsynt <- read.csv(args$links,
                         sep = "\t", header = TRUE)
links_ntsynt$seq_id <- factor(links_ntsynt$seq_id,
                              levels = input_chrom_order)
links_ntsynt <- links_ntsynt %>% arrange(factor(seq_id, levels = input_chrom_order))
links_ntsynt$seq_id2 <- as.character(links_ntsynt$seq_id2)
links_ntsynt$colour_block <- factor(links_ntsynt$colour_block,
                                    levels = input_chrom_order)

# Prepare scale bar data frame
scale <- args$scale

scale_bar <- tibble(x = c(0), xend = c(scale),
                    y = c(0), yend = c(0))

# Infer best units for scale bar
label <- paste(scale, "bp", sep = " ")
if (scale %% 1e9 == 0) {
  label <- paste(scale / 1e9, "Gbp", sep = " ")
} else if (scale %% 1e6 == 0) {
  label <- paste(scale / 1e6, "Mbp", sep = " ")
} else if (scale %% 1e3 == 0) {
  label <- paste(scale / 1e3, "kbp", sep = " ")
}

# Read in the data frame for chromosome painting features
painting <- read.csv(args$painting, sep = "\t", header = TRUE)

# Make the ribbon plot - these layers can be fully customized as needed!
make_plot <- function(links, sequences, painting, add_scale_bar = FALSE) {
  num_colours <- length(unique(links$colour_block))
  p <-  gggenomes(seqs = sequences, links = links, feats=painting)
  plot <- p + theme_gggenomes_clean(base_size = 15) +
    geom_link(aes(fill = colour_block), offset = 0, alpha = 0.5, size = 0.05) +
    geom_seq(size = 2, colour = "darkgrey") + # draw contig/chromosome lines
    geom_feat(aes(colour = colour_block), position = "identity", linewidth = 2) +
    geom_bin_label(aes(label = bin_id), size = 6, hjust = 0.9) + # label each bin
    #geom_seq_label(aes(label = seq_id), vjust = 1.1, size = 4) + # Can add seq labels if desired
    theme(axis.text.x = element_text(size = 25),
          legend.position = "bottom") +
    scale_fill_manual(values = hue_pal()(num_colours),
                      breaks = levels(links_ntsynt$colour_block)) +
    scale_colour_manual(values = hue_pal()(num_colours),
                      breaks = levels(links_ntsynt$colour_block)) +
    guides(fill = guide_legend(title = "", ncol = 10),
           colour = guide_legend(title = ""))
    xmax <- ggplot_build(plot)$layout$panel_params[[1]]$x.range[[2]]
    plot <- plot + xlim(0 - xmax*args$ratio, NA)

  if (add_scale_bar) {
    plot <- plot + geom_segment(data = scale_bar, aes(x = x, xend = xend, y = y, yend = yend),
                                linewidth = 1.5) +
      geom_text(data = scale_bar, aes(x = x + (xend / 2), y = y - 0.3, label = label)) +
      theme(axis.line.x = element_blank(), 
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank())
  }

  return(plot)

}

synteny_plot <- make_plot(links_ntsynt, sequences, painting, add_scale_bar = TRUE)

# Prepare the tree
ntsynt_tree <- treeio::read.newick(args$tree)
ntsynt_tree <- midpoint_root(ntsynt_tree)
ntsynt_ggtree <- ggtree(ntsynt_tree, branch.length = "none")

# Align the plots properly
synteny_y_range <- ggplot_build(synteny_plot)$layout$panel_params[[1]]$y.range

plots <- ggarrange(ntsynt_ggtree + scale_y_continuous(limits=synteny_y_range, expand=c(0, 0)), 
                   (synteny_plot %>% pick_by_tree(ntsynt_ggtree)),
                   common.legend = T, align = "hv",
                   widths = c(1, 10), legend = "bottom")


# Save the ribbon plot
ggsave(paste(args$prefix, ".png", sep = ""), plots,
       units = "cm", width = 50, height = 20, bg = "white")

cat(paste("Plot saved:", paste(args$prefix, ".png", sep = ""), "\n", sep = " "))
