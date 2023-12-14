#!/usr/bin/env Rscript
library(argparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)
library(scales)

# Example script for generating chromosome painting plots

# Parse the input arguments
parser <- ArgumentParser(description = "Generate chromosome painting plots to visualize ntSynt synteny")
parser$add_argument("-b", "--blocks", help = "ntSynt-formatted synteny blocks", required = TRUE)
parser$add_argument("-g", "--gaps",
                    help = "TSV file with gap coordinates. Required headers: chrom, chromStart, chromEnd (optional)",
                    required = FALSE)
parser$add_argument("-p", "--prefix", help = "Output prefix for PNG image (optional)", required = FALSE,
                    default = "synteny_chromosome-painting_plot")

args <- parser$parse_args()

# Read in the ntSynt blocks
ntsynt_blocks <- read.csv(args$blocks, sep = "\t", header = TRUE) %>%
  mutate(tool = "ntSynt")
block_rows_species <- tibble(other_species = unique(ntsynt_blocks$other_species))
ntsynt_blocks$target_chrom <- factor(ntsynt_blocks$target_chrom,
                                     levels = mixedsort(unique(ntsynt_blocks$target_chrom)))
ntsynt_blocks$other_chrom <- factor(ntsynt_blocks$other_chrom,
                                    levels = mixedsort(unique(ntsynt_blocks$other_chrom)))

# Read in gaps file if provided
if (!is.null(args$gaps)) {
  gaps <- read.csv(args$gaps, sep = "\t", header = TRUE) %>%
    mutate(target_chrom = chrom)
  gaps <- inner_join(gaps, ntsynt_blocks, by = "target_chrom") %>%
    select(chrom, chromStart, chromEnd, target_chrom)
  gaps <- crossing(gaps, block_rows_species)
  gaps$target_chrom <- factor(gaps$target_chrom,
                              levels = mixedsort(unique(gaps$target_chrom)))
}

# How much to nudge the blocks up or down based on their orientation
y_nudge <- 0.1

# Generate the chromosome painting plot
make_plot <- function(ntsynt_df, gaps = NULL) {
  my_plot <- ggplot() +
    geom_segment(data = ntsynt_df %>%
                   filter(relative_ori == "-"),
                 aes(x = target_start,
                     xend = target_end,
                     y = other_species,
                     yend = other_species,
                     colour = other_chrom),
                 position = position_nudge(y = -1 * y_nudge),
                 linewidth = 10) +
    geom_segment(data = ntsynt_df %>%
                   filter(relative_ori == "+"),
                 aes(x = target_start,
                     xend = target_end,
                     y = other_species,
                     yend = other_species,
                     colour = other_chrom),
                 position = position_nudge(y = y_nudge),
                 linewidth = 10) +
    theme_bw() +
    theme(text = element_text(size = 40),
          axis.text.y = element_text(margin = margin(l = 10),
                                     size = 30),
          axis.text.x = element_text(size = 28),
          strip.background = element_rect(fill = "gray95"),
          legend.position = "bottom") +
    ylab("") + xlab("Target chromosome position (Mbp)") +
    guides(colour = guide_legend(title = "Chromosome",
                                 override.aes = list(linewidth = 5),
                                 ncol = 10)) +
    facet_wrap(facets = vars(target_chrom)) +
    scale_x_continuous(labels = label_number(scale = 1e-6),
                       breaks = breaks_pretty(n = 3))

  if (!is.null(gaps)) {
    my_plot <- my_plot +
      geom_segment(data = gaps,
                   aes(x = chromStart, xend = chromEnd,
                       y = other_species,
                       yend = other_species),
                   linewidth = 5, colour = "black")
  }

  return(my_plot)
}

if (!is.null(args$gaps)) {
  painting_plot <- make_plot(ntsynt_blocks, gaps)
} else {
  painting_plot <- make_plot(ntsynt_blocks)
  cat("RAN without gaps")
}

ggsave(paste(args$prefix, ".chr-paint-plot.png", sep = ""), painting_plot,
       units = "cm", height = 40, width = 70, dpi = 300)
