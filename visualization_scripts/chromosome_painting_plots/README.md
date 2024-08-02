## Chromosome sequence painting plots
ggplot2 is used to generate chromosome painting plots to visualize the multi-genome synteny blocks. Here, the synteny blocks for the other species are compared to a selected 'target' species. The segments based on the coordinate system of the 'target' species are coloured based on the chromosome of the other species. In addition, the relative orientation between the 'target' species and the other species is indicated by nudging the boxes up (forward) or down (reverse).

### Steps:
1. Format the ntSynt synteny blocks for the chromosome painting R script using the provided script.
```
usage: format_blocks_chromosome_painting.py [-h] [--convert CONVERT] --target TARGET synteny_tsv

Formatting synteny blocks for chromosome painting

positional arguments:
  synteny_tsv        ntSynt-formatted synteny blocks TSV

optional arguments:
  -h, --help         show this help message and exit
  --convert CONVERT  TSV file with desired conversions for assembly names (optional)
  --target TARGET    Target assembly name
```

This script will output the results in standard out - pipe the output to a file.

2. Run the R script
* Required R packages: argparse, ggplot2, dplyr, tidyr, gtools, scales
```
usage: plot_synteny_blocks-chromosome-painting.R [-h] -b BLOCKS [-g GAPS] [-p PREFIX]

Generate chromosome painting plots to visualize ntSynt synteny

optional arguments:
  -h, --help            show this help message and exit
  -b BLOCKS, --blocks BLOCKS
                        ntSynt-formatted synteny blocks
  -g GAPS, --gaps GAPS  TSV file with gap coordinates. Required headers:
                        chrom, chromStart, chromEnd (optional)
  -p PREFIX, --prefix PREFIX
                        Output prefix for PNG image (optional, default
                        synteny_chromosome)
```
* These plots can also be customized as needed (plot size, scale, comparing additional tools, etc.)

Example:
![Example_gggenomes](https://github.com/bcgsc/ntSynt/blob/distance_est/visualization_scripts/example_chromosome-painting.png)