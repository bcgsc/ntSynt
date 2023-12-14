# Example visualization scripts for ntSynt

Here, we provide basic examples scripts for generating ribbon plots and chromosome sequence painting plots to visualize synteny blocks computed by ntSynt. Each R script script can be seen as a starting point, and can be customized as needed. 

## Ribbon plots
The R package gggenomes (https://thackl.github.io/gggenomes/) is used to generate ribbon plots to visualize multi-genome synteny blocks.

### Steps:
1. Format the ntSynt synteny blocks for the ribbon plots R script using the provided scripts. Before running `format_blocks_gggenomes.sh`, add the `visualization_scripts` directory to your PATH.
```
(btl) [lcoombe@hpce706 visualization_scripts]$ format_blocks_gggenomes.sh 
Usage: format_blocks_gggenomes.sh <synteny blocks TSV> <prefix> <length threshold> <assembly to use for colour> <FAI> <FAI> [FAI..]
```
This script will generate two TSV files: {prefix}.links.tsv  {prefix}.sequence_lengths.tsv

2. Run the R script
* Required R packages: argparse, gggenomes, gtools, scales
```
usage: plot_synteny_blocks_gggenomes.R [-h] -s SEQUENCES -l LINKS [--scale SCALE] [-p PREFIX]

Plot the ntSynt synteny blocks using gggenomes

optional arguments:
  -h, --help            show this help message and exit
  -s SEQUENCES, --sequences SEQUENCES
                        Input sequence lengths TSV
  -l LINKS, --links LINKS
                        Synteny block links
  --scale SCALE         Length of scale bar in bases (default 1 Gbp)
  -p PREFIX, --prefix PREFIX
                        Output prefix for PNG image
```
Example:
![Example_gggenomes](https://github.com/bcgsc/ntSynt/blob/main/visualization_scripts/example_gggenomes.png)

* These plots are highly customizable, so edit/adapt the script as needed!