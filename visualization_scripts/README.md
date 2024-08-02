# Visualizing ntSynt synteny blocks

Here, we provide an easy-to-use pipline for generating ribbon plots combined with chromosome painting to visualize the output synteny blocks from ntSynt.

This flexible pipeline implements numerous features, including:
* Option to normalize the strands of input chromosomes, based on a target assembly
* Evidence-guided ordering of assemblies from top-to-bottom, based on an input tree structure or distance estimates from the synteny blocks
* Sorting chromosomes right-to-left based on mappings to other assemblies
* Colouring both the ribbons and chromosomes based on the chromosomes in the target (top) assembly

These features ensure that the output ribbon plots (powered by [gggenomes](https://thackl.github.io/gggenomes/)) are as easily understandable and information-rich as possible.

### Dependencies
* quicktree
* R packages:
  * [gggenomes](https://github.com/thackl/gggenomes)
  * [treeio](https://www.bioconductor.org/packages/release/bioc/html/treeio.html)
  * [ggpubr](https://rpkgs.datanovia.com/ggpubr/)
  * [ggtree](https://github.com/YuLab-SMU/ggtree)
  * [phytools](https://cran.r-project.org/web/packages/phytools/index.html)
  * [dplyr](https://dplyr.tidyverse.org/)
  * [argparse](https://cran.r-project.org/web/packages/argparse/index.html)
  * [scales](https://scales.r-lib.org/)
  * [stringr](https://stringr.tidyverse.org/)

#### Installing dependencies using conda
```
conda install --yes -c conda-forge -c bioconda quicktree r-base bioconductor-treeio r-ggpubr bioconductor-ggtree r-phytools r-dplyr r-argparse r-scales r-stringr
R -e 'install.packages(c("gggenomes"), repos = "https://cran.r-project.org")'
```

### Usage
```
usage: plot_gggenomes.py [-h] --blocks BLOCKS --fais FAIS [FAIS ...] [--name_conversion NAME_CONVERSION] [--tree TREE] [--normalize] [--indel INDEL] [--length LENGTH]
                         [--centromeres CENTROMERES] [--prefix PREFIX] [--format {png,pdf}] [--scale SCALE] [--height HEIGHT] [--width WIDTH]
                         [--ribbon_adjust RIBBON_ADJUST] [-f] [-n]

Generate a ribbon plot to visualize ntSynt synteny blocks

optional arguments:
  -h, --help            show this help message and exit
  --blocks BLOCKS       ntSynt synteny blocks TSV
  --fais FAIS [FAIS ...]
                        FAI files for all input assemblies. Can be a list or a file with one FAI path per line.
  --name_conversion NAME_CONVERSION
                        TSV for converting names in the blocks TSV (old -> new). IMPORTANT: new names cannot have spaces. If you want to have spaces in the final ribbon
                        plot, use the special character '_'. All underscores in the new name will be converted to spaces.
  --tree TREE           User-input tree file in newick format. If specified, this tree will be plotted next to the output ribbon plot, and used for ordering the
                        assemblies. The names in the newick file must match the new names if --name_conversion is specified, or the genome file names in the synteny blocks
                        input file otherwise. If not specified, the synteny blocks will be used to estimate pairwise distances for the assembly ordering and associated
                        tree.
  --normalize           Normalize strand of chromosomes relative to the target (top) genome in the ribbon plots
  --indel INDEL         Indel size threshold [50000]
  --length LENGTH       Minimum synteny block length [50000]
  --centromeres CENTROMERES
                        TSV file with centromere positions. Must have the headers: bin_id,seq_id,start,end. bin_id must match the new names from --name_conversion or the
                        assembly names if --name_conversion is not specified. seq_id is the chromosome name.
  --prefix PREFIX       Prefix for output files [ntSynt_distance-est]
  --format {png,pdf}    Output format of ribbon plot [png]
  --scale SCALE         Length of scale bar in bases [1e9]
  --height HEIGHT       Height of plot in cm [20]
  --width WIDTH         Width of plot in cm [50]
  --ribbon_adjust RIBBON_ADJUST
                        Ratio for adjusting spacing beside ribbon plot. Increase if ribbon plot labels are cut off, and decrease to reduce the white space to the left of
                        the ribbon plot [0.1]
  -f, --force           Force a re-run of the entire pipeline
  -n                    Dry-run for snakemake pipeline

```
#### Example commands
All the files referenced in these commands can be found in the `tests` subfolder for you to use in testing.

##### Plot ribbon plots with an input cladogram in newick format, normalizing the strands of the assembly chromosomes
```
plot_gggenomes.py --blocks great-apes.ntSynt.synteny_blocks.tsv --fais fais.tsv --tree great-apes.mt-tree.nwk --name_conversion great-apes.name-conversions.tsv --normalize --prefix great-apes_ribbon-plots --ribbon_adjust 0.14
```
![Example_ribbon_plot](https://github.com/bcgsc/ntSynt/blob/distance_est/visualization_scripts/tests/great-apes_ribbon-plots.example1.png)

##### Plot ribbon plots without input cladogram, skipping normalization of the assembly chromosome strands, and changing scale size
```
plot_gggenomes.py --blocks great-apes.ntSynt.synteny_blocks.tsv --fais fais.tsv  --name_conversion great-apes.name-conversions.tsv  --prefix great-apes_ribbon-plots_no-tree --ribbon_adjust 0.15 --scale 500000000 
```
![Example_ribbon_plot](https://github.com/bcgsc/ntSynt/blob/distance_est/visualization_scripts/tests/great-apes_ribbon-plots.example2.png)

#### Other visualization options
For an example script for the chromosome painting plot shown in the ntSynt paper, see chromosome_painting_plots subdirectory

