![GitHub Release](https://img.shields.io/github/v/release/bcgsc/ntsynt)

![Logo](https://github.com/bcgsc/ntSynt/blob/main/logo/ntsynt-logo.png)

# ntSynt

Detecting synteny blocks between multiple genome assemblies using a dynamic minimizer graph approach.

## Contents
1. [Description of the algorithm](#description-of-the-algorithm)
2. [Credits](#credits)
3. [Usage](#usage)
4. [Installation instructions](#installation)
5. [Example](#example-command)
6. [Output files](#output-files)
7. [Asssessment](#basic-assessment-of-synteny-blocks)
8. [Tips](#tips)
9. [License](#license)

## Description of the algorithm

ntSynt can take multiple genome assemblies as input, and will compute synteny blocks that are in common with each of these input assemblies. ntSynt builds on the ntJoin codebase.

**Main steps in the algorithm:**
1. Generate ordered minimizer sketches for each of the input assemblies
    * By default, the minimizer sketches will be generated using a "common" Bloom filter, which contains all the k-mers that are in common between each assembly
    * Only minimizers from k-mers that are found in this common Bloom filter will be included in the output sketch
2. Filter the minimizers to retain only those that are found in all input assemblies and found in single copy in each assembly
3. Build an initial minimizer graph, where the nodes are minimizers and there are edges between adjacent minimizers
4. By default, run graph simplification to remove simple tangles in the graph
5. Find linear paths through the graph, and compute the initial synteny blocks
6. Compute minimizer sketches for regions not covered by the initial synteny blocks using a lower window size, and augment the minimizer graph
7. Repeat the graph simplification and filtering steps, and output refined synteny blocks. Repeat steps 6-7 for each desired lower window size.
8. Merge collinear synteny blocks within the specified range
9. Output the final synteny blocks

## Credits

Concept: Lauren Coombe and Rene Warren

Design and implementation: Lauren Coombe

## Usage

```
usage: ntSynt [-h] -d DIVERGENCE [-p PREFIX] [-k K] [-w W] [-t T] [--fpr FPR] [-b BLOCK_SIZE] [--merge MERGE] [--w_rounds W_ROUNDS [W_ROUNDS ...]]
                 [--indel INDEL] [-n] [--benchmark] [-f] [--dev] [-v]
                 fastas [fastas ...]

ntSynt: Multi-genome synteny detection using minimizer graphs

positional arguments:
  fastas                Input genome fasta files

optional arguments:
  -h, --help            show this help message and exit
  -d DIVERGENCE, --divergence DIVERGENCE
                        Approx. maximum percent sequence divergence between input genomes (Ex. -d 1 for 1% divergence).
                        This will be used to set --indel, --merge, --w_rounds, --block_size
                        See below for set values - You can also set any of those parameters yourself, which will override these settings.
  -p PREFIX, --prefix PREFIX
                        Prefix for ntSynt output files [ntSynt.k<k>.w<w>]
  -k K                  Minimizer k-mer size [24]
  -w W                  Minimizer window size [1000]
  -t T                  Number of threads [12]
  --fpr FPR             False positive rate for Bloom filter creation [0.025]
  -b BLOCK_SIZE, --block_size BLOCK_SIZE
                        Minimum synteny block size (bp)
  --merge MERGE         Maximum distance between collinear synteny blocks for merging (bp). 
                        Can also specify a multiple of the window size (ex. 3w)
  --w_rounds W_ROUNDS [W_ROUNDS ...]
                        List of decreasing window sizes for synteny block refinement
  --indel INDEL         Threshold for indel detection (bp)
  -n, --dry-run         Print out the commands that will be executed
  --benchmark           Store benchmarks for each step of the ntSynt pipeline
  -f, --force           Run all ntSynt steps, regardless of existing output files
  --dev                 Run in developer mode to retain intermediate files, log verbose output
  -v, --version         show program's version number and exit
```
Given the approximate maximum divergence between the supplied genomes, ntSynt will set these default parameters:
|Divergence range|Default parameters|
|----|----|
|< 1%|--block_size 500 --indel 10000 --merge 10000 --w_rounds 100 10|
|1% - 10%|--block_size 1000 --indel 50000 --merge 100000 --w_rounds 250 100|
|>10%|	--block_size 10000 --indel 100000 --merge 1000000 --w_rounds 500 250|

Any of these parameters can be overridden by specifying them in your command. While these settings work generally well for the associated divergence range, we highly recommend customizing them for your particular requirements.

## Installation

### Dependencies:
- python 3.9+ with modules:
  - [intervaltree](https://github.com/chaimleib/intervaltree)
  - [pybedtools](https://daler.github.io/pybedtools/)
  - [ncls](https://github.com/pyranges/ncls)
  - [python-igraph](https://python.igraph.org/en/stable/)
- GCC 6+ or Clang 5+ (with OpenMP and C++17 support)
- [btllib](https://github.com/bcgsc/btllib)
- [meson](https://mesonbuild.com/)
- [ninja](https://ninja-build.org/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/)
- [samtools](http://www.htslib.org/)
- [seqtk](https://github.com/lh3/seqtk)

**All ntSynt dependencies can be installed using conda:**
```
conda install -c bioconda -c conda-forge python intervaltree pybedtools ncls python-igraph compilers btllib meson ninja snakemake samtools
```

### Installing ntSynt from the source code
```
meson setup build --prefix=/path/to/desired/install/location
cd build
ninja install
```

### Testing ntSynt installation
Test your ntSynt installation using our provided demo:
```
cd tests
./run_ntSynt_demo.sh 
```
Once the script has executed successfully, you can compare the output files with those in tests/expected_results

### Example command
To compute the synteny blocks between 3 assemblies (assembly1.fa, assembly2.fa, assembly3.fa) with default parameters, where the maximum sequence divergence among these is ~5%, run:
```
ntSynt -d 5 assembly1.fa assembly2.fa assembly3.fa
```

### Output files
The main output file has the naming scheme `<prefix>.synteny_blocks.tsv`. This contains the synteny blocks computed in a TSV format.

The columns of this output synteny blocks TSV:
1. Synteny block ID - Lines with the same ID are part of the same synteny block
2. Genome file name
3. Genome chromosome/contig
4. Genome start coordinate
5. Genome end coordinate
6. Chromosome/contig strand
7. Number of mapped minimizers in this synteny block
8. Reason for discontinuity with previous synteny block

### Basic assessment of synteny blocks
For a basic summary of the statistics of the computed synteny blocks, you can run use the script `denovo_synteny_block_stats.py` found in `analysis_scripts`:
```
python3 denovo_synteny_block_stats.py -h
usage: denovo_synteny_block_stats.py [-h] --tsv TSV --fai FAI [FAI ...]

Compute de novo stats on synteny blocks

optional arguments:
  -h, --help           show this help message and exit
  --tsv TSV            ntSynt synteny block file
  --fai FAI [FAI ...]  FAI files for the compared genomes
```
More information can be found on our [wiki page](https://github.com/bcgsc/ntSynt/wiki/de-novo-statistics-summary)

### Tips
- To lower the peak memory usage, increase the false positive rate (--fpr) for the constructed Bloom filter
- Customize parameters such as --merge, --indel, --block_size and --w_rounds for your particular input data and research questions
- For example scripts for visualizing the output synteny blocks, see the sub-directory [visualization_scripts](https://github.com/bcgsc/ntSynt/tree/main/visualization_scripts)

## License
ntSynt Copyright (c) 2023 British Columbia Cancer Agency Branch. All rights reserved.

ntSynt is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein prebstein@bccancer.bc.ca
