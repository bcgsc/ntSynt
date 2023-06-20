![Logo](https://github.com/bcgsc/ntSynt/blob/main/logo/ntsynt-logo.png)

# ntSynt

Detecting synteny blocks between multiple genome assemblies using a dynamic minimizer graph approach.

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
usage: ntSynt.py [-h] [-k K] [-w W] [-t T] [--fpr FPR] [--no-solid] [--no-simplify-graph] [-p PREFIX] [--merge MERGE] [--w_rounds W_ROUNDS [W_ROUNDS ...]]
                 [--indel INDEL] [--dry-run] [-v]
                 fastas [fastas ...]

ntSynt: Genome synteny detection using dynamic minimizer graphs

positional arguments:
  fastas                Input genome fasta files

optional arguments:
  -h, --help            show this help message and exit
  -k K                  Minimizer kmer size [24]
  -w W                  Minimizer window size [1000]
  -t T                  Number of threads [4]
  --fpr FPR             False positive rate for Bloom filter creation [0.025]
  --no-solid            Do not use the solid BF for minimizer graph creation
  --no-simplify-graph   Do not run graph simplification on minimizer graph
  -p PREFIX, --prefix PREFIX
                        Prefix for ntSynt output files [ntSynt.k<k>.w<w>]
  --merge MERGE         Multiple of window size used for collinear synteny block merging [3]
  --w_rounds W_ROUNDS [W_ROUNDS ...]
                        List of window sizes for iterative rounds [100 10 5]
  --indel INDEL         Threshold for indel detection [500]
  --dry-run             Print out the commands that will be executed
  --benchmark           Store benchmarks for each step of the ntSynt pipeline
  -v, --version         show program's version number and exit
```

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

**All ntSynt dependencies can be installed using conda:**
```
conda install -c bioconda -c conda-forge python intervaltree pybedtools ncls python-igraph compilers btllib meson ninja snakemake samtools
```

### Installing ntSynt from the source code
```
meson build --prefix=/path/to/desired/install/location
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

## License
ntSynt Copyright (c) 2023 British Columbia Cancer Agency Branch. All rights reserved.

ntSynt is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein prebstein@bccancer.bc.ca
