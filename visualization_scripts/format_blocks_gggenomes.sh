#!/bin/bash

# Format the input synteny blocks for visualization with gggenomes

if [ $# -lt 6 ]; then
    echo "Usage: $(basename $0) <synteny blocks TSV> <prefix> <length threshold> <assembly to use for colour> <FAI> <FAI> [FAI..]"
    echo "NOTE: The order of FAI files specified will dictate the order of the genomes in the ribbon plot"
    exit 1
fi

synteny_tsv=$1; shift
prefix=$1; shift
length_threshold=$1; shift
target_colour=$1; shift
fais=$@

set -eux -o pipefail

# Sort the blocks based on the specified order
sort_ntsynt_blocks.py --synteny_blocks ${synteny_tsv} --sort_order ${fais} --fais > ${prefix}.synteny_blocks.sorted.tsv

# Generate the files needed for plotting with gggenomes
format_blocks_gggenomes.py --fai ${fais} --prefix ${prefix} --blocks ${prefix}.synteny_blocks.sorted.tsv \
    --length ${length_threshold} --colour ${target_colour}
cat ${prefix}.links.tsv  |mlr --tsv sort -f strand -n block_id > ${prefix}.links.sorted.tsv && mv ${prefix}.links.sorted.tsv ${prefix}.links.tsv
cat ${prefix}.sequence_lengths.tsv |mlr --tsv sort -f seq_id > ${prefix}.sequence_lengths.sorted.tsv && mv ${prefix}.sequence_lengths.sorted.tsv ${prefix}.sequence_lengths.tsv