#!/bin/bash

# Format the input synteny blocks for visualization with gggenomes

if [ $# -lt 6 ]; then
    echo "Usage: $(basename $0) <synteny blocks TSV> <prefix> <length threshold> <assembly to use for colour> <FAI> <FAI> [FAI..]"
    echo "NOTE: The order of FAI files will dictate the order of the genomes in the ribbon plot"
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

# Sort the sequences based on the synteny block mappings
target_fai=$(echo ${fais} |tr " " "\n" |grep "${target_colour}")
gggenomes_sort_sequences.py --target ${target_colour} --fai ${target_fai} --blocks ${synteny_tsv} \
    --lengths ${prefix}.sequence_lengths.tsv    > ${prefix}.sequence_lengths.updated.tsv
mv ${prefix}.sequence_lengths.updated.tsv ${prefix}.sequence_lengths.tsv 
