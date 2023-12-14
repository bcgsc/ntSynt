#!/bin/bash

# Format the input synteny blocks for visualization with gggenomes

if [ $# -lt 6 ]; then
    echo "Usage: $(basename $0) <synteny blocks TSV> <prefix> <length threshold> <assembly to use for colour> <FAI> <FAI> [FAI..]"
    exit 1
fi

synteny_tsv=$1; shift
prefix=$1; shift
length_threshold=$1; shift
target_colour=$1; shift
fais=$@

set -eux -o pipefail

format_for_visualization.py --fai ${fais} --prefix ${prefix} --blocks ${synteny_tsv} \
    --length ${length_threshold} --colour ${target_colour}
cat ${prefix}.links.tsv  |mlr --tsv sort -f strand -n block_id -f bin_id > ${prefix}.links.sorted.tsv && mv ${prefix}.links.sorted.tsv ${prefix}.links.tsv
cat ${prefix}.sequence_lengths.tsv |mlr --tsv sort -f bin_id,seq_id > ${prefix}.sequence_lengths.sorted.tsv && mv ${prefix}.sequence_lengths.sorted.tsv ${prefix}.sequence_lengths.tsv