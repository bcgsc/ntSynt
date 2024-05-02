#!/usr/bin/env python3
'''
Format input FAI files and synteny TSV for visualization with gggenomes in R
'''
import argparse
import re
from collections import namedtuple
import os

re_fai = re.compile(r'^(\S+).fai$')
SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand"])

def make_sequence_file(fai_files, prefix, name_conversion_file=None):
    "Generate a sequence length TSV from the input FAI files"
    if name_conversion_file:
        name_maps = {}
        with open(name_conversion_file, 'r', encoding="utf-8") as fin:
            for line in fin:
                line = line.strip().split('\t')
                old, new = line
                name_maps[old] = new

    with open(f"{prefix}.sequence_lengths.tsv", 'w', encoding="utf-8") as fout:
        fout.write("bin_id\tseq_id\tlength\n")
        for fai_file in fai_files:
            with open(fai_file, 'r', encoding="utf-8") as fai_in:
                for line in fai_in:
                    line = line.strip().split('\t')
                    chrom_name, length = line[:2]
                    file_basename = os.path.basename(fai_file)
                    if name_conversion_file:
                        fout.write(f"{name_maps[re.search(re_fai, file_basename).group(1)]}\t{chrom_name}\t{length}\n")
                    else:
                        fout.write(f"{re.search(re_fai, file_basename).group(1)}\t{chrom_name}\t{length}\n")

def make_links_file(synteny_file, prefix, valid_blocks_set, target_assembly):
    "Generate the links TSV from the input synteny blocks file"
    prev_line = None
    curr_block_inv = False
    lines_to_print = []
    target_assembly_chrom = None
    with open(f"{prefix}.links.tsv", 'w', encoding="utf-8") as fout:
        fout.write("block_id\tseq_id\tbin_id\tstart\tend\t"\
                    "seq_id2\tbin_id2\tstart2\tend2\tstrand\tblock_ori\tcolour_block\n")
        with open(synteny_file, 'r', encoding="utf-8") as fin:
            for line in fin:
                line = line.strip().split("\t")
                curr_block = SyntenyBlock(*line[:6])
                if prev_line is not None and prev_line.id == curr_block.id:
                    start_prev, end_prev = prev_line.start, prev_line.end
                    curr_block_inv = True if curr_block.strand == "-" else curr_block_inv
                    relative_ori = "-" if curr_block.strand != prev_line.strand else "+"
                    start_curr, end_curr = curr_block.start, curr_block.end
                    out_line = f"{curr_block.id}\t{prev_line.chrom}\t{prev_line.genome}\t{start_prev}\t{end_prev}\t"\
                                f"{curr_block.chrom}\t{curr_block.genome}\t{start_curr}\t{end_curr}\t{relative_ori}"
                    lines_to_print.append(out_line)

                if prev_line is not None and prev_line.id != curr_block.id:
                    block_type = "-" if curr_block_inv else "+"
                    if prev_line.id in valid_blocks_set:
                        for print_line in lines_to_print:
                            fout.write(f"{print_line}\t{block_type}\t{target_assembly_chrom}\n")
                    lines_to_print = []
                    curr_block_inv = False
                if curr_block.genome == target_assembly:
                    target_assembly_chrom = curr_block.chrom
                prev_line = curr_block

            if prev_line.id in valid_blocks_set:
                for print_line in lines_to_print:
                    fout.write(f"{print_line}\t{block_type}\t{target_assembly_chrom}\n")

def find_valid_block_ids(block_filename, length_threshold):
    "Return set of valid block IDs, where all extents are longer than the threshold"
    valid_blocks = set()
    with open(block_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            curr_block = SyntenyBlock(*line[:6])
            if int(curr_block.end) - int(curr_block.start) >= length_threshold:
                valid_blocks.add(curr_block.id)
    return valid_blocks


def main():
    "Parse the input files for ntSynt synteny visualization"
    parser = argparse.ArgumentParser(description="Format ntSynt synteny files for visualization in R with gggenomes")
    parser.add_argument("--fai", help="FAI files for genomes", required=True, nargs="+")
    parser.add_argument("--blocks", help="Synteny blocks from ntSynt (Must be in ntSynt format)",
                        required=True, type=str)
    parser.add_argument("-p", "--prefix", help="Prefix for output files [ntsynt_synteny_visuals]",
                        default="ntsynt_synteny_visuals", required=False, type=str)
    parser.add_argument("-l", "--length", help="Minimum block length [10kb]", required=False, type=int, default=10000)
    parser.add_argument("--colour", help="Add chromosome of specified assembly to a final column",
                        required=False, type=str)
    parser.add_argument("--name_conversion", help="Specified file with old -> new names in TSV format",
                        required=False, type=str)
    args = parser.parse_args()

    valid_blocks = find_valid_block_ids(args.blocks, args.length)

    colour_assembly = args.colour if args.colour else re.search(r'^(\S+).fai$', args.fai[0]).group(1)

    make_sequence_file(args.fai, args.prefix, args.name_conversion)
    make_links_file(args.blocks, args.prefix, valid_blocks, colour_assembly)

if __name__ == "__main__":
    main()
