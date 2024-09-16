#!/usr/bin/env python3
'''
Format input FAI files and synteny TSV for visualization with gggenomes in R
'''
import argparse
import re
from collections import namedtuple, defaultdict
import os

re_fai = re.compile(r'^(\S+).fai$')
SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand"])

REVERSE_COMP_CHAR = "*"
PLACEHOLDER_CHAR = "_"

def read_orientations_file(orientations):
    "Read in file with relative orientations"
    orientations_dict = defaultdict(dict)
    with open(orientations, 'r', encoding="utf-8") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            genome, chrom, relative_ori = line.strip().split("\t")
            orientations_dict[genome][chrom] = REVERSE_COMP_CHAR if relative_ori == "-" else PLACEHOLDER_CHAR
    return orientations_dict

def read_assembly_block_pairs(block_list):
    "Read a list of assembly:chromosome pairs, and load into a dictionary"
    block_dict = {}
    if not block_list:
        return block_dict

    for block in block_list:
        assembly, block = block.strip().split(":")
        if assembly not in block_dict:
            block_dict[assembly] = set()
        block_dict[assembly].add(block)

    return block_dict

def read_name_conversions(name_conversion_fin):
    "Read in the TSV file with name conversions"
    name_conversion_dict = {}
    with open(name_conversion_fin, 'r', encoding="utf-8") as fin:
        for line in fin:
            old_name, new_name = line.strip().split("\t")
            name_conversion_dict[old_name] = new_name
    return name_conversion_dict

def valid_print_sequence(genome_name, chrom_name, keep):
    "Check if a sequence should be printed"
    if not keep:
        return True
    if keep and genome_name in keep and chrom_name in keep[genome_name]:
        return True
    return False

def read_fais(fai_files, name_conversion_dict=None):
    "Read in the FAI file"
    fais = defaultdict(dict) # genome_name -> chrom_name -> length
    for fai_file in fai_files:
        with open(fai_file, 'r', encoding="utf-8") as fai_in:
            for line in fai_in:
                line = line.strip().split('\t')
                chrom_name, length = line[:2]
                file_basename = os.path.basename(fai_file)
                if name_conversion_dict:
                    genome_name = name_conversion_dict[re.search(re_fai, file_basename).group(1)]
                else:
                    genome_name = re.search(re_fai, file_basename).group(1)
                fais[genome_name][chrom_name] = int(length)
    return fais

def make_sequence_file(fais, prefix, orientations=None, keep=None):
    "Generate a sequence length TSV from the input FAI files"
    if orientations:
        orientations_dict = read_orientations_file(orientations) # assembly -> chr -> ori
    else:
        orientations_dict = None

    with open(f"{prefix}.sequence_lengths.tsv", 'w', encoding="utf-8") as fout:
        fout.write("bin_id\tseq_id\tlength\trelative_orientation\n")
        for genome_name in fais:
            for chrom_name, length in fais[genome_name].items():
                if valid_print_sequence(genome_name, chrom_name, keep):
                    relative_ori = f"\t{orientations_dict[genome_name][chrom_name]}" \
                        if orientations_dict and genome_name in orientations_dict else "\t_"
                    fout.write(f"{genome_name}\t{chrom_name}\t{length}{relative_ori}\n")

def make_links_file(synteny_file, prefix, valid_blocks_set, target_assembly):
    "Generate the links TSV from the input synteny blocks file"
    prev_line = None
    lines_to_print = []
    target_assembly_chrom = None
    with open(f"{prefix}.links.tsv", 'w', encoding="utf-8") as fout:
        fout.write("block_id\tseq_id\tbin_id\tstart\tend\t"\
                    "seq_id2\tbin_id2\tstart2\tend2\tstrand\tcolour_block\n")
        with open(synteny_file, 'r', encoding="utf-8") as fin:
            for line in fin:
                line = line.strip().split("\t")
                curr_block = SyntenyBlock(*line[:6])
                if prev_line is not None and prev_line.id == curr_block.id:
                    start_prev, end_prev = prev_line.start, prev_line.end
                    start_curr, end_curr = curr_block.start, curr_block.end
                    relative_ori = "-" if curr_block.strand != prev_line.strand else "+"
                    out_line = f"{curr_block.id}\t{prev_line.chrom}\t{prev_line.genome}\t{start_prev}\t{end_prev}\t"\
                                f"{curr_block.chrom}\t{curr_block.genome}\t{start_curr}\t{end_curr}\t{relative_ori}"
                    lines_to_print.append(out_line)

                if prev_line is not None and prev_line.id != curr_block.id:
                    if prev_line.id in valid_blocks_set:
                        for print_line in lines_to_print:
                            fout.write(f"{print_line}\t{target_assembly_chrom}\n")
                    lines_to_print = []
                if curr_block.genome == target_assembly:
                    target_assembly_chrom = curr_block.chrom
                prev_line = curr_block

            if prev_line.id in valid_blocks_set:
                for print_line in lines_to_print:
                    fout.write(f"{print_line}\t{target_assembly_chrom}\n")

def is_over_block_length_threshold(block, length_threshold):
    "Return True if the block is longer than the threshold, else False"
    return int(block.end) - int(block.start) > length_threshold

def is_over_seq_length_threshold(block, seq_length_threshold, fais):
    "Return True if the block is from a sequence longer than the threshold, else False"
    return fais[block.genome][block.chrom] > seq_length_threshold

def is_keep_block(synteny_block, keep_dict):
    "Return True if the block coordinates are in the keep dictionary or keep is empty, else False"
    if keep_dict:
        if synteny_block.genome in keep_dict and synteny_block.chrom in keep_dict[synteny_block.genome]:
            return True
        return False
    return True

def find_valid_blocks(block_filename, length_threshold, keep_dict, fais, seq_length_threshold):
    "Return set of valid block IDs, where all extents are longer than the threshold and respect the keep lists"
    omit_blocks = set()
    keep_blocks = set()
    all_blocks = set()
    with open(block_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            curr_block = SyntenyBlock(*line[:6])
            if not is_over_block_length_threshold(curr_block, length_threshold) or \
                not is_over_seq_length_threshold(curr_block, seq_length_threshold, fais):
                omit_blocks.add(curr_block.id)
            if is_keep_block(curr_block, keep_dict):
                keep_blocks.add(curr_block.id)
            all_blocks.add(curr_block)
    keep_blocks = keep_blocks - omit_blocks

    keep_seqs = {}
    for block in all_blocks:
        if block.id not in keep_blocks:
            continue
        if block.genome not in keep_seqs:
            keep_seqs[block.genome] = set()
        keep_seqs[block.genome].add(block.chrom)
    if not keep_seqs or not keep_blocks:
        raise ValueError("No valid sequences or links found - check that you are "
                         "using the correct assembly and chromosome naming in --keep.")

    return keep_blocks, keep_seqs


def main():
    "Parse the input files for ntSynt synteny visualization"
    parser = argparse.ArgumentParser(description="Format ntSynt synteny files for visualization in R with gggenomes")
    parser.add_argument("--fai", help="FAI files for genomes", required=True, nargs="+")
    parser.add_argument("--blocks", help="Synteny blocks from ntSynt (Must be in ntSynt format)",
                        required=True, type=str)
    parser.add_argument("-p", "--prefix", help="Prefix for output files [ntsynt_synteny_visuals]",
                        default="ntsynt_synteny_visuals", required=False, type=str)
    parser.add_argument("-l", "--length", help="Minimum block length (bp) [10000]",
                        required=False, type=int, default=10000)
    parser.add_argument("--colour", help="Add chromosome of specified assembly to a final column",
                        required=False, type=str)
    parser.add_argument("--name_conversion", help="Specified file with old -> new names in TSV format",
                        required=False, type=str)
    parser.add_argument("--orientations", help="File with relative orientations to target assembly in TSV format",
                        required=False, type=str)
    parser.add_argument("--keep", help="assembly:chromosome pairs to keep for ribbon plot output",
                        nargs="+", required=False, type=str)
    parser.add_argument("-s", "--seq-length", help="Minimum assembly sequence length [500000]",
                        required=False, type=float, default=500000)
    args = parser.parse_args()

    if args.name_conversion:
        name_conversion_dict = read_name_conversions(args.name_conversion)

    fais = read_fais(args.fai, name_conversion_dict if args.name_conversion else None)

    keep_blocks = read_assembly_block_pairs(args.keep)

    valid_block_ids, keep_seqs = find_valid_blocks(args.blocks, args.length, keep_blocks, fais, args.seq_length)

    colour_assembly = args.colour if args.colour else re.search(r'^(\S+).fai$', args.fai[0]).group(1)

    if args.name_conversion:
        make_sequence_file(fais, args.prefix, args.orientations, keep=keep_seqs)
    else:
        make_sequence_file(fais, args.prefix, orientations=args.orientations, keep=keep_seqs)

    make_links_file(args.blocks, args.prefix, valid_block_ids, colour_assembly)

if __name__ == "__main__":
    main()
