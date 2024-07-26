#!/usr/bin/env python3
'''
Normalize the strands of chromosomes in the ntSynt synteny blocks based on the given assembly order
The normalization will be based on the first listed assembly
'''
import argparse
import os
from collections import namedtuple

SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand", "num_mx", "reason"])

def read_fais(fais, name_convert=False):
    "Read in the FAI files, creating a dictionary with asm -> chr -> length"
    fai_dict = {}
    for fai in fais:
        base_fa = os.path.basename(os.path.realpath(fai)).removesuffix(".fai")
        if name_convert:
            base_fa = name_convert[base_fa]
        fai_dict[base_fa] = {}
        with open(fai, 'r', encoding="utf-8") as fin:
            for line in fin:
                chrom, length = line.strip().split("\t")[:2]
                fai_dict[base_fa][chrom] = int(length)
    return fai_dict


def read_name_conversions(name_conversion_fin):
    "Read in the TSV file with name conversions"
    name_conversion_dict = {}
    with open(name_conversion_fin, 'r', encoding="utf-8") as fin:
        for line in fin:
            old_name, new_name = line.strip().split("\t")
            name_conversion_dict[old_name] = new_name
    return name_conversion_dict

def tally_orientations(blocks_list, ori_dict):
    "Tally the lengths of same or different orientations of the blocks compared to the target (first)"
    target_block = blocks_list[0]
    for block in blocks_list[1:]:
        if block.genome not in ori_dict:
            ori_dict[block.genome] = {}
        if block.chrom not in ori_dict[block.genome]:
            ori_dict[block.genome][block.chrom] = {"+": 0, "-": 0}
        relative_strand = "+" if target_block.strand == block.strand else "-"
        ori_dict[block.genome][block.chrom][relative_strand] += int(block.end) - int(block.start)


def tally_orientation_blocks(synteny_blocks):
    "Tally the lengths of same or different orientations of the blocks compared to the target (first)"
    ori_dict = {} # asm -> chr -> {+: len, -: len}
    cur_block_id = None
    blocks = []
    with open(synteny_blocks, 'r', encoding="utf-8") as fin:
        for line in fin:
            block = SyntenyBlock(*line.strip().split("\t"))
            if cur_block_id != block.id:
                if cur_block_id is not None:
                    # Tally the information needed for all blocks
                    tally_orientations(blocks, ori_dict)
                blocks = [block]
                cur_block_id = block.id
            else:
                blocks.append(block)

    tally_orientations(blocks, ori_dict)
    return ori_dict

def assign_orientations(ori_dict):
    "Assign orientations based on the orientation dict length tallies"
    ori_choices = {} # asm -> chr -> + or -
    for asm in ori_dict:
        if asm not in ori_choices:
            ori_choices[asm] = {}
        for chrom in ori_dict[asm]:
            if ori_dict[asm][chrom]["+"] >= ori_dict[asm][chrom]["-"]:
                ori_choices[asm][chrom] = "+"
            else:
                ori_choices[asm][chrom] = "-"
    return ori_choices

def print_synteny_block(block, out_file):
    "Print the synteny block, adjusting the reason for the orientation if needed"
    out_file.write(f"{block.id}\t{block.genome}\t{block.chrom}\t{block.start}"
                   f"\t{block.end}\t{block.strand}\t{block.num_mx}\t{block.reason}\n")


def reverse_complement_coordinates(block, chrom_length):
    "Reverse complement the coordinates"
    start = chrom_length - int(block.end)
    end = chrom_length - int(block.start)
    strand = "+" if block.strand == "-" else "-"
    return SyntenyBlock(block.id, block.genome, block.chrom, start, end,
                        strand, block.num_mx, block.reason + "_rev")

def normalize_blocks(synteny_blocks_fin, fai_dict, orientations, prefix):
    "Normalize the synteny blocks"
    with open(synteny_blocks_fin, 'r', encoding="utf-8") as fin, \
            open(f"{prefix}.blocks.tsv", 'w', encoding="utf-8") as fout:
        for line in fin:
            block = SyntenyBlock(*line.strip().split("\t"))
            if block.genome not in orientations:
                print_synteny_block(block, fout)
            else:
                if orientations[block.genome][block.chrom] == "-":
                    block = reverse_complement_coordinates(block, fai_dict[block.genome][block.chrom])
                print_synteny_block(block, fout)


def main():
    "Normalize the synteny blocks based on the first (target) assembly"
    parser = argparse.ArgumentParser(description="Normalize the synteny blocks based on the first (target) assembly")
    parser.add_argument("-b", "--blocks", help="Input synteny blocks", required=True, type=str)
    parser.add_argument("-f", "--fais", help="FAI files for assemblies, in same sort order as --blocks", required=True,
                        type=str, nargs="+")
    parser.add_argument("-c", "--name_convert", help="Name conversion TSV", required=False, type=str)
    parser.add_argument("-p", "--prefix", help="Prefix for output files", default="normalize_strands", type=str)
    args = parser.parse_args()

    if args.name_convert:
        name_conversion_dict = read_name_conversions(args.name_convert)
        fai_dict = read_fais(args.fais, name_conversion_dict)
    else:
        fai_dict = read_fais(args.fais)

    orientation_tallies = tally_orientation_blocks(args.blocks)

    orientations = assign_orientations(orientation_tallies)

    with open(f"{args.prefix}.chrom-orientations.tsv", 'w', encoding="utf-8") as fout:
        for asm, chr_dict in fai_dict.items():
            if asm not in orientations:
                fout.write(f"#Orientations relative to {asm}\n")
                fout.write("genome\tchromosome\trelative_orientation\n")
        for asm, chr_dict in orientations.items():
            for chrom in chr_dict:
                fout.write(f"{asm}\t{chrom}\t{orientations[asm][chrom]}\n")

    normalize_blocks(args.blocks, fai_dict, orientations, args.prefix)

if __name__ == "__main__":
    main()
