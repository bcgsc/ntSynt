#!/usr/bin/env python3
'''
Format synteny blocks TSV (in ntSynt format) for chromosome painting plots
'''
import argparse
from collections import namedtuple

SyntenyBlock = namedtuple("SyntenyBlock", ["assembly", "chrom", "start", "end", "orientation"])

def read_conversions(conversion_filename):
    "Read in the desired conversions to a dictionary"
    convert_dict = {}
    with open(conversion_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            asm, new_name = line.strip().split("\t")
            convert_dict[asm] = new_name
    return convert_dict

def make_reformat_file(synteny_filename, target, convert=None):
    "Print out reformatted file"
    curr_block_id = "0"
    curr_blocks = []
    target_block = None

    print("block_id", "target_species", "target_chrom", "target_start", "target_end",
          "relative_ori", "other_species", "other_chrom", "other_start", "other_end", sep="\t")

    with open(synteny_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            block_id, asm, chrom, start, end, ori = line[:6]
            new_block = SyntenyBlock(asm, chrom, start, end, ori)
            if block_id != curr_block_id: # We have read in all the blocks
                if target_block is not None:
                    target_assembly = target_block.assembly if convert is None else convert[target_block.assembly]

                    for other_block in curr_blocks:
                        # Write a synteny block line for each 'other' species relative to the target species
                        other_assembly = other_block.assembly if convert is None else convert[other_block.assembly]
                        new_ori = "+" if other_block.orientation == target_block.orientation else "-"
                        print(curr_block_id, target_assembly, target_block.chrom, target_block.start, target_block.end,
                            new_ori, other_assembly, other_block.chrom, other_block.start, other_block.end, sep="\t")

                curr_block_id = block_id
                target_block = None
                curr_blocks = []
            # Store information about the new block
            if asm == target:
                target_block = new_block
            else:
                curr_blocks.append(new_block)

    # Don't forget last synteny block in the file
    target_assembly = target_block.assembly if convert is None else convert[target_block.assembly]

    for other_block in curr_blocks:
        # Write a synteny block line for each 'other' species relative to the target species
        other_assembly = other_block.assembly if convert is None else convert[other_block.assembly]
        ori = "+" if other_block.orientation == target_block.orientation else "-"
        print(curr_block_id, target_assembly, target_block.chrom, target_block.start, target_block.end,
                ori, other_assembly, other_block.chrom, other_block.start, other_block.end, sep="\t")



def main():
    "Reformat the synteny blocks TSV for chromosome sequence painting plots"
    parser = argparse.ArgumentParser(description="Formatting synteny blocks for chromosome painting")
    parser.add_argument("synteny_tsv", help="ntSynt-formatted synteny blocks TSV")
    parser.add_argument("--convert", help="TSV file with desired conversions for assembly names", required=False)
    parser.add_argument("--target", help="Target assembly name", required=True)
    args = parser.parse_args()

    if args.convert:
        conversion_dict = read_conversions(args.convert)
        make_reformat_file(args.synteny_tsv, args.target, convert=conversion_dict)
    else:
        make_reformat_file(args.synteny_tsv, args.target)


if __name__ == "__main__":
    main()
