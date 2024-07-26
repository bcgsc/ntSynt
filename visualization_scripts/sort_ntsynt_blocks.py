#!/usr/bin/env python3
'''
Sort the assemblies within each synteny block in the specified order.
Assumes ntSynt-formatted synteny blocks
'''
import argparse
import os
import re
from collections import namedtuple

SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand", "num_mx", "reason"])
SyntenyBlockComp = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand"])

faidx_re = re.compile(r'^(\S+)\.fai$')

def sort_blocks(synteny_blocks_fin, sorting_order):
    "Sort the assemblies within each synteny block"

    synteny_block_lines = []
    curr_block_id = None

    with open(synteny_blocks_fin, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            if len(line) == 8:
                block = SyntenyBlock(*line)
            else:
                block = SyntenyBlockComp(*line[:6])
            if curr_block_id is not None and block.id != curr_block_id:
                # Print out the sorted blocks, they are all read in
                for asm_block in sorted(synteny_block_lines, key=lambda x: sorting_order[x.genome]):
                    print(*asm_block, sep="\t")
                synteny_block_lines = []
            synteny_block_lines.append(block)
            curr_block_id = block.id

    # Print out the sorted blocks, they are all read in (making sure we have last block)
    for asm_block in sorted(synteny_block_lines, key=lambda x: sorting_order[x.genome]):
        print(*asm_block, sep="\t")


def main():
    "Sort the input ntSynt synteny blocks"
    parser = argparse.ArgumentParser(description="Sort the assemblies in the ntSynt synteny blocks in specified order")
    parser.add_argument("--synteny_blocks", help="Input synteny blocks", required=True, type=str)
    parser.add_argument("--sort_order", help="Desired assembly sort order", nargs="+",
                        required=True)
    parser.add_argument("--fais", help="The assembly sort order option lists the FAI files for the assemblies",
                        action="store_true")
    args = parser.parse_args()

    if args.fais:
        sort_order_dict = {re.search(faidx_re,os.path.basename(os.path.realpath(asm))).group(1): i
            for i, asm in enumerate(args.sort_order)
        }
    else:
        sort_order_dict = {asm: i for i, asm in enumerate(args.sort_order)}

    sort_blocks(args.synteny_blocks, sort_order_dict)


if __name__ == "__main__":
    main()
