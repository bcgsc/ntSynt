#!/usr/bin/env python3
'''
Represents an ntSynt synteny block
'''
from collections import namedtuple
import re
from assembly_block import AssemblyBlock
from ntjoin_utils import Minimizer

# NamedTuples
SyntenyBlockNode = namedtuple("SyntenyBlockNode", ["mx", "positions"])

# Regexes
fa_tsv_re = re.compile(r'^(\S+)\.k\d+\.w\d+.tsv')

class SyntenyBlock:
    "A Synteny Block between the input assemblies"
    def __init__(self, k, m, *assemblies):
        "Instantiate a dictionary to keep track of assembly blocks for this synteny block"
        self.assembly_blocks = {assembly: AssemblyBlock(k) for assembly in assemblies}
        self.m = m
        self.broken_reason = None


    def assign_block(self, assembly, assembly_block):
        "Over-write any existing assembly block for the specified assembly"
        if assembly not in self.assembly_blocks:
            raise ValueError(f"{assembly} not found in this Synteny Block!")
        self.assembly_blocks[assembly] = assembly_block

    def continue_block(self, mx, list_mx_info):
        "Given minimizer and preliminary blocks, return if synteny block should extend, else False"
        return all(mx_dict[mx][0] == self.assembly_blocks[assembly].contig_id \
            for assembly, mx_dict in list_mx_info.items())

    def extend_block(self, mx, list_mx_jnfo):
        "Extend the synteny block by extending each assembly block"
        for assembly, mx_dict in list_mx_jnfo.items():
            self.assembly_blocks[assembly].minimizers.append(Minimizer(mx, mx_dict[mx][1]))

    def start_block(self, mx, list_mx_info):
        "Start the new synteny block"
        for assembly, mx_dict in list_mx_info.items():
            ctg, pos = mx_dict[mx]
            self.assembly_blocks[assembly].contig_id = ctg
            self.assembly_blocks[assembly].minimizers.append(Minimizer(mx, int(pos)))

    def determine_orientations(self):
        "Determine the orientations of each assembly block"
        for _, assembly_block in self.assembly_blocks.items():
            positions = [mx.position for mx in assembly_block.minimizers]
            if all(x < y for x, y in zip(positions, positions[1:])):
                assembly_block.ori = "+"
            elif all(x > y for x, y in zip(positions, positions[1:])):
                assembly_block.ori = "-"
            else:
                tally = [x < y for x, y in zip(positions, positions[1:])]
                positive_perc = tally.count(True)/float(len(positions)-1)*100
                negative_perc = 100 - positive_perc
                if positive_perc >= self.m:
                    assembly_block.ori = "+"
                elif negative_perc >= self.m:
                    assembly_block.ori = "-"
                else:
                    assembly_block.ori = "?"


    def all_oriented(self):
        "Return true if all of the assembly blocks in the synteny block are oriented"
        return all(assembly_block.ori in ["+", "-"] for _, assembly_block in self.assembly_blocks.items())

    def get_block_string(self, num, verbose=False):
        "Given the specified synteny block ID, print the synteny blocks"
        return_str = ""
        for assembly, assembly_block in sorted(self.assembly_blocks.items()):
            if assembly_match := re.search(fa_tsv_re, assembly):
                assembly = assembly_match.group(1)
            start_pos = assembly_block.get_block_start()
            end_pos = assembly_block.get_block_end()
            block_string = f"{num}\t{assembly}\t{assembly_block.contig_id}\t{start_pos}" \
                f"\t{end_pos}\t{assembly_block.ori}\t{len(assembly_block.minimizers)}\n"
            if verbose:
                block_string = f"{block_string.strip()}\t{self.broken_reason}\n"
            return_str += block_string
        return return_str

    def get_node(self, i):
        "Given an index into the minimizers, return a synteny block node"
        mxs = set()
        positions = []
        for _, assembly_block in sorted(self.assembly_blocks.items()):
            mx, pos = assembly_block.minimizers[i]
            mxs.add(mx)
            positions.append(pos)
        return SyntenyBlockNode(mxs.pop(), positions)

    def get_number_of_minimizers(self):
        "Return the number of minimizers that comprise this synteny block"
        rep_assembly = list(self.assembly_blocks.keys()).pop()
        return len(self.assembly_blocks[rep_assembly].minimizers)

    def __lt__(self, other):
        lexico_smallest_assembly = sorted(self.assembly_blocks.keys())[0]
        if self.assembly_blocks[lexico_smallest_assembly].contig_id == \
            other.assembly_blocks[lexico_smallest_assembly].contig_id:
            return self.assembly_blocks[lexico_smallest_assembly].get_block_start() < \
                    other.assembly_blocks[lexico_smallest_assembly].get_block_start()
        return self.assembly_blocks[lexico_smallest_assembly].contig_id < \
            other.assembly_blocks[lexico_smallest_assembly].contig_id

    def __eq__(self, other):
        is_equal = True
        for assembly, block in self.assembly_blocks.items():
            if block != other.self.assembly_blocks[assembly]:
                is_equal = False
        return is_equal
