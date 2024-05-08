#!/usr/bin/env python3
'''
Estimates pairwise distances given a synteny block file from ntSynt
'''
import sys
from collections import defaultdict
import itertools



class SyntenyBlock:
    "Represents a synteny block"
    def __init__(self, id, genome, chrom, start, end, strand, num_mx, reason):
        self.id = int(id)
        self.genome = genome
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.num_mx = int(num_mx)
        self.reason = reason

def get_block_i(block_id, block_ids):
    "Generate the new block ID"
    if block_id not in block_ids:
            block_ids[block_id] = len(block_ids)
    return block_ids[block_id]

def load_blocks(blocks_filename):
    "Read in the synteny blocks"
    blocks = {} # key: synteny block id, value: dictionary with (key: assembly; value: SyntenyBlock)
    block_ids = set() # set of block ids
    assemblies = set()
    block_i_ids = {}
    with open(blocks_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            block = SyntenyBlock(*line)
            block_ids.add(block.id)
            block_i = get_block_i(block.id, block_i_ids)
            if block_i not in blocks:
                blocks[block_i] = {}
            blocks[block_i][block.genome] = block
            assemblies.add(block.genome)
    return blocks, len(block_ids), list(assemblies)

def get_difference_between_blocks(blocki, blockj):
    "Returns the gap between the blocks on the assembly, as well as whether the positions are increasing"
    if blockj.start > blocki.start: # Increasing positions
        return min(blockj.start, blockj.end) - max(blocki.start, blocki.end), False
    return min(blocki.start, blocki.end) - max(blockj.start, blockj.end), True


def asm_blocks_consistent(block_i1, block_i2, block_j1, block_j2, indel_threshold):
    "Check if the two blocks from the  assemblies are consistent"
    asm1_transition = block_i1.strand == block_j1.strand
    asm2_transition = block_i2.strand == block_j2.strand
    if asm1_transition != asm2_transition:
        return False
    len_diff_i, is_decreasing_i = get_difference_between_blocks(block_i1, block_j1)
    len_diff_j, is_decreasing_j = get_difference_between_blocks(block_i2, block_j2)
    valid_transition = False
    if block_i1.strand == "+" and block_i2.strand == "+" and not is_decreasing_i and not is_decreasing_j:
        valid_transition = True
    elif block_i1.strand == "-" and block_i2.strand == "-" and is_decreasing_i and is_decreasing_j:
        valid_transition = True
    elif block_i1.strand == "+" and block_i2.strand == "-" and not is_decreasing_i and is_decreasing_j:
        valid_transition = True
    elif block_i1.strand == "-" and block_i2.strand == "+" and is_decreasing_i and not is_decreasing_j:
        valid_transition = True
    if not valid_transition:
        return False    

    if abs(len_diff_j - len_diff_i) > indel_threshold:
        return False
    return True


def compare_block_consistencies(blocks, num_blocks, assembly_pairs, indel_threshold):
    "Given the synteny blocks, compare consecutive pairs to see if they are consistent"
    i = 0
    assembly_distances = defaultdict(dict)
    for asm1, asm2 in assembly_pairs:
        assembly_distances[asm1][asm2] = []

    while i < num_blocks - 1:
        j = i + 1
        for asm1, asm2 in assembly_pairs: 
            if asm_blocks_consistent(blocks[i][asm1], blocks[i][asm2],
                                     blocks[j][asm1], blocks[j][asm2], indel_threshold):
                assembly_distances[asm1][asm2].append(True)
            else:
                assembly_distances[asm1][asm2].append(False)
        i += 1

    return assembly_distances

def calculate_distance(pairwise_estimates):
    "Calculate the distance based on the consistency bit vector"
    bit_dist = 1 - (sum(1 for is_consistent in pairwise_estimates if is_consistent)/len(pairwise_estimates))
    return bit_dist


def main():
    "Estimate pairwise distances from ntSynt blocks"
    if len(sys.argv[1:]) != 2:
        print(f'Usage: {sys.argv[0]} <Synteny blocks TSV> <Indel threshold>', file=sys.stderr, flush=True)
        sys.exit(1)

    blocks_filename = sys.argv[1]
    indel_threshold = int(sys.argv[2])
    blocks, num_blocks, assemblies = load_blocks(blocks_filename)

    all_assembly_pairs = sorted([(asm1, asm2) if asm1 < asm2 else (asm2, asm1) \
                                for asm1, asm2 in itertools.combinations(assemblies, 2)])

    pairwise_vectors = compare_block_consistencies(blocks, num_blocks, all_assembly_pairs, indel_threshold)

    print("asm1", "asm2", "est_distance", sep="\t")
    # Print out each pairwise vector list with the assembly names
    for asm1, asm2 in all_assembly_pairs:
        edit_distance = calculate_distance(pairwise_vectors[asm1][asm2])
        print(asm1, asm2, round(edit_distance, 7),
              sep="\t")


if __name__ == "__main__":
    main()
