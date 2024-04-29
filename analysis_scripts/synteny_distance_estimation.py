#!/usr/bin/env python3
'''
Estimates pairwise distances given a synteny block file from ntSynt
'''
import sys
from collections import defaultdict, namedtuple
import itertools
import numpy as np

PairConsistency = namedtuple("PairConsistency", ["max_length", "is_consistent", "block_i_avg_length", "sum_block_lengths", "min_length"])

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

MERGE_THRESHOLD = 100000 # TODO: Update

def load_blocks(blocks_filename):
    "Read in the synteny blocks"
    blocks = {} # key: synteny block id, value: dictionary with key: assembly value: SyntenyBlock
    block_ids = set() # set of block ids
    assemblies = set()
    block_mappings = {}
    cur_block = 0
    with open(blocks_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            block = SyntenyBlock(*line)
            if block.id in block_mappings:
                block_id_int = block_mappings[block.id]
            else:
                block_id_int = cur_block
                block_mappings[block.id] = block_id_int
                cur_block += 1
            block_ids.add(block_id_int)
            if block_id_int not in blocks:
                blocks[block_id_int] = {}
            blocks[block_id_int][block.genome] = block
            assemblies.add(block.genome)
    return blocks, len(block_ids), list(assemblies)

def get_difference_between_blocks(block1, block2, diff_strands=True):
    "Returns the gap between the blocks on the assembly"
    if diff_strands and block1.strand == "-" and block2.strand == "-":
        return  block1.start - block2.end
    return block2.start - block1.end #!! Making assumption about block ordering


def asm_blocks_consistent(block_i1, block_i2, block_j1, block_j2, indel_threshold):
    "Check if the two blocks from the  assemblies are consistent"
    asm1_transition = block_i1.strand == block_i2.strand
    asm2_transition = block_j1.strand == block_j2.strand
    if asm1_transition != asm2_transition:
        return False
    # if block_i1.strand != block_j1.strand or block_i2.strand != block_j2.strand:
    #     return False
    diff_strands = block_i1.strand != block_i2.strand
    len_diff_i = get_difference_between_blocks(block_i1, block_j1, diff_strands=diff_strands)
    len_diff_j = get_difference_between_blocks(block_i2, block_j2, diff_strands=diff_strands)
    if abs(len_diff_j - len_diff_i) > indel_threshold:
        return False
    # if max(len_diff_i, len_diff_j) > MERGE_THRESHOLD:
    #     return False
    if len_diff_i < 0 or len_diff_j < 0:
        return False
    return True

def get_block_stats(block_i1, block_j1, block_i2, block_j2):
    "Given the block pair, return block max block length, average length of block i, average of all block lengths"
    max_length = max([block_i1.end - block_i1.start, block_j1.end - block_j1.start,
                     block_i2.end - block_i2.start, block_j2.end - block_j2.start])
    min_length = min([block_i1.end - block_i1.start, block_j1.end - block_j1.start,
                     block_i2.end - block_i2.start, block_j2.end - block_j2.start])
    avg_len = get_block_avg(block_i1, block_i2)
    total_length = sum([block_i1.end - block_i1.start, block_j1.end - block_j1.start,
                        block_i2.end - block_i2.start, block_j2.end - block_j2.start])
    return max_length, avg_len, total_length, min_length

def get_block_avg(block_i1, block_i2):
    "Given the block pair, return the minimum block size"
    avg_length = (block_i1.end - block_i1.start + block_i2.end - block_i2.start)/2
    return avg_length

def compare_block_consistencies(blocks, num_blocks, assembly_pairs, indel_threshold):
    "Given the synteny blocks, compare consecutive pairs to see if they are consistent"
    i = 0
    assembly_distances = defaultdict(dict)
    for asm1, asm2 in assembly_pairs:
        assembly_distances[asm1][asm2] = []

    while i < num_blocks - 1:
        j = i + 1
        for asm1, asm2 in assembly_pairs:
            max_length, avg_len_i, total_length, min_length = get_block_stats(blocks[i][asm1], blocks[i][asm2],
                                             blocks[j][asm1], blocks[j][asm2])            
            if asm_blocks_consistent(blocks[i][asm1], blocks[i][asm2],
                                     blocks[j][asm1], blocks[j][asm2], indel_threshold):
                assembly_distances[asm1][asm2].append(PairConsistency(max_length, True, avg_len_i, total_length, min_length))
                #print(blocks[i][asm1].id, True)
            else:
                assembly_distances[asm1][asm2].append(PairConsistency(max_length, False, avg_len_i, total_length, min_length))
                #print(blocks[i][asm1].id, False)
        i += 1

    return assembly_distances

def calculate_distance(pairwise_estimates):
    "Calculate the distances based on the consistency bit vector"
    dist_len_max = 1 - (sum(max_length for max_length, is_consistent, _, _, _ in pairwise_estimates if is_consistent)/sum(max_length for max_length, _, _, _, _ in pairwise_estimates))
    dist_len_max_unnorm = sum(max_length for max_length, is_consistent, _, _, _ in pairwise_estimates if is_consistent)
    dist_len_min = 1 - (sum(min_length for _, is_consistent, _, _, min_length in pairwise_estimates if is_consistent)/sum(min_length for _, _, _, _, min_length in pairwise_estimates))
    dist_len_min_unnorm = sum(min_length for _, is_consistent, _, _, min_length in pairwise_estimates if is_consistent)
    bit_dist = 1 - (sum(1 for _, is_consistent, _, _, _ in pairwise_estimates if is_consistent)/len(pairwise_estimates))
    bit_vector = [1 if is_consistent else 0 for _, is_consistent, _, _, _ in pairwise_estimates]
    list_runs = [(category, list(g)) for category, g in itertools.groupby(bit_vector)]
    weighted_sums_list = [category * len(list_run) + 0.5 * len(list_run) if category else 0 for category, list_run in list_runs]
    weighted_sums_geometric = sum([category * len(list_run) * len(list_run) if category else 0 for category, list_run in list_runs])
    weighted_sums_norm_denom = sum([len(list_run) * len(list_run) for category, list_run in list_runs])

    return bit_dist, dist_len_max, dist_len_max_unnorm, sum(weighted_sums_list), weighted_sums_geometric, dist_len_min, dist_len_min_unnorm


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

    print("asm1", "asm2", "edit_distance", "edit_distance_len_max", "edit_distance_len_max_unnorm", "edit_distance_len_min", "edit_distance_len_min_unnorm",
          "block_lengths", "weighted_sum", "weighted_sum_geometric", sep="\t")
    # Print out each pairwise vector list with the assembly names
    for asm1, asm2 in all_assembly_pairs:
        edit_distance, edit_distance_len_max, edit_distance_len_max_unnorm, weighted_sum, weighted_sum_geometric, edit_distance_len_min, edit_distance_len_min_unnorm = calculate_distance(pairwise_vectors[asm1][asm2])
        block_lengths = sum(length for _, consistent, length, _, _ in pairwise_vectors[asm1][asm2] if consistent)
        print(asm1, asm2, round(edit_distance, 7), round(edit_distance_len_max, 7), edit_distance_len_max_unnorm,
              edit_distance_len_min, edit_distance_len_min_unnorm,
              block_lengths, weighted_sum, weighted_sum_geometric,
              sep="\t")


if __name__ == "__main__":
    main()
