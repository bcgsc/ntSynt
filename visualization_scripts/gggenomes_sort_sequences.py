#!/usr/bin/env python3
'''
Use synteny block mappings to sort the sequences in the assemblies
'''
import argparse
import itertools
from collections import namedtuple, defaultdict
import intervaltree

SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand", "num_mx", "reason"])
MapRegion = namedtuple("MapRegion", ["chrom", "length"])

def insert_interval_tree(tree, blocks) -> None:
    "Insert block info into the intervaltree"
    for block_target, block_next in zip(blocks, blocks[1:]):
        if block_target.genome not in tree:
            tree[block_target.genome] = {}
        if block_target.chrom not in tree[block_target.genome]:
            tree[block_target.genome][block_target.chrom] = intervaltree.IntervalTree()
        block_data = MapRegion(block_next.chrom, int(block_next.end) - int(block_next.start))
        tree[block_target.genome][block_target.chrom][int(block_target.start): int(block_target.end)] = block_data


def load_synteny_blocks(blocks_file):
    "Load the synteny blocks into an interval tree"
    blocks_tree = {}

    cur_block_id = None
    synteny_blocks = []
    with open(blocks_file, 'r', encoding="utf-8") as fin:
        for line in fin:
            block = SyntenyBlock(*line.strip().split("\t"))
            if block.id != cur_block_id:
                if cur_block_id is not None:
                    # At the next block. insert info into intervaltree
                    insert_interval_tree(blocks_tree, synteny_blocks)
                # Reset variables
                cur_block_id = block.id
                synteny_blocks = []
            synteny_blocks.append(block)

    # Don't forget last one
    if cur_block_id is not None:
        # At the next block. insert info into intervaltree
        insert_interval_tree(blocks_tree, synteny_blocks)
    asm_orders = [block.genome for block in synteny_blocks]

    return blocks_tree, asm_orders

def tally_chromosome_hits(tree, asm, chrom, length, tile):
    "Tally chromosome hits per assembly"
    tiles = []
    for start in range(0, length, tile):
        asm_tallies = defaultdict(defaultdict) # chrom -> length
        end = min(start + tile, length)
        for interval_hit in tree[asm][chrom][start: end]:
            map_region = interval_hit.data
            if map_region.chrom not in asm_tallies:
                asm_tallies[map_region.chrom] = 0
            asm_tallies[map_region.chrom] += map_region.length
        if asm_tallies:
            best_hit_chrom, best_hit_length = sorted([(chrom, length) for chrom, length in asm_tallies.items()],
                                                    key=lambda x: x[1],
                                                    reverse=True)[0]

            tiles.append((best_hit_chrom, best_hit_length))
    return tiles

def get_chrom_orders(tile_list):
    "Given the tile list for a given assembly, return the chromosome orders"  
    collapsed_list = []
    for chrom, group in itertools.groupby(tile_list, key=lambda x: x[0]):
        total = sum(c[1] for c in group)
        collapsed_list.append((chrom, total))
    seq_order = []
    max_indexes = {} # chrom -> (max_index, max_val)
    for i, tup in enumerate(collapsed_list):
        chrom, length = tup
        if chrom not in max_indexes or length > max_indexes[chrom][1]:
            max_indexes[chrom] = (i, length)
    seq_order = [tile_tup[0] for i, tile_tup in enumerate(collapsed_list) if max_indexes[tile_tup[0]][0] == i]
    seq_order = {chrom: i for i, chrom in enumerate(seq_order)}
    return seq_order

def get_mapped_tiles(tree, fai_filenames, tile, asm_orders):
    "Get an ordered list of chrom, length per assembly using tiles over the target"
    all_asm_seq_orders = {} # asm -> chrom -> index
    fai_num = 0
    for fai_filename in fai_filenames:
        if fai_num >= len(asm_orders) - 1:
            break
        asm_base, asm_target = asm_orders[fai_num], asm_orders[fai_num + 1]
        tiles = [] # [list of chromosome hits]
        stored_fai_lines = []
        with open(fai_filename, 'r', encoding="utf-8") as fin:
            for line in fin:
                chrom, length = line.strip().split("\t")[:2]
                length = int(length)
                stored_fai_lines.append((chrom, length))
                if asm_base in all_asm_seq_orders and chrom not in all_asm_seq_orders[asm_base]:
                    all_asm_seq_orders[asm_base][chrom] = len(all_asm_seq_orders[asm_base])
        if asm_base in all_asm_seq_orders:
            stored_fai_lines = sorted(stored_fai_lines, key=lambda x: all_asm_seq_orders[asm_base][x[0]])
        for chrom, length in stored_fai_lines:
            tiles += tally_chromosome_hits(tree, asm_base, chrom, length, tile)

        all_asm_seq_orders[asm_target] = get_chrom_orders(tiles)
        fai_num += 1

    return all_asm_seq_orders


def main():
    "Sort the sequences based on tiles on the target assembly"
    parser = argparse.ArgumentParser(description="Sort chromosomes based on tiles in target assembly")
    parser.add_argument("--fais", help="FAI files for assemblies, in same sort order as --blocks", required=True,
                        type=str, nargs="+")
    parser.add_argument("--blocks", help="Synteny blocks TSV - sorted by order for ribbon plots (if applicable)",
                        required=True, type=str)
    parser.add_argument("--tile", help="Tile size in bp [1 Mbp]", default=1000000, type=int)
    parser.add_argument("--lengths", help="Sequences lengths gggenomes TSV", required=True, type=str)

    args = parser.parse_args()

    blocks_tree, asm_orders = load_synteny_blocks(args.blocks)

    asm_seq_orders = get_mapped_tiles(blocks_tree, args.fais, args.tile, asm_orders)

    stored_lines = {} # asm -> stored_lines
    with open(args.lengths, 'r', encoding='utf-8') as fin:
        for line in fin:
            asm_name, chrom, length = line.strip().split("\t")
            if asm_name == "bin_id" or asm_name == asm_orders[0]:
                print(asm_name, chrom, length, sep="\t")
            else:
                if asm_name not in stored_lines:
                    stored_lines[asm_name] = []
                stored_lines[asm_name].append((asm_name, chrom, length))
                if chrom not in asm_seq_orders[asm_name]:
                    asm_seq_orders[asm_name][chrom] = len(asm_seq_orders[asm_name])

    for asm, lines_list in stored_lines.items():
        asm_orders_asm = asm_seq_orders[asm]
        for line in sorted(lines_list, key=lambda x: asm_orders_asm[x[1]]):
            print(*line, sep="\t")


if __name__ == "__main__":
    main()
