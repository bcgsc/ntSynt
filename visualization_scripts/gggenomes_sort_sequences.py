#!/usr/bin/env python3
'''
Use syntent block mappings to sort the sequences in the assemblies
'''
import argparse
import intervaltree
from collections import namedtuple, defaultdict
import itertools

SyntenyBlock = namedtuple("SyntenyBlock", ["id", "genome", "chrom", "start", "end", "strand", "num_mx", "reason"])
MapRegion = namedtuple("MapRegion", ["chrom", "length"])

def insert_interval_tree(tree, target_block, other_blocks) -> None:
    "Insert block info into the intervaltree"
    if target_block.chrom not in tree:
        tree[target_block.chrom] = intervaltree.IntervalTree()
    block_data = {}
    for block in other_blocks:
        block_data[block.genome] = MapRegion(block.chrom, int(block.end) - int(block.start))
    tree[target_block.chrom][int(target_block.start): int(target_block.end)] = block_data


def load_synteny_blocks(blocks_file, target_asm):
    "Load the synteny blocks into an interval tree"
    blocks_tree = {}

    cur_block_id = None
    target_synteny_block = None
    other_synteny_blocks = []
    with open(blocks_file, 'r', encoding="utf-8") as fin:
        for line in fin:
            block = SyntenyBlock(*line.strip().split("\t"))
            if block.id != cur_block_id:
                if cur_block_id is not None:
                    # At the next block. insert info into intervaltree
                    insert_interval_tree(blocks_tree, target_synteny_block, other_synteny_blocks)
                # Reset variables
                cur_block_id = block.id
                other_synteny_blocks = []
                target_synteny_block = None
            if block.genome == target_asm:
                target_synteny_block = block
            else:
                other_synteny_blocks.append(block)

    # Don't forget last one
    if cur_block_id is not None:
        # At the next block. insert info into intervaltree
        insert_interval_tree(blocks_tree, target_synteny_block, other_synteny_blocks)

    return blocks_tree

def get_mapped_tiles(tree, fai_filename, tile):
    "Get an ordered list of chrom, length per assembly using tiles over the target"
    tiles = {} # asm -> [list of chromosomes]
    with open(fai_filename, 'r', encoding="utf-8") as fin:
        for line in fin:
            chrom, length = line.strip().split("\t")[:2]
            length = int(length)
            for start in range(0, length, tile):
                asm_tallies = defaultdict(defaultdict) # asm -> chrom -> length
                end = min(start + tile, length)
                for interval_hit in tree[chrom][start: end]:
                    for asm, map_region in interval_hit.data.items():
                        if asm not in asm_tallies or map_region.chrom not in asm_tallies[asm]:
                            asm_tallies[asm][map_region.chrom] = 0
                        asm_tallies[asm][map_region.chrom] += map_region.length
                for asm in asm_tallies:
                    best_hit_chrom, best_hit_length = sorted([(chrom, length) for chrom, length in asm_tallies[asm].items()], key=lambda x: x[1],
                                                            reverse=True)[0]
                    if asm not in tiles:
                        tiles[asm] = []
                    tiles[asm].append((best_hit_chrom, best_hit_length))
    return tiles


def main():
    "Sort the sequences based on tiles on the target assembly"
    parser = argparse.ArgumentParser(description="Sort chromosomes based on tiles in target assembly")
    parser.add_argument("--target", help="Name of target assembly file", required=True, type=str)
    parser.add_argument("--fai", help="FAI file for target assembly", required=True, type=str)
    parser.add_argument("--blocks", help="Synteny blocks TSV", required=True, type=str)
    parser.add_argument("--tile", help="Tile size in bp [1 Mbp]", default=1000000, type=int)
    parser.add_argument("--lengths", help="Sequences lengths gggenomes TSV", required=True, type=str)

    args = parser.parse_args()

    blocks_tree = load_synteny_blocks(args.blocks, args.target)

    tile_list = get_mapped_tiles(blocks_tree, args.fai, args.tile)

    collapsed_list = {}

    for asm, tiles in tile_list.items():
        collapsed_list[asm] = []
        for chrom, group in itertools.groupby(tiles, key=lambda x: x[0]):
            total = sum(c[1] for c in group)
            collapsed_list[asm].append((chrom, total))
    asm_orders = {} # asm -> [sequence order]
    for asm, tiles in collapsed_list.items():
        max_indexes = {} # chrom -> (max_index, max_val)
        for i, tup in enumerate(tiles):
            chrom, length = tup
            if chrom not in max_indexes or length > max_indexes[chrom][1]:
                max_indexes[chrom] = (i, length)
        asm_order = [tile_tup[0] for i, tile_tup in enumerate(tiles) if max_indexes[tile_tup[0]][0] == i]
        asm_orders[asm] = {chrom: i for i, chrom in enumerate(asm_order)}


    stored_lines = {} # asm -> stored_lines
    with open(args.lengths, 'r', encoding='utf-8') as fin:
        for line in fin:
            asm_name, chrom, length = line.strip().split("\t")
            if asm_name == "bin_id" or asm_name == args.target:
                print(asm_name, chrom, length, sep="\t")
            else:
                if asm_name not in stored_lines:
                    stored_lines[asm_name] = []
                stored_lines[asm_name].append((asm_name, chrom, length))
                if chrom not in asm_orders[asm_name]:
                    asm_orders[asm_name][chrom] = len(asm_orders[asm_name])

    for asm, lines_list in stored_lines.items():
        asm_orders_asm = asm_orders[asm]
        for line in sorted(lines_list, key=lambda x: asm_orders_asm[x[1]]):
            print(*line, sep="\t")


if __name__ == "__main__":
    main()