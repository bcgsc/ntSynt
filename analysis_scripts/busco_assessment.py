#!/usr/bin/env python3
'''
Using miniprot mappings of BUSCO proteins to assess the quality of detected synteny blocks
'''
import argparse
from collections import namedtuple, defaultdict
import re
import itertools
import sys
import intervaltree

Busco = namedtuple("Busco", ["id", "strand", "ms_score", "base_id"])
BuscoList = namedtuple("BuscoList", ["buscos", "strand", "start_reserve", "end_reserve", "assembly"])
BlockIdBuscoList = namedtuple("BlockIdBuscoList", ["block_id", "busco_list"])

# Matching same BUSCO, potentially different representative sequences
busco_re = re.compile(r'^(\S+)_\d+$')

def read_mappings(assembly: str, busco_fin: str, busco_index: dict) -> set:
    "Read the filtered, parsed mapping file"
    stored_buscos = set()
    with open(busco_fin, 'r', encoding="utf-8") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            busco_id, chrom, start, end, strand, ms_score = line[:6]
            start, end = int(start), int(end)
            if assembly not in busco_index or chrom not in busco_index[assembly]:
                busco_index[assembly][chrom] = intervaltree.IntervalTree()
            busco_base = re.search(busco_re, busco_id).group(1) if re.search(busco_re, busco_id) else busco_id
            busco_index[assembly][chrom][start:end] = Busco(busco_id, strand, int(ms_score), busco_base)
            stored_buscos.add(busco_base)
    return stored_buscos


def get_list_of_buscos(busco_index: dict, assembly: str, chrom: str, start: int,
                       end: int, strand: str, common_buscos: set) -> BuscoList:
    "Given a synteny block, find all BUSCOs that map to that region, return as Busco_list namedtuple"
    buscos = []
    start_reserve = [] # Extend the coordinates if needed - accounting for mappings near block edges
    end_reserve = []
    for busco in sorted(busco_index[assembly][chrom].envelop(start, end), key=lambda x: (x.begin, x.end)):
        if busco.data.base_id in common_buscos:
            buscos.append(busco.data)
    for busco in sorted(busco_index[assembly][chrom][start: end], key=lambda x: (x.begin, x.end)):
        if busco.data.base_id in common_buscos:
            if busco.begin <= start: # Overlaps the start
                start_reserve.append(busco.data)
            elif busco.end > end: # Overlaps the end
                end_reserve.append(busco.data)
    return BuscoList(buscos, strand, start_reserve, end_reserve, assembly)


def extend_busco_list(buscos: list) -> None:
    "Given the list of BUSCO lists, use the reserves to extend the lists if needed"
    present_buscos = {busco_id.base_id for busco_asm_list in buscos for busco_id in busco_asm_list.buscos}
    for busco_asm_list in buscos:
        for busco in reversed(busco_asm_list.start_reserve):
            if busco.base_id in present_buscos:
                busco_asm_list.buscos.insert(0, busco)
        for busco in busco_asm_list.end_reserve:
            if busco.base_id in present_buscos:
                busco_asm_list.buscos.append(busco)


def merge_representative(buscos: list) -> list:
    "When there are adjacent representative proteins of the same BUSCO, merge them"
    return_buscos = []
    for busco_asm_list in buscos:
        out_list = []
        grouped_list = [(base_id, list(group_iterator)) \
            for base_id, group_iterator in itertools.groupby(busco_asm_list.buscos, key=lambda x: x.base_id)]
        for _, group_iterator in grouped_list:
            best_busco = sorted(list(group_iterator), key=lambda x: x.ms_score, reverse=True)[0]
            out_list.append(best_busco)
        return_buscos.append(BuscoList(out_list, busco_asm_list.strand,
                                       busco_asm_list.start_reserve, busco_asm_list.end_reserve,
                                       busco_asm_list.assembly))
    return return_buscos


def filter_repeat(blocks_list: list) -> list:
    "Filter out BUSCOs that (after merging representative) are found more than once in a given assembly"
    dup_buscos = {busco_list.assembly: set() for busco_list in blocks_list[0].busco_list}
    asm_buscos = {busco_list.assembly: set() for busco_list in blocks_list[0].busco_list}
    for _, buscos_list in blocks_list:
        for busco_asm_list in buscos_list:
            for busco in busco_asm_list.buscos:
                if busco.base_id in asm_buscos[busco_asm_list.assembly]:
                    dup_buscos[busco_asm_list.assembly].add(busco.base_id)
                else:
                    asm_buscos[busco_asm_list.assembly].add(busco.base_id)

    new_list = []
    dup_buscos = {busco for _, busco_list in dup_buscos.items() for busco in busco_list}
    for busco_id, buscos_list in blocks_list:
        new_entry = []
        for busco_asm_list in buscos_list:
            new_entry.append(BuscoList(buscos=[busco for busco in busco_asm_list.buscos if \
                                    busco.base_id not in dup_buscos],
                                    strand=busco_asm_list.strand, start_reserve=busco_asm_list.start_reserve,
                                    end_reserve=busco_asm_list.end_reserve, assembly=busco_asm_list.assembly))
        new_list.append(BlockIdBuscoList(busco_id, new_entry))
    return new_list


def check_buscos(buscos_list: list) -> list:
    "Given a list of buscos per synteny block, assess if they are consistent"
    # Try to extend the lists if the start/ends are not the same/compatible
    extend_busco_list(buscos_list)
    buscos_list = merge_representative(buscos_list)
    return buscos_list

def print_different_buscos(buscos_list: BuscoList, block_id: int) -> None:
    "Verbose printing of BUSCOs that differed between within block in assemblies"
    common_buscos = set.intersection(*[{busco_id_strand.base_id for busco_id_strand in busco.buscos} \
        for busco in buscos_list])
    for busco in buscos_list:
        for my_busco in busco.buscos:
            if my_busco.base_id not in common_buscos:
                print(f"{block_id}\tdifferent_busco\t{my_busco.base_id}\t{busco.assembly}", file=sys.stderr)

def print_wrong_order(id_list1: list, id_list2: list, block_id: int) -> None:
    "Verbose printing of block lists found in different orders"
    print(f"{block_id}\twrong_order\t{id_list1}\t{id_list2}")

def print_wrong_ori(busco1: Busco, busco2: Busco, block_id: int) -> None:
    "Verbose printing of BUSCOs from blocks with wrong orientation"
    print(f"{block_id}\twrong_ori\t{busco1[0].id}{busco1[0].strand}\t{busco1[1]}\t"
          f"{busco2[0].id}{busco2[0].strand}\t{busco2[1]}")

def assess_buscos(buscos_list: list, block_id: str, error_counter: dict, verbose: bool) -> None:
    "Assess the BUSCO lists after filtering, etc"
    total_buscos = len(set.intersection(*[{busco_id_strand.base_id for busco_id_strand in busco.buscos} \
        for busco in buscos_list]))
    error_counter["total_buscos_assessed"] += total_buscos

    # Check for synteny blocks where the BUSCO IDs are not the same
    if len({len(busco.buscos) for busco in buscos_list}) != 1:
        error_counter["different_buscos"].add(block_id)
        if verbose:
            print_different_buscos(buscos_list, block_id)
        return

    num_buscos = len(buscos_list[0].buscos)
    if len(set.union(*[{busco_id_strand.base_id for busco_id_strand in busco.buscos} \
            for busco in buscos_list])) != num_buscos:
        error_counter["different_buscos"].add(block_id)
        if verbose:
            print_different_buscos(buscos_list, block_id)
        return

    # Check for synteny blocks where the order of BUSCOs is not the same
    for j in range(len(buscos_list) - 1):
        base_ids_1 = [busco_id_strand.base_id for busco_id_strand in buscos_list[j].buscos]
        base_ids_2 = [busco_id_strand.base_id for busco_id_strand in buscos_list[j+1].buscos]
        if not (base_ids_1 == base_ids_2 or base_ids_1 == list(reversed(base_ids_2))):
            # they are truly different orders, not just strand
            error_counter["wrong_order"].add(block_id)
            if verbose:
                print_wrong_order(base_ids_1, base_ids_2, block_id)
            return

    # Check for synteny blocks with orientation issues
    i = 0
    while i < len(buscos_list[0].buscos):
        curr_buscos = []
        for assembly_buscos in buscos_list:
            if assembly_buscos.strand == "-":
                curr_buscos.append((assembly_buscos.buscos[num_buscos - 1 - i], assembly_buscos.strand))
            else:
                curr_buscos.append((assembly_buscos.buscos[i], assembly_buscos.strand))
        for busco1, busco2 in zip(curr_buscos, curr_buscos[1:]):
            if (busco1[0].base_id != busco2[0].base_id) or \
                (busco1[0].strand != busco2[0].strand and busco1[1] == busco2[1]) or \
                (busco1[0].strand == busco2[0].strand and busco1[1] != busco2[1]):
                error_counter["wrong_ori"].add(block_id)
                if verbose:
                    print_wrong_ori(busco1, busco2, block_id)
        i += 1


def assess_synteny_blocks(synteny_blocks: str, common_buscos: set, busco_idx: dict, verbose: bool) -> tuple:
    "Assess the synteny blocks based on the mappings"
    curr_block_id = None
    total_blocks = 0
    all_blocks = [] # list of BlockIdBuscoList
    curr_block_buscos = [] # list of Busco_list
    error_counter = {"different_buscos": set(), "wrong_order": set(), "wrong_ori": set(),
                     'total_buscos_assessed': 0}
    with open(synteny_blocks, 'r', encoding="utf-8") as fin:
        for line in fin:
            block_id, assembly, chrom, start, end, strand = line.strip().split("\t")[:-1]
            start, end = int(start), int(end)
            if block_id == curr_block_id:
                curr_block_buscos.append(get_list_of_buscos(busco_idx, assembly, chrom, start, end, strand,
                                                            common_buscos))
            else:
                if curr_block_id is not None:
                    curr_block_buscos = check_buscos(curr_block_buscos)
                    all_blocks.append(BlockIdBuscoList(curr_block_id, curr_block_buscos))
                    total_blocks += 1
                curr_block_id = block_id
                curr_block_buscos = [get_list_of_buscos(busco_idx, assembly, chrom, start, end, strand,
                                                            common_buscos)]
    if curr_block_id is not None:
        curr_block_buscos = check_buscos(curr_block_buscos)
        all_blocks.append(BlockIdBuscoList(curr_block_id, curr_block_buscos))
        total_blocks += 1

    all_blocks = filter_repeat(all_blocks)

    for block_id, busco_list in all_blocks:
        assess_buscos(busco_list, block_id, error_counter, verbose)

    return error_counter, total_blocks


def main():
    "Run miniprot mapping of BUSCO proteins-based synteny block assessment"
    parser = argparse.ArgumentParser(description="Synteny block assessment using miniprot mappings of BUSCO proteins")
    parser.add_argument("--blocks", help="Synteny blocks TSV", required=True)
    parser.add_argument("--mappings",
                        help="TSV file with two columns: assembly filename, path to filtered miniprot mappings",
                        required=True)
    parser.add_argument("-p", "--prefix", help="Prefix of output files", default="ntsynt_busco_assessment",
                        type=str)
    parser.add_argument("--verbose", help="Verbose output to standard error", action="store_true")
    args = parser.parse_args()

    busco_index = defaultdict(dict) # assembly -> chromosome -> IntervalTree (data: Busco)
    common_proteins = None
    with open(args.mappings, 'r', encoding="utf-8") as fin:
        for line in fin:
            assembly, mappings_fin = line.strip().split("\t")
            present_buscos = read_mappings(assembly, mappings_fin, busco_index)
            if common_proteins is None:
                common_proteins = present_buscos
            else:
                common_proteins = common_proteins.intersection(present_buscos)

    error_counter, total_blocks = assess_synteny_blocks(args.blocks, common_proteins, busco_index, args.verbose)

    with open(args.prefix + ".summary.tsv", 'w', encoding="utf-8") as fout:
        fout.write('Total_assessed_buscos\tTotal_blocks\tCorrect_blocks\twrong_orientation\t'\
            'wrong_order\tdifferent_buscos\n')
        correct_blocks = total_blocks - sum({len(blocks) \
            for stat, blocks in error_counter.items() if stat != "total_buscos_assessed"})
        outstr = f"{error_counter['total_buscos_assessed']}\t{total_blocks}\t"\
            f"{correct_blocks}\t{len(error_counter['wrong_ori'])}\t"\
            f"{len(error_counter['wrong_order'])}\t{len(error_counter['different_buscos'])}\n"
        fout.write(outstr)

    with open(args.prefix + ".error-blocks.tsv", 'w', encoding="utf-8") as fout:
        fout.write("Error_type\tBlock_number\n")
        for error_type, error_set in error_counter.items():
            if error_type == "total_buscos_assessed":
                continue
            for block in error_set:
                fout.write(f"{error_type}\t{block}\n")

if __name__ == "__main__":
    main()
