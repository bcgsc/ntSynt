#!/usr/bin/env python3
'''
Using BUSCO runs to assess the quality of detected synteny blocks
'''
import argparse
from collections import namedtuple, defaultdict
import collections
from functools import total_ordering
import re
import sys
import intervaltree
import itertools

Busco = namedtuple("Busco", ["id", "strand", "ms_score", "base_id"])
BuscoList = namedtuple("BuscoList", ["buscos", "strand", "start_reserve", "end_reserve"])

# Matching same BUSCO, potentially different representative sequences
busco_re = re.compile(r'^(\S+)_\d+$')

def read_buscos(assembly: str, busco_fin: str, busco_index: dict) -> set:
    "Read the BUSCO file, storing all Complete BUSCOs"
    stored_buscos = set()
    with open(busco_fin, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            if line.startswith("#"):
                continue
            line = line.split("\t")
            if line[1] != "Complete":
                continue
            busco_id = line[0]
            chrom, start, end, strand, ms_score = line[2:7]
            start, end = int(start), int(end)
            if assembly not in busco_index or chrom not in busco_index[assembly]:
                busco_index[assembly][chrom] = intervaltree.IntervalTree()
            busco_base = re.search(busco_re, busco_id).group(1) if re.search(busco_re, busco_id) else busco_id
            busco_index[assembly][chrom][start: end] = Busco(busco_id, strand, int(ms_score), busco_base)
            stored_buscos.add(busco_id)
    return stored_buscos


def get_list_of_buscos(busco_index: dict, assembly: str, chrom: str, start: int,
                       end: int, strand: str, common_buscos: set) -> BuscoList:
    "Given a synteny block, find all BUSCOs that map to that region, return as Busco_list namedtuple"
    busco_list = []
    start_reserve = []
    end_reserve = []
    if chrom not in busco_index[assembly]:
        # Check for substring in the other direction - BUSCO can sometimes shorten the sequence names
        possible_match = []
        for busco_chrom in busco_index[assembly].keys():
            if busco_chrom in chrom:
                possible_match.append(busco_chrom)
        if len(possible_match) == 1:
            chrom = possible_match.pop()
        elif len(possible_match) != 1:
            sys.exit(f"No match for chrom {chrom} in BUSCO TSV")
    for busco in sorted(busco_index[assembly][chrom].envelop(start, end), key=lambda x: (x.begin, x.end, x.data.id)):
        if busco.data.id in common_buscos:
            busco_list.append(busco.data)
    busco_set = set(busco_list)
    for busco in sorted(busco_index[assembly][chrom][start: end], key=lambda x: (x.begin, x.end, x.data.id)):
        if busco.data.id in common_buscos and busco.data.id not in busco_set:
            if busco.begin <= start: # Overlaps the start
                start_reserve.append(busco.data)
            elif busco.begin > start and busco.end > end:
                end_reserve.append(busco.data)
    return BuscoList(busco_list, strand, start_reserve, end_reserve)


def extend_busco_list(buscos: list) -> None:
    "Given the list of BUSCO lists, use the reserves to extend the lists if needed"
    present_buscos = {busco_id for busco_asm_list in buscos for busco_id in busco_asm_list.buscos}
    for busco_asm_list in buscos:
        for busco in reversed(busco_asm_list.start_reserve):
            if busco in present_buscos:
                busco_asm_list.buscos.insert(0, busco)
        for busco in busco_asm_list.end_reserve:
            if busco in present_buscos:
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
                                       busco_asm_list.start_reserve, busco_asm_list.end_reserve))
    return return_buscos


def filter_repeat(buscos_list: list) -> list:
    "Filter out BUSCOs that (after merging representative) are not found once per assembly"
    dup_buscos = set()
    for busco_asm_list in buscos_list:
        asm_buscos = set()
        for busco in busco_asm_list.buscos:
            if busco.base_id in asm_buscos:
                dup_buscos.add(busco.base_id)
            else:
                asm_buscos.add(busco.base_id)

    new_list = []
    for busco_asm_list in buscos_list:
        new_list.append(BuscoList(buscos=[busco for busco in busco_asm_list.buscos if busco.base_id not in dup_buscos],
                                  strand=busco_asm_list.strand, start_reserve=busco_asm_list.start_reserve,
                                  end_reserve=busco_asm_list.end_reserve))
    return new_list


def check_buscos(buscos_list: list, block_id: str, error_counter: dict) -> None:
    "Given a list of buscos per synteny block, assess if they are consistent"
    i = 0
    # Try to extend the lists if the start/ends are not the same/compatible
    extend_busco_list(buscos_list)
    buscos_list = merge_representative(buscos_list)
    buscos_list = filter_repeat(buscos_list)

    total_buscos = len(set.intersection(*[{busco_id_strand.base_id for busco_id_strand in busco.buscos} \
        for busco in buscos_list]))
    error_counter["total_buscos_assessed"] += total_buscos

    if len({len(busco.buscos) for busco in buscos_list}) != 1:
        error_counter["different_buscos"].add(block_id)
        return
    num_buscos = len(buscos_list[0].buscos)
    if len(set.union(*[{busco_id_strand.base_id for busco_id_strand in busco.buscos} \
            for busco in buscos_list])) != num_buscos:
        error_counter["different_buscos"].add(block_id)

        return
    while i < len(buscos_list[0].buscos):
        curr_buscos = []
        for assembly_buscos in buscos_list:
            if assembly_buscos.strand == "-":
                curr_buscos.append((assembly_buscos.buscos[num_buscos - 1 - i], assembly_buscos.strand))
            else:
                curr_buscos.append((assembly_buscos.buscos[i], assembly_buscos.strand))
        for j, busco_tup in enumerate(zip(curr_buscos, curr_buscos[1:])):
            busco1, busco2 = busco_tup
            if busco1[0].base_id != busco2[0].base_id:
                if (busco1[1] == "+" and buscos_list[j].buscos[num_buscos - 1 - i].base_id == busco2[0].base_id) or \
                    (busco1[1] == "-" and buscos_list[j].buscos[i].base_id == busco2[0].base_id) :
                    busco1 = (buscos_list[j].buscos[num_buscos - 1 - i], buscos_list[j].strand) \
                        if busco1[1] == "+" else (buscos_list[j].buscos[i], buscos_list[j].strand)
            if busco1[0].base_id != busco2[0].base_id:
                error_counter["wrong_order"].add(block_id)
            elif busco1[0].strand != busco2[0].strand and busco1[1] == busco2[1]:
                error_counter["wrong_ori"].add(block_id)
            elif busco1[0].strand == busco2[0].strand and busco1[1] != busco2[1]:
                error_counter["wrong_ori"].add(block_id)
        i += 1


def assess_synteny_blocks(synteny_blocks: str, common_buscos: set, busco_idx: dict) -> tuple:
    "Assess the synteny blocks based on the BUSCO output"
    curr_block_id = None
    total_blocks = 0
    curr_block_buscos = [] # List Busco_list
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
                    check_buscos(curr_block_buscos, curr_block_id, error_counter)
                    total_blocks += 1
                curr_block_id = block_id
                curr_block_buscos = [get_list_of_buscos(busco_idx, assembly, chrom, start, end, strand,
                                                            common_buscos)]
    if curr_block_id is not None:
        check_buscos(curr_block_buscos, curr_block_id, error_counter)
        total_blocks += 1

    return error_counter, total_blocks


def main():
    "Run BUSCO-based synteny block assessment"
    parser = argparse.ArgumentParser(description="BUSCO-based synteny block assessment")
    parser.add_argument("--blocks", help="Synteny blocks TSV", required=True)
    parser.add_argument("--busco", help="TSV file with two columns: assembly filename, full path to BUSCO full TSV",
                        required=True)
    parser.add_argument("-p", "--prefix", help="Prefix of output files", default="ntsynt_busco_assessment",
                        type=str)
    args = parser.parse_args()

    busco_index = defaultdict(dict) # assembly -> chromosome -> IntervalTree (data: BUSCO_ID, strand)
    common_buscos = None
    with open(args.busco, 'r', encoding="utf-8") as fin:
        for line in fin:
            assembly, busco_fin = line.strip().split("\t")
            present_buscos = read_buscos(assembly, busco_fin, busco_index)
            if common_buscos is None:
                common_buscos = present_buscos
            else:
                common_buscos = common_buscos.intersection(present_buscos)

    error_counter, total_blocks = assess_synteny_blocks(args.blocks, common_buscos, busco_index)

    with open(args.prefix + ".summary.tsv", 'w', encoding="utf-8") as fout:
        fout.write('Total_assessed_buscos\tTotal_blocks\tCorrect_blocks\twrong_orientation\twrong_order\tdifferent_buscos\n')
        correct_blocks = total_blocks - sum({len(blocks) for stat, blocks in error_counter.items() if stat != "total_buscos_assessed"})
        outstr = f"{error_counter['total_buscos_assessed']}\t{total_blocks}\t{correct_blocks}\t{len(error_counter['wrong_ori'])}\t"\
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
