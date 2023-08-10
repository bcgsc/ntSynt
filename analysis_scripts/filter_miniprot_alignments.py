#!/usr/bin/env python3
'''
Given miniprot mappings in PAF format, filter to retain unique mappings based on a coverage threshold
'''
import argparse
from collections import namedtuple
import re

PafEntry = namedtuple("PafEntry", ["busco", "query_len", "query_start", "query_end", "ref", "ref_start",
                                   "ref_end", "strand", "ms_score"])

ms_re = re.compile(r'ms:i:(\d+)') # ms score - DP alignment score excluding introns

def filter_paf(paf_file: str, cov_threshold: float) -> None:
    "Filter the given PAF file using the coverage threshold for the query"
    current_query = None
    print_alignment = True
    with open(paf_file, 'r', encoding="utf-8") as fin:
        for line in fin:
            ms_score = int(re.search(ms_re, line).group(1))
            line = line.strip().split("\t")
            busco, query_len, query_start, query_end, strand, ref, _, ref_start, ref_end = line[:9]
            query_len, query_start, query_end = map(int, [query_len, query_start, query_end])
            if current_query is not None and busco != current_query.busco and print_alignment: # Ensure mapping unique
                # Print the current alignment if over coverage threshold
                query_cov = (current_query.query_end - current_query.query_start)/current_query.query_len
                if query_cov >= cov_threshold:
                    print(current_query.busco, current_query.ref, current_query.ref_start,
                        current_query.ref_end, current_query.strand, current_query.ms_score, sep="\t")
            elif current_query is not None and busco != current_query.busco:
                print_alignment = True # Refresh the print status - print_alignment must have been false before
            elif current_query is not None and busco == current_query.busco:
                print_alignment = False # This busco has multiple mappings
            current_query = PafEntry(busco, query_len, query_start, query_end,
                                     ref, ref_start, ref_end, strand, ms_score)

    if current_query is not None and print_alignment: # Catch last line
        # Print the current alignment
        query_cov = (current_query.query_end - current_query.query_start)/current_query.query_len
        if query_cov >= cov_threshold:
            print(current_query.busco, current_query.ref, current_query.ref_start,
                    current_query.ref_end, current_query.strand, current_query.ms_score, sep="\t")

def main():
    "Filter the input PAF"
    parser = argparse.ArgumentParser(description="Filter miniprot PAF mappings, and output in BUSCO-like TSV format")
    parser.add_argument("PAF", help="miniprot mappings in PAF format")
    parser.add_argument("--cov", help="Required query coverage [0.95]", default=0.95, type=float,
                        required=False)

    args = parser.parse_args()

    filter_paf(args.PAF, args.cov)

if __name__ == "__main__":
    main()
