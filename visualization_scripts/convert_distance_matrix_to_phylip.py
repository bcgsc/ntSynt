#!/usr/bin/env python3
'''
Convert a pairwise distance matrix to a phylip format
'''
import argparse
from collections import defaultdict
from gggenomes_normalize_strands import read_name_conversions

def load_distance_matrix(tsv_in):
    "Load the distance matrix into a dictionary"
    distances = defaultdict(dict)
    with open(tsv_in, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            asm1, asm2, distance = line[:3]
            try:
                distance = float(distance)
            except ValueError:
                continue
            distances[asm1][asm2] = distance
            distances[asm2][asm1] = distance
    return distances

def print_phylip_lower_triangle(distances_dict, conversion_dict):
    "Print the distances in phylip lower triangle format"
    ordered_assemblies = sorted(list(distances_dict.keys()))
    print(len(ordered_assemblies))
    for i, asm1 in enumerate(ordered_assemblies):
        dist_str = [conversion_dict[asm1]] if conversion_dict else [asm1]
        for j in range(0, i):
            asm2 = ordered_assemblies[j]
            dist_str.append(distances_dict[asm1][asm2])
        print(" ".join(map(str,dist_str)))

def main():
    "Convert pairwise distances to phylip format matrix (lower triangle)"
    parser = argparse.ArgumentParser(description='Convert a pairwise distance matrix to a phylip format')
    parser.add_argument('distance_tsv', help='TSV file of pairwise distances')
    parser.add_argument("--convert", help="TSV file to convert assembly names as desired. " \
                        "Expects the existing name in first column, desired rename in second column.",
                        required=False, type=str)
    args = parser.parse_args()

    distances = load_distance_matrix(args.distance_tsv)

    conversion_dict = None
    if args.convert:
        conversion_dict = read_name_conversions(args.convert)

    print_phylip_lower_triangle(distances, conversion_dict)


if __name__ == "__main__":
    main()
