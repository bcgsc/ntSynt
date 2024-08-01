#!/usr/bin/env python3
'''
Rename the second columns of the ntSynt synteny blocks (assembly name)
Useful when naming synteny blocks based on some other name is clearer than assembly names
'''
import sys
from format_blocks_gggenomes import read_name_conversions

def rename_synteny_tsv(tsv_in, names_dict):
    "Rename the first columns of ntSynt synteny blocks"
    with open(tsv_in, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[1] in names_dict:
                line[1] = names_dict[line[1]]
            print("\t".join(line))


def main():
    "Rename the second columns (assembly name) of the ntSynt synteny blocks"
    if len(sys.argv[1:]) != 2:
        print(f"Usage: {sys.argv[0]} <synteny TSV file> <old name -> new name TSV>")
        sys.exit()

    synteny_tsv = sys.argv[1]
    correspondence_tsv = sys.argv[2]

    names_dict = read_name_conversions(correspondence_tsv)
    rename_synteny_tsv(synteny_tsv, names_dict)

if __name__ == "__main__":
    main()
