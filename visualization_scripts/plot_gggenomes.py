#!/usr/bin/env python3
'''
Driver script for running distance estimation + gggenomes snakemake pipeline
'''
import argparse
import shlex
import subprocess
import os

def read_fai_files(fai_file):
    "Read the FAI files from the file of files"
    fai_list = []
    with open(fai_file, 'r', encoding="utf-8") as fin:
        for line in fin:
            fai_list.append(line.strip())
    return fai_list

def main():
    "Run distance estimation and gggenomes for ntSynt results"
    parser = argparse.ArgumentParser(
        description="Running ntSynt synteny block distance estimation and generating a ribbon plot")
    parser.add_argument("--blocks", help="ntSynt synteny blocks TSV", required=True, type=str)
    parser.add_argument("--name_conversion", help="TSV for converting names in the blocks TSV (old -> new). NOTE: new names cannot have spaces.",
                        required=False, type=str)
    parser.add_argument("--fais", 
                        help="FAI files for all input assemblies. Can be a list of a file with one FAI path per line.",
                        nargs="+", required=True, type=str)
    parser.add_argument("--prefix", help="Prefix for output files [ntSynt_distance-est]", required=False, type=str,
                        default="ntSynt_distance-est")
    parser.add_argument("--scale", help="Length of scale bar in bases (default 1Gbp)", required=False, type=int,
                        default=1e9)
    parser.add_argument(
        "--ribbon_adjust",
        help="Ratio for adjusting spacing between ribbon plot and cladogram. "
        "Increase if ribbon plot labels are cut off, and decrease to reduce  "
        "the white space between the ribbon plot and cladogram",
        default=0.1, type=float, required=False
    )
    parser.add_argument(
        "--cladogram_adjust",
        help="Ratio for adjusting xlimit of the cladogram. "
        "Increase if cladogram labels are cut off, and decrease to reduce  "
        "the white space on the right side of the cladogram",
        default=0.1, type=float, required=False
    )
    parser.add_argument("-f", "--force", help="Force a re-run of the script", action="store_true")
    parser.add_argument("-n", help="Dry-run for snakemake pipeline", action="store_true")
    
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.realpath(__file__))

    if len(args.fais) == 1:
        args.fais = read_fai_files(args.fais[0])

    cmd = f"snakemake -s {base_dir}/plot_gggenomes.smk --cores 2 " \
           f"-p --cores 2 " \
           f"--config " \
           f"prefix={args.prefix} " \
           f"blocks={args.blocks} " \
           f"name_conversion={args.name_conversion} " \
           f"fai='{args.fais}' " \
           f"ribbon_adjust={args.ribbon_adjust} " \
           f"cladogram_adjust={args.cladogram_adjust} " \
           f"scale={args.scale} "
    if args.force:
        cmd += " -F "
    if args.n:
        cmd += " -n "
    print(cmd)

    subprocess.check_call(shlex.split(cmd))

if __name__ == "__main__":
    main()
