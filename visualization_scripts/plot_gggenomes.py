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

def check_name_conversion(name_conversion, parser):
    "Check if the name conversion file is valid"
    with open(name_conversion, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip().split("\t")
            _, new = line
            if " " in new:
                raise parser.error(
                    f"New name {new} cannot have spaces. "
                    f"Use the special character '_' instead. "
                    f"All underscores in the new name will be "
                    f"converted to spaces for the final output.")

def main():
    "Run distance estimation and gggenomes for ntSynt results"
    parser = argparse.ArgumentParser(
        description="Generate a ribbon plot to visualize ntSynt synteny blocks")
    required_group = parser.add_argument_group('required arguments')
    main_formatting_group = parser.add_argument_group('main plot formatting arguments')
    block_filter_group = parser.add_argument_group('block filtering arguments')
    output_group = parser.add_argument_group('output arguments')
    execution_group = parser.add_argument_group('execution arguments')

    required_group.add_argument("--blocks", help="ntSynt synteny blocks TSV", required=True, type=str)
    required_group.add_argument("--fais",
                        help="FAI files for all input assemblies. Can be a list or a file with one FAI path per line.",
                        nargs="+", required=True, type=str)
    main_formatting_group.add_argument("--name_conversion",
                        help="TSV for converting names in the blocks TSV (old -> new). "
                        "IMPORTANT: new names cannot have spaces. If you want to have spaces in the final ribbon plot, "
                        "use the special character '_'. All underscores in the new name will be converted to spaces.",
                        required=False, type=str)
    main_formatting_group.add_argument("--tree", help="User-input tree file in newick format. "
                        "If specified, this tree will be plotted next to the output ribbon plot, "
                        "and used for ordering the assemblies. "
                        "The names in the newick file must match the new names if --name_conversion is specified, "
                        "or the genome file names in the synteny blocks input file otherwise. "
                        "If not specified, the synteny blocks will be used to estimate pairwise distances "
                        "for the assembly ordering and associated tree.",
                        required=False, type=str)
    main_formatting_group.add_argument("--target-genome", help="Target genome. If specified, "
                        "this genome will be at the top of the ribbon plot, "
                        "with ribbons coloured based on its chromosomes and "
                        "(if applicable) other chromosomes normalized to it. "
                        "If not specified, the top genome will be arbitrary.", required=False, type=str)
    main_formatting_group.add_argument("--normalize", help="Normalize strand of chromosomes relative to the "
                        "target (top) genome in the ribbon plots",
                        action="store_true")
    block_filter_group.add_argument("--indel", help="Indel size threshold [50000]", default=50000, type=int)
    block_filter_group.add_argument("--length", help="Minimum synteny block length [100000]", default=100000, type=int)
    block_filter_group.add_argument("--seq_length", help="Minimum sequence length [500000]", default=500000, type=float)
    block_filter_group.add_argument("--keep", help="List of assembly_name:chromosome to show in visualization. "
                                    "All chromosomes with links to the specified chromosomes will also be shown.",
                                    nargs="+", required=False, type=str)
    main_formatting_group.add_argument("--centromeres",
                        help="TSV file with centromere positions. Must have the headers: bin_id,seq_id,start,end. "\
                            "bin_id must match the new names from --name_conversion or "
                            "the assembly names if --name_conversion is not specified. "\
                            "seq_id is the chromosome name.", required=False, type=str)
    main_formatting_group.add_argument("--haplotypes", help="File listing haplotype assembly names: TSV, "
                        "maternal/paternal assemblies separated by tabs.",
                        required=False, type=str)
    output_group.add_argument("--prefix", help="Prefix for output files [ntSynt_ribbon-plot]", required=False, type=str,
                        default="ntSynt_ribbon-plot")
    output_group.add_argument("--format", help="Output format of ribbon plot [png]",
                        required=False, choices=["png", "pdf"], default="png")
    output_group.add_argument("--scale", help="Length of scale bar in bases [1e9]", required=False, type=float,
                        default=1e9)
    output_group.add_argument("--height", help="Height of plot in cm [20]", required=False, type=int, default=20)
    output_group.add_argument("--width", help="Width of plot in cm [50]", required=False, type=int, default=50)
    main_formatting_group.add_argument("--no-arrow", help="Only used with --normalize; "
                        "do not draw arrows indicating reverse-complementation",
                        action="store_true")
    output_group.add_argument(
        "--ribbon_adjust",
        help="Ratio for adjusting spacing beside ribbon plot. "
        "Increase if ribbon plot labels are cut off, and decrease to reduce  "
        "the white space to the left of the ribbon plot [0.1]",
        default=0.1, type=float, required=False
    )
    execution_group.add_argument("-f", "--force", help="Force a re-run of the entire pipeline", action="store_true")
    execution_group.add_argument("-n", help="Dry-run for snakemake pipeline", action="store_true")

    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.realpath(__file__))

    if args.name_conversion:
        check_name_conversion(args.name_conversion, parser)

    if len(args.fais) == 1:
        args.fais = read_fai_files(args.fais[0])

    cmd = f"snakemake -s {base_dir}/plot_gggenomes.smk " \
            f"--cores 2 " \
            f"--config " \
            f"prefix={args.prefix} " \
            f"blocks={args.blocks} " \
            f"name_conversion={args.name_conversion} " \
            f"fai='{args.fais}' " \
            f"ribbon_ratio={args.ribbon_adjust} " \
            f"scale={args.scale} " \
            f"indel_threshold={args.indel} " \
            f"min_length={args.length} " \
            f"format={args.format} " \
            f"height={args.height} " \
            f"width={args.width} " \
            f"min_seq_length={args.seq_length} "

    if args.centromeres:
        cmd += f"centromeres={args.centromeres} "
    if args.normalize:
        cmd += "normalize=True "
    if args.no_arrow:
        cmd += "no_arrow=True "
    if args.target_genome:
        cmd += f"target_genome={args.target_genome} "
    if args.haplotypes:
        cmd += f"haplotypes={args.haplotypes} "
    if args.keep:
        cmd += f"keep=\"{args.keep}\" "
    if args.tree:
        cmd += f"tree={args.tree} "
        target = "gggenomes_ribbon_plot_tree"
    else:
        target = "all"
    if args.force:
        cmd += " -F "
    if args.n:
        cmd += " -n "
    cmd += f"-p {target}"
    print(cmd)

    subprocess.check_call(shlex.split(cmd))

if __name__ == "__main__":
    main()
