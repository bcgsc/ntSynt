#!/usr/bin/env python3
'''
Wrapper script for running miniprot-based de novo assessment of synteny blocks
'''
import argparse
import filecmp
import os
import subprocess
import shlex
import sys


def main():
    "Run the busco assessment snakemake file"
    parser = argparse.ArgumentParser(
        description="Run the miniprot-based de novo synteny block assessment using busco lineage proteins")
    parser.add_argument("--assemblies", help="Input assembly files (matching file names in synteny blocks TSV)",
                        required=True, nargs="+")
    parser.add_argument("--blocks", help="ntSynt-formatted synteny blocks", required=True, type=str)
    parser.add_argument("--lineage", help="BUSCO lineage name", required=True, type=str)
    parser.add_argument("--db", help="Path to BUSCO lineage faa file", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use per process [8]", default=8, type=int)
    parser.add_argument("--processes", help="Number of parallel processes to allow at a time [5]", default=5, type=int)
    parser.add_argument("-p", "--prefix", help="Prefix for output files [ntsynt_busco_assessment]",
                        type=str, default="ntsynt_busco_assessment")
    parser.add_argument("--outs", help="Value for --outs parameter of miniprot [0.9]", default=0.9, type=float)
    parser.add_argument("-n", help="Run dry-run of snakemake file", action="store_true")
    parser.add_argument("-F", help="Force all steps to run", action="store_true")
    parser.add_argument("--verbose", help="Verbose assesment output to stderr", action="store_true")

    args = parser.parse_args()

    assembly_basenames = []

    for assembly in args.assemblies:
        #Check if need to make soft links to bring the assembly files into the current WD
        assembly_basename = os.path.basename(assembly)
        assembly_basenames.append(assembly_basename)
        if assembly_basename == assembly:
            continue
        if os.path.exists(assembly_basename):
            if not filecmp.cmp(assembly_basename, assembly):
                sys.exit(f"ERROR: {assembly_basename} exists in the current working directory, "\
                    f"and is different from {assembly}. Thus, cannot make the soft link")
        else:
            cmd = f"ln -s {assembly}"
            ret_code = subprocess.call(shlex.split(cmd))
            assert ret_code == 0


    smk_command = f"snakemake -s /projects/btl/lcoombe/git/ntSynt/analysis_scripts/busco_assessment.smk "\
                    f"--config blocks={args.blocks} assemblies='{assembly_basenames}' threads={args.threads} "\
                    f"prefix={args.prefix} lineage={args.lineage} db={args.db} outs={args.outs}"

    if args.verbose:
        smk_command += " verbose=True"
    if args.n:
        smk_command += " -n"
    if args.F:
        smk_command += " -F"

    smk_command += f" -p --cores={args.processes}"

    print("Running command: ", smk_command)

    smk_command = shlex.split(smk_command)
    ret_code = subprocess.call(smk_command)
    assert ret_code == 0

if __name__ == "__main__":
    main()
