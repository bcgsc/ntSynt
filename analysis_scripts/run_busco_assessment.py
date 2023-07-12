#!/usr/bin/env python3
'''
Wrapper script for running miniprot-based de novo assessment of synteny blocks
'''
import argparse
import os
import subprocess
import shlex


def main():
    "Run the busco assessment snakemake file"
    parser = argparse.ArgumentParser(
        description="Run the miniprot-based de novo synteny block assessment using busco lineage faas")
    parser.add_argument("--assemblies", help="Input assembly files (matching file names in synteny blocks TSV)",
                        required=True, nargs="+")
    parser.add_argument("--blocks", help="ntSynt-formatted synteny blocks", required=True, type=str)
    parser.add_argument("--lineage", help="BUSCO lineage name", required=True, type=str)
    parser.add_argument("--db", help="BUSCO lineage faa file to use", required=True)
    parser.add_argument("-t", "--threads", help="Number of threads to use per process [8]", default=8, type=int)
    parser.add_argument("--processes", help="Number of parallel processes to allow at a time [5]", default=5, type=int)
    parser.add_argument("-p", "--prefix", help="Prefix for output files [ntsynt_busco_assessment]",
                        type=str, default="ntsynt_busco_assessment")
    parser.add_argument("--outdir", help="Name of output directory for results [ntsynt_busco_assessment_dir]",
                        default="ntsynt_busco_assessment_dir")
    parser.add_argument("-n", help="Run dry-run of snakemake file", action="store_true")
    parser.add_argument("-F", help="Force all steps to run", action="store_true")

    args = parser.parse_args()

    assembly_basenames = [os.path.basename(assembly) for assembly in args.assemblies]
    blocks_fullpath = os.path.abspath(args.blocks)

    smk_command = f"snakemake -s /projects/btl/lcoombe/git/ntSynt/analysis_scripts/busco_assessment.smk "\
                    f"--config blocks={blocks_fullpath} assemblies='{assembly_basenames}' threads={args.threads} "\
                    f"prefix={args.prefix} lineage={args.lineage} db={args.db} -p --cores={args.processes}"
    if args.n:
        smk_command += " -n"
    if args.F:
        smk_command += " -F"

    print("Running command: ", smk_command)

    smk_command = shlex.split(smk_command)
    ret_code = subprocess.call(smk_command)
    assert ret_code == 0

if __name__ == "__main__":
    main()
