#!/usr/bin/env python3
'''
ntSynt: Genome synteny detection using dynamic undirected minimizer graphs
Written by Lauren Coombe (@lcoombe)
'''
import argparse
import os
import shlex
import subprocess

def main():
    "Run ntSynt snakemake file"
    parser = argparse.ArgumentParser(description="ntSynt: Genome synteny detection using dynamic minimizer graphs")
    parser.add_argument("fastas", help="Input genome fasta files", nargs="+")
    parser.add_argument("-k", help="Minimizer kmer size [24]", type=int, required=False, default=24)
    parser.add_argument("-w", help="Minimizer window size [1000]", type=int, required=False, default=1000)
    parser.add_argument("-t", help="Number of threads [4]", type=int, default=4)
    parser.add_argument("--fpr", help="False positive rate for Bloom filter creation [0.025]",
                        default=0.025, type=float)
    parser.add_argument("--no-solid", help="Do not use the solid BF for minimizer graph creation",
                        action="store_true")
    parser.add_argument("--no-simplify-graph", help="Do not run graph simplification on minimizer graph",
                        action="store_true")
    parser.add_argument("-p", "--prefix", help="Prefix for ntSynt output files [ntSynt.k<k>.w<w>]",
                        required=False)
    parser.add_argument("--merge", help="Multiple of window size used for collinear synteny block merging [3]",
                        default=3, type=int)
    parser.add_argument("--w_rounds", help="List of window sizes for iterative rounds [100 10 5]",
                        nargs="+", default=[100, 10, 5])
    parser.add_argument("--indel", help="Threshold for indel detection [500]", default=500, type=int)
    parser.add_argument("--dry-run", help="Print out the commands that will be executed", action="store_true")
    parser.add_argument("--benchmark", help="Store benchmarks for each step of the ntSynt pipeline",
                        action="store_true")
    parser.add_argument("-v", "--version", action="version", version="ntSynt v0.0.1")

    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.realpath(__file__))
    if not args.prefix:
        args.prefix = f"ntSynt.k{args.k}.w{args.w}"

    args.w_rounds = " ".join(map(str, args.w_rounds))
    command = f"snakemake -s {base_dir}/ntsynt_run_pipeline.smk -p --cores 1 " \
                f"--config references='{args.fastas}' k={args.k} w={args.w} t={args.t} fpr={args.fpr} " \
                f"prefix={args.prefix} w_rounds='{args.w_rounds}' indel_merge={args.indel} " \
                f"collinear_merge={args.merge}w "
    command += "solid=False " if args.no_solid else "solid=True "
    command += "simplify_graph=False " if args.no_simplify_graph else "simplify_graph=True "
    command += "benchmark=True " if args.benchmark else "benchmark=False "

    if args.dry_run:
        command += " -n"

    print(f"Running {command}", flush=True)

    command = shlex.split(command)

    ret = subprocess.call(command)
    if ret != 0:
        raise subprocess.SubprocessError("ntSynt failed - check the logs for the error.")

if __name__ == "__main__":
    main()
