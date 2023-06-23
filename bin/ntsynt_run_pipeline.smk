#!/usr/bin/env snakemake -s
import shutil
import sys
import os

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

# Read in parameters
references = config["references"] if "references" in config else "Must specify 'references'"
k = config["k"] if "k" in config else 24
w = config["w"] if "w" in config else 1000
fpr = config["fpr"] if "fpr" in config else 0.025
threads = config["t"] if "t" in config else 4
prefix = config["prefix"] if "prefix" in config else "ntSynt_out"
solid = config["solid"] if "solid" in config else True
repeat = config["repeat"] if "repeat" in config else False
w_rounds = config["w_rounds"] if "w_rounds" in config else [100, 10, 5]
indel_merge = config["indel_merge"] if "indel_merge" in config else 500
collinear_merge = config["collinear_merge"] if "collinear_merge" in config else "3w"
simplify_graph = config["simplify_graph"] if "simplify_graph" in config else True
benchmark = config["benchmark"] if "benchmark" in config else False

# If want to benchmark, use memusg or /usr/bin/time
# The snakemake benchmarking is quite inaccurate for shorter runtimes
benchmark_path = "" 
if benchmark:
    if benchmark_memusg := shutil.which("memusg"):
        benchmark_path = f"{benchmark_memusg} -t"
    elif benchmark_time := shutil.which("time"):
        benchmark_path = f"{benchmark_time} -v"
    else:
        print("WARNING: memusg and time not found in PATH. Benchmarks will not be tallied.")


# Get path to the base directory (where snakemake file is located)
script_path = workflow.basedir

rule all:
    input: expand("{ref}.k{k}.w{w}.tsv", ref=references, k=k, w=w),
            expand("{prefix}.synteny_blocks.tsv", prefix=prefix)

rule faidx:
    input: fa="{file}"
    output: "{file}.fai"
    params: benchmarking=expand("{benchmark_path} -o {{file}}.faidx.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} samtools faidx {input.fa}"

rule make_solid_bf:
    input: refs=references
    output: expand("{prefix}.solid.bf", prefix=prefix)
    params: options=expand("-p {prefix}.solid --fpr {fpr} -k {k} -t {threads}", prefix=prefix, fpr=fpr, k=k, threads=threads),
            path_to_script=expand("{base_path}/ntsynt_make_solid_bf", base_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.make_solid_bf.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} {params.path_to_script} --genome {input.refs} {params.options}"

# For posterity, included but this is experimental
rule make_repeat_bf:
    input: refs=references
    output: expand("{prefix}.repeat.bf", prefix=prefix)
    params: options=expand("-p {prefix}.repeat --fpr {fpr} -k {k} -t {threads}", prefix=prefix, fpr=fpr, k=k, threads=threads),
            path_to_script=expand("{base_path}/ntsynt_make_repeat_bfs.py", base_path=script_path),
            benchmarking=expand("{benchmark_path} -o {prefix}.make_repeat_bf.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else []
    shell: "{params.benchmarking} {params.path_to_script} --genome {input.refs} {params.options}"

rule indexlr:
    input: fa="{fasta}",
            solid=expand("{prefix}.solid.bf", prefix=prefix) if solid is True else [],
            repeat=expand("{prefix}.repeat.bf", prefix=prefix) if repeat is True else []
    output: expand("{{fasta}}.k{k}.w{w}.tsv", k=k, w=w)
    params: options=expand("-k {k} -w {w} -t {t} --long --seq --pos", k=k, w=w, t=threads),
            bf="-s" if solid is True else [],
            repeat="-r" if repeat is True else [],
            benchmarking=expand("{benchmark_path} -o {{fasta}}.indexlr.time", benchmark_path=benchmark_path) if benchmark else []
    shell: "{params.benchmarking} indexlr {params.options} {params.bf} {input.solid} {params.repeat} {input.repeat} {input.fa} > {output}"

rule ntsynt_synteny:
    input: mx=expand("{fasta}.k{k}.w{w}.tsv", fasta=references, k=k, w=w),
            solid=expand("{prefix}.solid.bf", prefix=prefix) if solid is True else [],
            repeat=expand("{prefix}.repeat.bf", prefix=prefix) if repeat is True else [],
            fais=expand("{fasta}.fai", fasta=references)
    output: expand("{prefix}.synteny_blocks.tsv", prefix=prefix)
    params: path_to_script=expand("{base_path}/ntsynt_run.py", base_path=script_path),
            options=expand("-k {k} -w {w} --w-rounds {w_rounds} -p {prefix} --bp {indel_merge} --collinear-merge {collinear_merge}",
                            k=k, w=w, w_rounds=[w_rounds], prefix=prefix, indel_merge=indel_merge, collinear_merge=collinear_merge),
            solid_bf="--solid" if solid is True else [],
            repeat_bf="--repeat" if repeat is True else [],
            simplify_graph="--simplify-graph" if simplify_graph is True else [],
            benchmarking=expand("{benchmark_path} -o {prefix}.synteny_blocks.time", benchmark_path=benchmark_path, prefix=prefix) if benchmark else [] 
    shell: "{params.benchmarking} python3 {params.path_to_script} {input.mx} {params.options} {params.solid_bf} {input.solid}  {params.simplify_graph} \
             {params.repeat_bf} {input.repeat}"
