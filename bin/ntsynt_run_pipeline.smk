#!/usr/bin/env snakemake -s
import shutil
import sys
import os

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

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


# Get path to the base directory (where snakemake file is located)
script_path = workflow.basedir

rule all:
    input: expand("{ref}.k{k}.w{w}.tsv", ref=references, k=k, w=w),
            expand("{prefix}.merge_collinear.synteny_blocks.tsv", prefix=prefix)

rule faidx:
    input: fa="{file}"
    output: "{file}.fai"
    shell: "samtools faidx {input.fa}"

rule make_solid_bf:
    input: refs=references
    output: expand("{prefix}.solid.bf", prefix=prefix)
    params: options=expand("-p {prefix}.solid --fpr {fpr} -k {k} -t {threads}", prefix=prefix, fpr=fpr, k=k, threads=threads),
            path_to_script=expand("{base_path}/ntsynt_make_solid_bf", base_path=script_path)
    shell: "{params.path_to_script} --genome {input.refs} {params.options}"

# For posterity, included but this is experimental
rule make_repeat_bf:
    input: refs=references
    output: expand("{prefix}.repeat.bf", prefix=prefix)
    params: options=expand("-p {prefix}.repeat --fpr {fpr} -k {k} -t {threads}", prefix=prefix, fpr=fpr, k=k, threads=threads),
            path_to_script=expand("{base_path}/ntsynt_make_repeat_bfs.py", base_path=script_path)
    shell: "{params.path_to_script} --genome {input.refs} {params.options}"

rule indexlr:
    input: fa="{fasta}",
            solid=expand("{prefix}.solid.bf", prefix=prefix) if solid is True else [],
            repeat=expand("{prefix}.repeat.bf", prefix=prefix) if repeat is True else []
    output: expand("{{fasta}}.k{k}.w{w}.tsv", k=k, w=w)
    params: options=expand("-k {k} -w {w} -t {t} --long --seq --pos", k=k, w=w, t=threads),
            bf="-s" if solid is True else [],
            repeat="-r" if repeat is True else []
    shell: "indexlr {params.options} {params.bf} {input.solid} {params.repeat} {input.repeat} {input.fa} > {output}"

rule ntsynt_synteny:
    input: mx=expand("{fasta}.k{k}.w{w}.tsv", fasta=references, k=k, w=w),
            solid=expand("{prefix}.solid.bf", prefix=prefix) if solid is True else [],
            repeat=expand("{prefix}.repeat.bf", prefix=prefix) if repeat is True else [],
            fais=expand("{fasta}.fai", fasta=references)
    output: expand("{prefix}.merge_collinear.synteny_blocks.tsv", prefix=prefix)
    params: path_to_script=expand("{base_path}/ntsynt_run.py", base_path=script_path),
            options=expand("-k {k} -w {w} --w-rounds {w_rounds} -p {prefix} --bp {indel_merge} --collinear-merge {collinear_merge} --simplify-graph",
                            k=k, w=w, w_rounds=[w_rounds], prefix=prefix, indel_merge=indel_merge, collinear_merge=collinear_merge),
            solid_bf="--solid" if solid is True else [],
            repeat_bf="--repeat" if repeat is True else []        
    shell: "python3 {params.path_to_script} celegansN2_referencegenome.fna.k24.w1000.tsv survivor_F.fa.k24.w1000.tsv {params.options} {params.solid_bf} {input.solid} \
             {params.repeat_bf} {input.repeat}" 