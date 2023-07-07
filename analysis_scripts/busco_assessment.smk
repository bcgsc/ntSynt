#!/usr/bin/env snakemake -s
import shutil
import sys
import re
import os

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

# Read in parameters
assemblies = config["assemblies"] if "assemblies" in config else None # Full paths to assemblies
assemblies_base = [os.path.basename(assembly) for assembly in assemblies]
synteny_blocks = config["blocks"] if "blocks" in config else None
threads = config["t"] if "t" in config else 8
lineage = config["lineage"] if "lineage" in config else None
lineage_name = os.path.basename(lineage) if os.path.isdir(lineage) else lineage
prefix = config["prefix"] if "prefix" in config else "ntsynt_busco_assessment"

rule all:
    input: expand("{prefix}.summary.tsv", prefix=prefix)

rule busco:
    input: "{assembly}"
    output: expand("busco_{{assembly}}/run_{lineage}/full_table.tsv", lineage=lineage_name)
    params: options=expand("-l {lineage} -c {threads} -m genome -f --tar", lineage=lineage, threads=threads)
    shell:
        "busco -i {input} -o busco_{wildcards.assembly} {params.options} &> busco_{wildcards.assembly}.busco.log"

rule make_config:
    input: expand("busco_{assembly}/run_{lineage}/full_table.tsv", assembly=assemblies_base, lineage=lineage_name)
    output: out_file=expand("{prefix}_config.tsv", prefix=prefix)
    run:
        print(output.out_file[0])
        with open(output.out_file[0], 'w') as fout:
            for busco_file in input:
                if assembly_match := re.search(r'^busco_(\S+)\/run', busco_file):
                    fout.write(f"{assembly_match.group(1)}\t{busco_file}\n")

rule run_assessment:
    input: config=expand("{prefix}_config.tsv", prefix=prefix), blocks=synteny_blocks
    output: expand("{prefix}.summary.tsv", prefix=prefix)
    params: options=expand("--prefix {prefix}", prefix=prefix)
    shell:
        "python3 /projects/btl/lcoombe/git/ntSynt/analysis_scripts/busco_assessment.py  --blocks {input.blocks} --busco {input.config} {params.options}"
