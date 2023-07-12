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
prefix = config["prefix"] if "prefix" in config else "ntsynt_busco_assessment"
ref_db = config["db"] if "db" in config else None

rule all:
    input: expand("{prefix}.summary.tsv", prefix=prefix)

rule miniprot:
    input: fa="{assembly}", query=expand("{ref_db}", ref_db=ref_db)
    output: expand("{{assembly}}.{lineage}.paf", lineage=lineage)
    params: options=expand("-t {t}", t=threads)
    log: stderr="logs/{assembly}.miniprot.log"
    shell: "miniprot {params.options} {input.fa} {input.query} > {output} 2> {log.stderr}"

rule filter_miniprot:
    input: rules.miniprot.output
    output: expand("{{assembly}}.{lineage}.paf.tsv", lineage=lineage)
    shell: "python3 /projects/btl/lcoombe/git/ntSynt/analysis_scripts/filter_miniprot_alignments.py {input} > {output}"

rule make_config_miniprot:
    input:  expand("{assembly}.{lineage}.paf.tsv", lineage=lineage, assembly=assemblies)
    output: out_file=expand("{prefix}_config.tsv", prefix=prefix)
    run:
        with open(output.out_file[0], 'w') as fout:
            for busco_file in input:
                if assembly_match := re.search(rf'^(\S+).{lineage}', busco_file):
                    fout.write(f"{assembly_match.group(1)}\t{busco_file}\n")

rule run_assessment:
    input: config=rules.make_config_miniprot.output, blocks=synteny_blocks
    output: expand("{prefix}.summary.tsv", prefix=prefix)
    params: options=expand("--prefix {prefix}", prefix=prefix)
    shell:
        "python3 /projects/btl/lcoombe/git/ntSynt/analysis_scripts/busco_assessment.py  --blocks {input.blocks} --busco {input.config} {params.options}"
