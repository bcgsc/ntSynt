#!/usr/bin/env snakemake -s
import shutil
import sys
import os

onsuccess:
    shutil.rmtree(".snakemake", ignore_errors=True)

# Parameters
synteny_blocks = config["blocks"] if "blocks" in config else None
name_conversion = config["name_conversion"] if "name_conversion" in config else None
prefix = config["prefix"] if "prefix" in config else "ntSynt_gggenomes"
indel_threshold = config["indel_threshold"] if "indel_threshold" in config else 50000
min_len = config["min_length"] if "min_length" in config else 50000
fais = config["fai"] if "fai" in config else None
ribbon_ratio = config["ribbon_ratio"] if "ribbon_ratio" in config else 0.1
cladogram_ratio = config["cladogram_adjust"] if "cladogram_adjust" in config else 0.1
scale = config["scale"] if "scale" in config else 1e9
blocks_no_suffix = os.path.basename(synteny_blocks).removesuffix(".tsv")

def sort_fais(fai_list, name_conversion, orders):
    "Based on the name conversion TSV, sort the FAIs based on orders"
    # Read in name conversions
    names_dict = {} # new name -> old name
    with open(name_conversion, 'r') as fin:
        for line in fin:
            old, new = line.strip().split("\t")
            names_dict[new] = old

    # Make order dictionary, converting new names to FAI-compatible
    order_dict = {}
    with open(orders, 'r') as fin:
        i = 0
        for line in fin:
            new_name = line.strip()
            order_dict[names_dict[new_name]] = i
            i += 1

    return " ".join(sorted(fai_list, key=lambda x: order_dict[os.path.basename(x).removesuffix(".fai")]))

def get_first_asm_and_fai(orders, name_conversion, fai_list):
    "Given the orders files, and a list of corresponding FAI, return the first file name and the FAI"
    target = None
    with open(orders, 'r', encoding="utf-8") as fin:
        for line in fin:
            line = line.strip()
            target = line
            break

    asm_name = None
    with open(name_conversion, 'r', encoding="utf-8") as fin:
        for line in fin:
            old_name, new_name = line.strip().split("\t")
            if new_name == target:
                asm_name = old_name
                break

    for fai in fai_list:
        if asm_name in fai:
            return target, fai
    raise ValueError(f"Expected FAI for assembly {asm_name} not found. FAI list: {fai_list}")

rule all:
    input: f"{prefix}_ribbon-plot.png"

rule renaming:
    input: 
        blocks=f"{synteny_blocks}"
    output:
        renamed_blocks=temp(f"{blocks_no_suffix}.renamed.tsv")
    run:
        if name_conversion is not None:
            shell(f"rename_synteny_blocks.py {input.blocks} {name_conversion} > {output.renamed_blocks}")
        else:
            shell(f"ln -s {input.blocks} {output.renamed_blocks}")

rule compute_distances:
    input: 
        rules.renaming.output
    output: 
        distance_tsv = f"{prefix}_est-distances.tsv"
    params: indel_threshold=indel_threshold
    shell:
        "synteny_distance_estimation.py {input} {params.indel_threshold} > {output.distance_tsv}"

rule convert_phylip:
    input: 
        rules.compute_distances.output
    output:
        phylip_format = f"{prefix}_est-distances.phylip"
    shell:
        "convert_distance_matrix_to_phylip.py {input}  > {output.phylip_format}"

rule make_nj_tree:
    input:
        rules.convert_phylip.output
    output:
        nwk = f"{prefix}_est-distances.nwk"
    shell:
        "quicktree -in m -out t {input} > {output.nwk}"

rule cladogram:
    input:
        rules.make_nj_tree.output
    output:
        cladogram = f"{prefix}_est-distances.cladogram.png",
        orders = f"{prefix}_est-distances.order.tsv"
    params:
        prefix = f"{prefix}_est-distances",
        ratio = cladogram_ratio
    shell:
        "distance_cladogram.R --nwk {input} -p {params.prefix} --lim {params.ratio}"

rule sort_blocks:
    input: 
        blocks = f"{blocks_no_suffix}.renamed.tsv",
        orders = rules.cladogram.output.orders
    output:
        sorted_blocks = f"{blocks_no_suffix}.renamed.sorted.tsv"
    run:
        order_blocks = []
        with open(input.orders, 'r') as fin:
            for line in fin:
                order_blocks.append(line.strip())
        order_blocks_str = " ".join(order_blocks)
        shell(f"sort_ntsynt_blocks.py --synteny_blocks {input.blocks} --sort_order {order_blocks_str} > {output.sorted_blocks}")

rule gggenomes_files:
    input: 
        fais = fais,
        orders = rules.cladogram.output.orders,
        blocks = rules.sort_blocks.output.sorted_blocks
    output:
        links = f"{prefix}.links.tsv",
        sequences = f"{prefix}.sequence_lengths.tsv"
    params:
        prefix = prefix
    run:
        first_block = None
        with open(input.orders, 'r') as fin:
            for line in fin:
                first_block = line.strip()
                break
        if name_conversion:
            shell(f"format_blocks_gggenomes.py --fai {sort_fais(input.fais, name_conversion, input.orders)} --prefix {params.prefix} --blocks {input.blocks} --length {min_len} --colour {first_block} --name_conversion {name_conversion}")
        else:
            print("TODO: implement this")
            pass
            # fais_str = " ".join(input.fais)
            # shell(f"format_blocks_gggenomes.py --fai {fais_str} --prefix {params.prefix} --blocks {input.blocks} --length {min_len} --colour {first_block}")

rule chrom_sorting:
    input: 
        fais = fais,
        orders = rules.cladogram.output.orders,
        sequences = rules.gggenomes_files.output.sequences,
        blocks = rules.sort_blocks.output
    output:
        sorted_seqs = f"{prefix}.sequence_lengths.sorted.tsv"
    run:
        if name_conversion:
            fais = sort_fais(input.fais, name_conversion, input.orders)
            target, fai = get_first_asm_and_fai(input.orders, name_conversion, input.fais)
            shell(f"gggenomes_sort_sequences.py --fai {fais} --blocks {input.blocks} --lengths {input.sequences} > {output.sorted_seqs}")
        else:
            fai = input.fais[0]
            target = fai.removesuffix(".fai")
            shell(f"gggenomes_sort_sequences.py --fai {fai} --blocks {input.blocks} --lengths {input.sequences} > {output.sorted_seqs}")


rule ribbon_plot:
    input: 
        links = rules.gggenomes_files.output.links,
        sequences = rules.chrom_sorting.output.sorted_seqs,
        tree = rules.make_nj_tree.output
    output:
        png = f"{prefix}_ribbon-plot.png"
    params:
        prefix = f"{prefix}_ribbon-plot",
        ratio = ribbon_ratio,
        scale = scale
    shell:
        "plot_synteny_blocks_gggenomes_ggtree.R -s {input.sequences} -l {input.links} -p {params.prefix} --tree {input.tree} --ratio {params.ratio} --scale {params.scale}"
