#!/usr/bin/env snakemake -s
import shutil
import sys
import os


# Parameters
synteny_blocks = config["blocks"] if "blocks" in config else None
name_conversion = config["name_conversion"] if "name_conversion" in config else None
name_conversion = None if name_conversion == "None" else name_conversion
prefix = config["prefix"] if "prefix" in config else "ntSynt_gggenomes"
indel_threshold = config["indel_threshold"] if "indel_threshold" in config else 50000
min_len = config["min_length"] if "min_length" in config else 50000
fais = config["fai"] if "fai" in config else None
ribbon_ratio = config["ribbon_ratio"] if "ribbon_ratio" in config else 0.1
cladogram_ratio = config["cladogram_adjust"] if "cladogram_adjust" in config else 0.1
scale = config["scale"] if "scale" in config else 1e9
centromeres = config["centromeres"] if "centromeres" in config else None
normalize = config["normalize"] if "normalize" in config else False
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

def sort_fais_no_name_conversion(fai_list, orders):
    "Sort the FAIS when have no name conversion file"
    # Make order dictionary
    order_dict = {}
    with open(orders, 'r') as fin:
        i = 0
        for line in fin:
            order_dict[line.strip()] = i
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
        orders = rules.cladogram.output.orders,
        fais = fais
    output:
        sorted_blocks = f"{blocks_no_suffix}.renamed.sorted.blocks.tsv"
    params:
        intermediate_blocks = f"{blocks_no_suffix}.renamed.sorted-tmp.tsv",
        name_conversion = f"-c {name_conversion}" if name_conversion else "",
        prefix = f"{blocks_no_suffix}.renamed.sorted"
    run:
        order_blocks = []
        with open(input.orders, 'r') as fin:
            for line in fin:
                order_blocks.append(line.strip())
        order_blocks_str = " ".join(order_blocks)
        if normalize:
            shell(f"sort_ntsynt_blocks.py --synteny_blocks {input.blocks} --sort_order {order_blocks_str} > {params.intermediate_blocks}")
            if name_conversion:
                shell(f"gggenomes_normalize_strands.py --blocks {params.intermediate_blocks} --fais {sort_fais(input.fais, name_conversion, input.orders)} -p {params.prefix} {params.name_conversion}")
            else:
                shell(f"gggenomes_normalize_strands.py --blocks {params.intermediate_blocks} --fais {sort_fais_no_name_conversion(input.fais, input.orders)} -p {params.prefix} {params.name_conversion}")
        else:
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
            shell(f"format_blocks_gggenomes.py --fai {sort_fais_no_name_conversion(input.fais, input.orders)} --prefix {params.prefix} --blocks {input.blocks} --length {min_len} --colour {first_block}")


rule chrom_sorting:
    input: 
        fais = fais,
        orders = rules.cladogram.output.orders,
        sequences = rules.gggenomes_files.output.sequences,
        blocks = rules.sort_blocks.output.sorted_blocks
    output:
        sorted_seqs = f"{prefix}.sequence_lengths.sorted.tsv"
    run:
        if name_conversion:
            fais = sort_fais(input.fais, name_conversion, input.orders)
        else:
            fais = sort_fais_no_name_conversion(input.fais, input.orders)
        shell(f"gggenomes_sort_sequences.py --fai {fais} --blocks {input.blocks} --lengths {input.sequences} > {output.sorted_seqs}")

rule chrom_paint:
    input: links = rules.gggenomes_files.output.links
    output: colour_feats = f"{prefix}.chrom-paint-feats.tsv"
    shell:
        '''cat {input.links} |perl -ne 'chomp; @a=split("\t"); if(!defined $ct){{print "block_id\tseq_id\tbin_id\tstart\tend\tcolour_block\n"; $ct=1}} else {{print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[11]\n"; print "$a[0]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[11]\n";}}' > {output.colour_feats}'''

rule ribbon_plot:
    input: 
        links = rules.gggenomes_files.output.links,
        sequences = rules.chrom_sorting.output.sorted_seqs,
        tree = rules.make_nj_tree.output,
        colour_feats = rules.chrom_paint.output.colour_feats
    output:
        png = f"{prefix}_ribbon-plot.png"
    params:
        prefix = f"{prefix}_ribbon-plot",
        ratio = ribbon_ratio,
        scale = scale,
        centromeres = f"--centromeres {centromeres}" if centromeres is not None else ""
    shell:
        "plot_synteny_blocks_gggenomes_ggtree.R -s {input.sequences} -l {input.links} -p {params.prefix} --tree {input.tree} --ratio {params.ratio} --scale {params.scale} -c {input.colour_feats} {params.centromeres}"
