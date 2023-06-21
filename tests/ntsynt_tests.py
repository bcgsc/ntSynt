#!/usr/bin/env python3
'''
Tests for ntSynt
'''
import shlex
import subprocess

def launch_ntSynt(*genomes, prefix="ntSynt", k=24, w=1000):
    "Launch ntSynt"
    genome_list = " ".join(genomes)
    cmd = f"ntSynt.py --force {genome_list} -k{k} -w {w} --prefix {prefix}"
    cmd = shlex.split(cmd)
    ret_code = subprocess.call(cmd)
    assert ret_code == 0

def are_expected_blocks(block_file1, block_file2):
    "Check that the given synteny blocks are as expected"
    with open(block_file1, 'r', encoding="utf-8") as fin1:
        with open(block_file2, 'r', encoding="utf-8") as fin2:
            while (line1 := fin1.readline()) and (line2 := fin2.readline()):
                assert line1 == line2

def test_prep_files():
    "Prep the test files"
    infiles = "celegans-chrII-III.fa.gz", "celegans-chrII-III.A.fa.gz", "celegans-chrII-III.B.fa.gz"
    for infile in infiles:
        cmd = shlex.split(f"unpigz {infile}")
        ret_code = subprocess.call(cmd)
        assert ret_code == 0

def test_2_genomes():
    "Testing ntSynt with two input genomes"
    genome1, genome2 = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa"
    launch_ntSynt(genome1, genome2, k=24, prefix="celegans-A-ntSynt")
    are_expected_blocks("celegans-A-ntSynt.synteny_blocks.tsv",
                        "expected_result/celegans-A-ntSynt.synteny_blocks.tsv")

def test_3_genomes():
    "Testing ntSynt with three input genomes"
    genome1, genome2, genome3 = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa", "celegans-chrII-III.B.fa"
    launch_ntSynt(genome1, genome2, genome3, k=20, prefix="celegans-A-B-ntSynt")
    are_expected_blocks("celegans-A-B-ntSynt.synteny_blocks.tsv",
                        "expected_result/celegans-A-B-ntSynt.synteny_blocks.tsv")

def test_cleanup():
    "Re-compress the files"
    infiles = "celegans-chrII-III.fa", "celegans-chrII-III.A.fa", "celegans-chrII-III.B.fa"
    for infile in infiles:
        cmd = shlex.split(f"pigz {infile}")
        ret_code = subprocess.call(cmd)
        assert ret_code == 0
