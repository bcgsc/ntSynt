#!/usr/bin/env python3
'''
Make a BF of k-mers with multiplicity of 2+ in either input genome
'''
import argparse
import btllib
import ntjoin_utils

def main():
    "Generate the Bloom filters"
    parser = argparse.ArgumentParser(description="Generating BF of k-mer 2+ multiplicities")
    parser.add_argument("--genome", help="Input genome file(s)", nargs="+")
    parser.add_argument("-k", help="K-mer size (bp)", required=True, type=int)
    parser.add_argument("--bf", help="Bloom filter size [accepted units: B (bytes), k (kilobytes)," \
                        " M (megabytes), G (gigabytes)]",
                        type=str)
    parser.add_argument("-t", help="Number of threads [4]", required=False, type=int, default=4)
    parser.add_argument("-p", help="Prefix for output BF", default="out.bf", type=str)
    parser.add_argument("--fpr", help="False positive rate for Bloom filter." \
                    "Only used if --bf is not specified. [0.01]",
                    default=0.01, type=float)

    args = parser.parse_args()

    if not args.bf:
        bf_bytes = ntjoin_utils.approximate_bf_size(args.genome[0], args.fpr, args.t)
    else:
        bf_bytes = ntjoin_utils.parse_bf_size(args.bf, parser)

    rep_bf = btllib.KmerBloomFilter(bf_bytes, 1, args.k)

    for genome in args.genome:
        genome_bf = btllib.KmerBloomFilter(bf_bytes, 1, args.k)
        with btllib.SeqReader(genome, btllib.SeqReaderFlag.LONG_MODE, args.t) as reader:
            for record in reader:
                hasher = btllib.NtHash(record.seq, 1, args.k)
                while hasher.roll():
                    if genome_bf.contains(hasher.hashes()):
                        rep_bf.insert(hasher.hashes())
                    else:
                        genome_bf.insert(hasher.hashes())

    rep_bf.save(f"{args.p}.bf")


if __name__ == "__main__":
    main()
