#!/usr/bin/env python3 
import argparse
import sys


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("-m", "--mirnas", type=str, required=True, help="path to the mirnas file")
    args.add_argument("-o", "--output", type=str, required=True, help="path to the output file")
    args = args.parse_args()

    with open(args.mirnas, "r") as f:
        mirna_counts = {mirna: 0 for mirna in f.read().split()}

    total_reads = 0
    for line in sys.stdin:
        cols = line.split('\t')
        if len(cols) >= 15:
            col = cols[14]
            mirnas = [mirna.strip() for mirna in col.split(":")[-1].split(",")]
            for mirna in mirnas:
                total_reads += 1
                if mirna in mirna_counts:
                    mirna_counts[mirna] += 1
    print("Total reads: ", total_reads)
    print("Matched reads: ", sum(mirna_counts.values()))
    print("Matched mirnas: ", sum([value != 0 for key, value in mirna_counts.items()]))
    print("Unmatched mirnas: ", sum([value == 0 for key, value in mirna_counts.items()]))
    if not any(mirna_counts.values()):
        print("No matches!")
    with open(args.output, "w+") as f:
        for mirna, count in mirna_counts.items():
            f.write(f"{mirna}\n{count}\n")
