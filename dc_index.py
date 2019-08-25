#!/usr/bin/env python
import argparse

from dupcheck_utilities.utilities import gc_content, log
from dupcheck_utilities.SeqReader import SeqReader


def main():
    parser = argparse.ArgumentParser(description='Creat and index for validating dups')
    parser.add_argument("reference", metavar="<in.fasta>", type=str, help="Reference genome in fasta format.")
    parser.add_argument("-b", metavar="<bin_size>", type=int, help="Bin size")

    args = parser.parse_args()
    reference_file = args.reference
    bin_size = args.b

    # Parse the reference and make the GC bed file
    seq_lens = dict()
    x = SeqReader(reference_file)
    all_bed_intervals = list()
    for header, seq in x.parse_fasta():
        seq_lens[header] = len(seq)

        for i in range(0, seq_lens[header], bin_size):
            if i%1000000 == 0:
                log("Processed %r bp" %i)

            this_line = (header, i, i+bin_size, int(gc_content(seq[i:i + bin_size])*100))
            all_bed_intervals.append(this_line)

        # Remove the last bin if it is not complete
        if all_bed_intervals[-1][2] > seq_lens[header]:
            all_bed_intervals = all_bed_intervals[:-1]

    all_bed_intervals = sorted(all_bed_intervals)
    all_bed_str = []
    for i in all_bed_intervals:
        all_bed_str.append("\t".join([str(j) for j in i]))

    with open('index.bed', 'w') as f:
        f.write("\n".join(all_bed_str) + "\n")


if __name__ == "__main__":
    main()
