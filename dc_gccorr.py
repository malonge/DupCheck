#!/usr/bin/env python
import argparse


def main():
    bed_help = """
    Output bed file from bedtools coverage. The input bedfile for bedtools should be that produced by dc_index.py.
    """
    parser = argparse.ArgumentParser(description='Perform GC bias correction on observed raw read coverage')
    parser.add_argument("bed", metavar="<in.bed>", type=str, help=bed_help)

    args = parser.parse_args()
    bed_file = args.bed

    # First pass through the bed file
    # Calculate global average and GC specific average
    gc_total = [0 for i in range(101)]  # Sum of all coverage for the 100 possible GC values
    gc_count = [0 for i in range(101)]  # Count of all bins with a given GC value
    total = 0
    count = 0
    with open(bed_file, "r") as f:
        for line in f:
            bed_fields = line.rstrip().split("\t")
            total += int(bed_fields[4])
            count += 1

            # GC specific
            gc_total[int(bed_fields[3])] += int(bed_fields[4])
            gc_count[int(bed_fields[3])] += 1

    bin_avg = total/count
    gc_avg = []
    for i, j in zip(gc_total, gc_count):
        if not j:
            gc_avg.append(0)
        else:
            gc_avg.append(i/j)

    # Second pass through the bed file
    # For each coverage value, calculate the GC corrected coverage
    with open(bed_file, "r") as f:
        for line in f:
            bed_fields = line.rstrip().split("\t")
            bed_fields.append(str(int(bed_fields[4]) * (bin_avg/gc_avg[int(bed_fields[3])])))
            print("\t".join(bed_fields))


if __name__ == "__main__":
    main()
