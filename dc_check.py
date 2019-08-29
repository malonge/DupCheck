#!/usr/bin/env python
import sys
import argparse

import numpy as np
from scipy import stats

from dupcheck_utilities.utilities import log


def remove_iqr_outliers(x):
    """ From https://medium.com/datadriveninvestor/finding-outliers-in-dataset-using-python-efc3fce6ce32"""
    x = np.sort(x)
    q1, q3 = np.percentile(x,[25,75])
    iqr = q3 - q1
    lower_bound = q1 - (1.5 * iqr)
    upper_bound = q3 + (1.5 * iqr)
    x = x[x > lower_bound]
    return x[x < upper_bound]


def fit_gaussian(cov_dict):
    """ Fit a gaussian distribution to the coverage across bins genome-wide. """
    x = np.concatenate(list(cov_dict.values()))
    x = remove_iqr_outliers(x)
    return stats.norm.fit(x)


def get_adj_dup_len(dup_len, bin_size):
    new_dup_len = (dup_len//bin_size)*bin_size
    if dup_len%bin_size > bin_size//2:
        new_dup_len += bin_size

    return new_dup_len


def get_lop_size(dup_start, adj_dup_len, bin_size):
    lop_size = (dup_start % adj_dup_len)//bin_size

    # If we are more than half of the length away, better to extend past the start.
    if ((dup_start % adj_dup_len) % bin_size) > (bin_size//2):
        lop_size += 1
    return lop_size


def main():
    parser = argparse.ArgumentParser(description='Check DUPs in a VCF file for increased short-read coverage')
    parser.add_argument("vcf", metavar="<vars.vcf>", type=str, help="VCF file containing DUPs to check. May contain non-dups.")
    parser.add_argument("cov", metavar="<coverage.bed>", type=str, help="Bed file from dc_index-bedtools-dc_gccorr")

    # Get command line args
    args = parser.parse_args()
    vcf_file = args.vcf
    cov_file = args.cov
    bin_size = 200

    # Read the coverage file
    log("Reading {}".format(cov_file))
    chr_covs = dict()
    with open(cov_file, "r") as f:
        for line in f:
            header, start, end, gc, a, b, c, d, corr_cov = line.rstrip().split("\t")
            if header not in chr_covs:
                chr_covs[header] = []
            chr_covs[header].append(int(round(float(corr_cov))))

    for i in chr_covs:
        chr_covs[i] = np.asarray(chr_covs[i], dtype=np.int32)

    # Get the average coverage across the initial bins
    log("Calculating global average coverage.")
    glob_mean, glob_std = fit_gaussian(chr_covs)
    log("After removing outliers, the global average coverage is {} X.".format(glob_mean))

    # Process each line of the vcf file. If its a comment or non dup, just print it.
    with open(vcf_file, "r") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("#"):
                print(line)
            else:
                fields = line.split("\t")
                sv_type = ''
                sv_len = 0
                tags = fields[7].split(";")
                for tag in tags:
                    if tag.startswith("SVTYPE="):
                        sv_type = tag[7:]
                    if tag.startswith("SVLEN="):
                        sv_len = int(tag[6:])

                assert sv_type
                assert sv_len
                if sv_type != "DUP":
                    print(line)
                else:
                    # This is a duplication, so lets try to verify it.
                    dup_header = fields[0]
                    dup_start = int(fields[1])
                    dup_end = dup_start + sv_len
                    dup_len = sv_len
                    if dup_len < 1000:
                        continue

                    cov_arr = chr_covs[dup_header]
                    log("Processing a DUP at {}:{}-{}.".format(dup_header, dup_start, dup_end))

                    # Get the adjusted dup length which is as close to the real dup length as possible
                    # but divisible by the bin size.
                    adj_dup_len = get_adj_dup_len(dup_len, bin_size)
                    len_diff = abs(dup_len - adj_dup_len)
                    log("The original DUP length is {} and the DUP bin size is {}.".format(dup_len, adj_dup_len))
                    log("The difference in length is {} bp".format(len_diff))

                    # Decide how many windows we should shave off the start to get close to the
                    # real starting point.
                    lop_size = get_lop_size(dup_start, adj_dup_len, bin_size)

                    # Lop off the required number of bins
                    lop_cov_arr = cov_arr[lop_size:]

                    # Merge bins into the necessary adjusted dup size
                    factor = adj_dup_len//bin_size

                    # Truncate the coverage array so that it is divisible by the factor
                    trunc_cov_arr = lop_cov_arr[:(len(lop_cov_arr)//factor)*factor]
                    cov_mat = trunc_cov_arr.reshape(len(trunc_cov_arr)//factor, factor)
                    merged_cov_arr = np.average(cov_mat, axis=1)
                    dup_coverage = merged_cov_arr[dup_start//adj_dup_len]

                    # Now our coverage is merged and aligned. Let's see the coverage of our duplication window
                    log("The DUP bin coverage = {}".format(dup_coverage))

                    # And the coverage of all of the bins that size
                    log("This DUP is {}X the average coverage.".format((dup_coverage/glob_mean)))
                    if dup_coverage/glob_mean > 1.75:
                        print(line)
                        log("Keeping DUP.")
                    else:
                        log("Discarding DUP.")


if __name__ == "__main__":
    main()
