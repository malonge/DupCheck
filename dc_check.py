#!/usr/bin/env python
import sys
import argparse

import numpy as np
from scipy import stats


def fit_gaussian(cov_dict, trim=0):
    """ Fit a gaussian distribution to the coverage across bins genome-wide. """
    # Later will trim off outliers
    return stats.norm.fit(np.asarray(list(cov_dict.values())).flatten())


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
    glob_mean, glob_std = fit_gaussian(chr_covs)

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
                for tag in fields[8]:
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
                    cov_arr = chr_covs[dup_header]
    
                    # Get the adjusted dup length which is as close to the real dup length as possible
                    # but divisible by the bin size.
                    adj_dup_len = get_adj_dup_len(dup_len, bin_size)
                    print(dup_len, adj_dup_len)

                    # Decide how many windows we should shave off the start to get close to the
                    # real starting point.
                    lop_size = get_lop_size(dup_start, adj_dup_len, bin_size)
                    print("lop_size = %r" % lop_size)

                    # Lop off the required number of bins
                    lop_cov_arr = cov_arr[lop_size:]

                    # Merge bins into the necessary adjusted dup size
                    factor = adj_dup_len//bin_size
                    print("loped coverage array has length= %r" %(len(lop_cov_arr)))

                    # Truncate the coverage array so that it is divisible by the factor
                    trunc_cov_arr = lop_cov_arr[:(len(lop_cov_arr)//factor)*factor]
                    cov_mat = trunc_cov_arr.reshape(len(trunc_cov_arr)//factor, factor)
                    merged_cov_arr = np.average(cov_mat, axis=1)
                    print(len(trunc_cov_arr), len(merged_cov_arr))

                    # Now our coverage is merged and aligned. Let's see the coverage of our duplication window
                    print("DUP bin coverage = %f" %(merged_cov_arr[dup_start//adj_dup_len]))

                    # And the coverage of all of the bins that size
                    print("average coverage of original bins = %f" % glob_mean)
                    print("\n\nmoving on, and getting over\n\n")


if __name__ == "__main__":
    main()
