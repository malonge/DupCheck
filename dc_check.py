#!/usr/bin/env python
import sys
import argparse

import numpy as np

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
    #parser = argparse.ArgumentParser(description='Check DUPs in a VCF file for increased short-read coverage')
    #parser.add_argument("vcf", metavar="<vars.vcf>", type=str, help="VCF file containing DUPs to check. May contain non-dups.")

    #args = parser.parse_args()
    #vcf_file = args.vcf

    dup_start = int(sys.argv[2])
    dup_end = int(sys.argv[3])
    cov_file = sys.argv[1] # For now, assumes its just for a single chromosomes

    dup_len = dup_end - dup_start
    bin_size = 200

    # Read in the array of bin coverages
    covs = []
    with open(cov_file, 'r') as f:
        for line in f:
            L1 = line.rstrip().split("\t")
            covs.append(int(round(float(L1[8]))))
    cov_arr = np.asarray(covs, dtype=np.int32)
    
    # Get the adjusted dup length which is as close to the real dup length as possible
    # but divisible by the bin size.
    adj_dup_len = get_adj_dup_len(dup_len, bin_size)
    print(dup_len, adj_dup_len)

    # Decide how many windows we should shave off the start to get close to the
    # real starting point.
    lop_size = get_lop_size(dup_start, adj_dup_len, bin_size)
    print("lop_size = %r" %lop_size)

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
    print("average coverage of bins of the same size = %f" %(np.average(merged_cov_arr)))
if __name__ == "__main__":
    main()
