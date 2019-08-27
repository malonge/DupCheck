#!/usr/bin/env python
import sys
import argparse


def main():
    #parser = argparse.ArgumentParser(description='Check DUPs in a VCF file for increased short-read coverage')
    #parser.add_argument("vcf", metavar="<vars.vcf>", type=str, help="VCF file containing DUPs to check. May contain non-dups.")

    #args = parser.parse_args()
    #vcf_file = args.vcf

    chr = sys.argv[1]
    dup_start = int(sys.argv[2])
    dup_end = int(sys.argv[3])


if __name__ == "__main__":
    main()