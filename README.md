# DupCheck

Given a vcf file and a short-read alignment bam file, identify duplications
associated with increased short read coverage.

## Steps
1. Index

`$ python3 dc_index.py reference.fasta -b <bin_size>`

2. Get coverage

`$ bedtools coverage -abam alns.bam -b bins.bed`

3. Validate dups
`$ python3 dc_validate.py`