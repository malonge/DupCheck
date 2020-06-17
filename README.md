# DupCheck

Given a vcf file and a short-read alignment bam file, identify duplications
associated with increased short read coverage.

This work is a part of a larger panSV-genome analysis. The code is written expressly for 
that analysis, and will not be updated or maintained for future use.

## Dependencies
- python3
- numpy (developed and run on v1.17.2)
- scipy (developed and run on v1.3.1)
- bedtools (v2.17.0)

## Steps
1. Index

```
$ python3 dc_index.py reference.fasta -b 200
```

2. Get coverage
```
$ bedtools coverage -abam alns.bam -b index.bed > index.cov.bed
$ sort -k1,1 -k2,2n -k3,3n index.cov.bed > index.cov.srt.bed
```

3. GC correct coverage
```
$ python3 dc_gccorr.py index.cov.srt.bed > index.cov.corr.bed
```

4. Validate dups
```
$ python3 dc_check.py vars.vcf index.cov.corr.bed
```
