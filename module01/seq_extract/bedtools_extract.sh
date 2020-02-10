#!/usr/bin/env bash

# arguments for this script
# bedtools_extract <snps.bed> <input file> <range> <reference.fa>
snps_bed=$1     # get bed file
infile=$2       # get input file
ntrange=$3      # get nucleotide range
reference=$4    # get reference genome

# if you want to manually set the aboves, reassign below.
# snps_bed="input.bed"
# infile="blah.txt"
# ntrange=20
# reference="reference.fa"

# get base name for infile
outbase="${infile%.*}"

# extract sequences for each significant SNP
# step 1) grep significant SNPs from our SNPs.bed file
# step 2) add +/- our nucleotide range, ntrange
# step 3) sort the bed file just in case
# step 4) merge overlapping regions so we don't double sample
# step 5) replace 'chr' with 'Chr' because reference genome is dumb.
# step 6) extract sequences from the reference genome and output to an outfile.
grep -f "${infile}" "${snps_bed}" | \
    awk -v b=$ntrange '{printf("%s\t%d\t%d\t%s\n",$1,$2-b,$3+b,$4)}' | \
    sort -k1,1 -k2,2g -k3,3g | \
    bedtools merge -i stdin | \
    sed 's/chr/Chr/g' | \
    bedtools getfasta -name -fi "${reference}" -bed stdin -fo "${outbase}.fa"
