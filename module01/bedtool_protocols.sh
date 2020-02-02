#!/usr/bin/env bash

# name of the gene expression data file; has gene positions and names.
gene_features='942_FPKM_B73_genes_w_feature.txt'
# name of genome annotation file
genome_features='Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.1.gff3.gz'
# name of snp locations file
snp_features='snp_loc.txt'

# name of output bed files
gene_bed='genes.bed'
exon_bed='exons.bed'
exon_merge_bed='exons.merge.bed'
snp_bed='snps.bed'
overlap_bed='overlap.bed'

###############################################################################

# extract gene names and locations
# step 1) remove header from file using tail
# step 2) extract features, apply transformations if needed
# step 3) pipe to disk
tail -n +2 "${gene_features}" | \
    awk '{printf("%s\t%d\t%d\t%s\n","chr"$2,$4-1,$5-1,$1)}' \
    > "${gene_bed}"

# extract exon features
# step 1) extract genome features from archive
# step 2) grep all exons
# step 3) remove all contig features
# step 4) extract features, apply transformations if needed
# step 5) sort the exon file
# step 6) pipe to disk
zcat "${genome_features}" | \
    grep 'exon' | \
    grep -v 'ctg' | \
    awk '{printf("%s\t%d\t%d\n","chr"$1,$4-1,$5-1)}' | \
    sort -k1,1 -k2,2g -k3,3g \
    > "${exon_bed}"

# merge overlapping exons
# step 1) merge overlapping exons
# step 2) pipe to disk
bedtools merge -i "${exon_bed}" \
    > "${exon_merge_bed}"

# extract SNP locations
# step 1) remove the header from the SNP feature file
# step 2) extract features, apply transformations if needed
# step 3) sort the snp file
# step 4) pipe to disk
tail -n +2 "${snp_features}" | \
    awk '{printf("%s\t%d\t%d\t%s\n",$2,$3-1,$3,$1)}' | \
    sort -k1,1 -k2,2g -k3,3g \
    > "${snp_bed}"

# overlap exon and snp bed file to get where genes are
# step 1) intersect preserving all data
# step 2) pipe to disk
bedtools intersect -a "${exon_merge_bed}" -b "${snp_bed}" -wa -wb \
    > "${overlap_bed}"
