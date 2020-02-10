import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
%matplotlib inline

#read in the file generated from bedtools on the command line
#this analysis includes introns
#this file includes the following information: chromosome, gene start and stop positions, gene ID, SNP start and stop positions, and SNP ID (there are 2 columns of chromosome number, one for genes and one for SNPs)
overlap = pd.read_csv('overlap.bed',sep="\t",header=None)
#view the first 5 lines of the file
overlap.head()
#group by the fourth column, which is gene ID
genes=overlap.groupby(3)

#loop for finding the distance between SNPs
#first define an empty list to contain the distances
snpdist = []
#now iterate through the file
for name, gene in genes:
    #sort SNPs by position
    gene.sort_values(by=[5],inplace=True)
    #list comprehension loop for calculating the distances. The sixth column (5) is SNP start, and the seventh column (6) is SNP stop position.
    snpdist += [gene.iloc[i,5]-gene.iloc[i+1,5] for i in range(0,len(gene[6])-1)] #len(gene[6])-1 is used because we're subtracting, so if there are n SNPs, n-1 distances will be calculated
#view resulting list
snpdist

#plot histogram of values
plt.hist(snpdist,bins=1000)
plt.xlim(-2000,0)

#descriptive statistics
np.median(snpdist)
np.std(stpdist) #standard deviation
np.percentile(snpdist,75)-np.percentile(snpdist,25) #interquartile range

#now repeat the analysis with only exons
#alternative splice forms for each gene were merged using bedtools on the command line
exons = pd.read_csv('overlap_exons.bed',sep="\t",header=None)
exons.head()
exongenes = exons.groupby([0,1])

exon_snpdist = []
for name, gene in exongenes:
    #sort SNPs by position
    gene.sort_values(by=[4],ascending=False,inplace=True)
    exon_snpdist += [gene.iloc[i,4]-gene.iloc[i+1,4] for i in range(0,len(gene[4])-1)]
#plot the new histogram
plt.hist(exon_snpdist,bins=200)
#descriptive statistics
np.median(exon_snpdist)
np.std(exon_snpdist)
np.percentile(exon_snpdist,75)-np.percentile(exon_snpdist,25)

#find out how many SNPs are within 4 bp of each other
absdist=np.abs(exon_snpdist)
(absdist<=4).sum()
#find out how many total SNP distances within the same exon were calculated
len(absdist)
