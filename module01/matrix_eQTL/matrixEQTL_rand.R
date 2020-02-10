library(MatrixEQTL)
library(data.table)
args <- commandArgs(T)
print(args)

#4 parameters
#1: cis threhold
#2: trans thresold
#3: cis distance
#4 mode (qt_norm)

#iter <- as.numeric(args[5])
#print("randomizing genotype")
#g <- fread("/mnt/research/compbio/wanglab/haowang/stuff/LaPres/MatrixEQTL/genotype")
#g <- as.data.frame(g)
#g[,-1] <- g[,-1][,sample(1:dim(g[,-1])[2])]
#temp_genotype_pth <- paste0("/mnt/ls15/scratch/users/wangha73/LaPres-PanelSeq/intermediate/randomized_genotype")
#fwrite(g,temp_genotype_pth,sep="\t")
#rm(g)
#gc()
#print("radomize finish")
#genotype file
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile("../input_data_matrixeqtl/Genotype.txt");






rand_ind <- sample(1:942)

for(sl in 1:length(snps)){

  snps[[sl]] <- snps[[sl]][,rand_ind]


}


#expression file
expr = SlicedData$new();
expr$fileDelimiter = "\t";      # the TAB character
expr$fileOmitCharacters = "NA"; # denote missing values;
expr$fileSkipRows = 1;          # one row of column labels
expr$fileSkipColumns = 1;       # one column of row labels
expr$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
expr$LoadFile("../input_data_matrixeqtl/gene_exp.txt")
mode <- args[4]
if(mode=="qt_norm"){
for( sl in 1:length(expr) ) {
  mat = expr[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  mat = qnorm(mat / (ncol(expr)+1));
  expr[[sl]] = mat;
}
rm(sl, mat);
}
#covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile("../pca/pc11.txt")

cisDist = as.numeric(args[3]);

pvOutputThreshold_cis = as.numeric(args[1]);
pvOutputThreshold_tra = as.numeric(args[2]);

output_file_name_cis = paste0("../matrixeqtl_output/cis_pc11_",mode,"_",pvOutputThreshold_cis,"_",cisDist,"_",".txt")
output_file_name_tra = paste0("../matrixeqtl_output/trans_pc11_",mode,"_",pvOutputThreshold_tra,"_",cisDist,"_",".txt")

snpspos = read.table("../input_data_matrixeqtl/snp_loc.txt", header = TRUE, stringsAsFactors = FALSE);


genepos = read.table("../input_data_matrixeqtl/GeneLocation.tsv", header = TRUE, stringsAsFactors = FALSE);

useModel = modelLINEAR

errorCovariance <- numeric()

me = Matrix_eQTL_main(
  snps = snps, 
  gene = expr, 
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = F, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra)
unlink(output_file_name_cis)
save(me,file=paste0("../matrixeqtl_output/cis_pc11_",mode,"_",pvOutputThreshold_cis,"_",pvOutputThreshold_tra,"_",cisDist,"_",args[5],".Rdata"))


#save(me,file=paste0("/mnt/gs18/scratch/users/wangha73/LaPres-PanelSeq/output/rand/correct_pc3_me_out_",mode,"_",pvOutputThreshold_cis,"_",pvOutputThreshold_tra,"_",cisDist,"_",args[5]))
