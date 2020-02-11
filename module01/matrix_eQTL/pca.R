library(prcomp)

data <- read.table("../input_data_matrixeqtl/Genotype.txt")

prcom <- prcomp(data)

save(prcom,file="../snp_prcom.Rdata")