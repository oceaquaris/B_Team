
file <- dir()

file <- file[grep("pc11_qt_norm_0.001_0",file)]


bg_p_val <- list()

for(i in 1:1000){
  print(i)
  load(file[i])
  
  bg_p_val[[i]] <- me$cis$min.pv.gene
  
  
}

bg_p_val <- do.call(rbind,bg_p_val)

#calculate pt
# find_threshold <- function(x){
#   
#   id <- which.min(abs(x$FDR-0.1))
#   
#   pt <- x$pval[id]
#   
#   return(pt)
#   
# }
# 
# 
# load(file[1001])
# 
# cis_eqtl <- me$cis$eqtls
# 
# gene_wise_pt <- split(cis_eqtl,cis_eqtl$gene)
# 
# gene_wise_pt <- sapply(gene_wise_pt,find_threshold)
# 
load(file[1001])

cis_eqtl <- me$cis$eqtls

id <- match(cis_eqtl$gene,colnames(bg_p_val))

n <- dim(cis_eqtl)[1]

p_adj <- c()

for(i in 1:n){
  
  p_adj[i] <- (1+sum(bg_p_val[,id[i]]<cis_eqtl$pval[i]))/1001
  
  
  
}
