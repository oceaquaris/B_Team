#scan n mers
#x input sequence
#n : width
args <- commandArgs(T)

n_length <- as.numeric(args[2])

file <- read.table(args[1],colClass="character")

kmers <- function(x,n){
  substring(x,1:(nchar(x)-n+1),n:nchar(x))
}

kmers <- Vectorize(kmers,vectorize.args="x")

dis_kmers <- kmers(file,n_length)

write.table(dis_kmers,"./kmers_discovered.txt",col.names=F,row.names=F,sep="\t",quote=F)


