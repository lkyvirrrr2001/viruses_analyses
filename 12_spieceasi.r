args <- commandArgs(T)
library(SpiecEasi)
df <- read.csv(args[1], sep = '\t', row.names = 1) # input abundance matrix
data=t(as.matrix(df))
pargs=list(seed=1001,ncores=20,rep.num=30)
spiec = spiec.easi(as.matrix(data), method='glasso',sel.criterion = "stars",lambda.min.ratio=1e-3, nlambda=50,pulsar.params=pargs)
cor=cov2cor(apply(getOptCov(spiec),2,as.numeric))
row.names(cor)=colnames(data)
colnames(cor)=colnames(data)
cor2 <- as.matrix(cor*getRefit(spiec))
cor2[lower.tri(cor2, diag = T)]=NA
ind = which(is.na(cor2)==F&cor2!=0,arr.ind = T)
result = data.frame(start = colnames(cor2)[ind[,1]],end = colnames(cor2)[ind[,2]],cor = cor2[ind],stringsAsFactors = F)
write.table(result,args[2], sep = '\t', quote = F, row.names = F)
