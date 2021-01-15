library(reshape2)
args <- commandArgs(T)
Coabundance<-read.table(args[1],stringsAsFactors=F,quote="",sep="\t")
colname <- as.character(Coabundance[1,])
names(Coabundance) <- colname
n <- nrow(Coabundance)
Coabundance <- Coabundance[2:n,]
colnames(Coabundance)[1]<-c("OTU_id")
Coabundance<-melt(Coabundance,id.vars="OTU_id",variable.name="OTU2",value.name="SparCC")
names(Coabundance)[1]<-"OTU1"

Pvalues<-read.table(args[2],stringsAsFactors=F,quote ="",sep="\t")
colname <- as.character(Pvalues[1,])
names(Pvalues) <- colname
n <- nrow(Pvalues)
Pvalues <- Pvalues[2:n,]
colnames(Pvalues)[1]<-c("OTU_id")
Pvalues<-melt(Pvalues,id.vars="OTU_id",variable.name="OTU2",value.name="p")
names(Pvalues)[1]<-"OTU1"

Coabundance<-merge(Coabundance,Pvalues,by=c("OTU1","OTU2"),all.x=TRUE)
Coabundance$SparCC <- as.numeric(Coabundance$SparCC)
Coabundance$p <- as.numeric(Coabundance$p)
p<-Coabundance$p
qobj<-p.adjust(p,"fdr")
Coabundance$qvalue<-qobj
Coabundance<-Coabundance[Coabundance$qvalue<0.05,]
Coabundance<-subset(Coabundance,abs(Coabundance$SparCC)>0.6)

write.table(Coabundance,args[3],sep='\t',row.names=F,quote = F)
