library(reshape2)
args <- commandArgs(T)
name_rows <- function(df){
  tmp <- t(apply(df[,1:2],1,sort))
  index <- !duplicated(tmp)
  df <- df[index,]
  a <- apply(tmp,1,function(x)paste(x[1],x[2], sep = '&'))
  row.names(df) <- a[index]
  return(df)
}
int_df <- name_rows(read.csv(args[1], sep = '\t'))  #control interaction dataframe
pairs <- row.names(int_df)
r0 <- int_df[pairs,3]
diff_mx <- matrix(0, nrow = 100, ncol = length(pairs))
colnames(diff_mx) <- pairs
for (i in 0:99){
  boot_file <- paste(paste('con_combined_genus_cor',i,sep = '_'),'txt',sep = '.')
##boot_file is generated in SPARCC
  co_abundance <- read.csv(boot_file, sep = '\t',stringsAsFactors = F) #attention!blank will replaced with '.' in colnames in read table/csv
  colnames(co_abundance) <- c('OTU',co_abundance[,1])
  co_abundance2 <- melt(co_abundance) 
  names(co_abundance2) <- c('OTU1', 'OTU2', 'SparCC')
  co_abundance2 <- name_rows(co_abundance2) #with duplications removed
  r1 <- co_abundance2[pairs,3]
  diff <- abs(r0-r1)
  diff_mx[i+1,] <- diff
}
write.table(diff_mx,'diff_permutation.txt', sep = '\t', row.names = F, quote = F)

int_df2 <- name_rows(read.csv(args[2], sep = '\t')) #case interaction dataframe
pairs2 <- row.names(int_df2)

p_value <- function(pair){
  h0 <- as.numeric(diff_mx[,pair])
  r0 <- int_df[pair,3]
  r1 <- int_df2[pair,3]
  diff <- abs(r0-r1)
  p <- length(h0[h0>diff])/length(h0)
  return(p)
}

specific_int1 <- setdiff(pairs2,pairs)
union_int <- intersect(pairs, pairs2)
p_union_int <- sapply(union_int,p_value)
specific_int2 <- names(p_union_int[p_union_int<=0.05])
out_df <- int_df2[c(specific_int1,specific_int2),]
write.table(out_df,args[3], sep = '\t', quote = F, row.names = F)



