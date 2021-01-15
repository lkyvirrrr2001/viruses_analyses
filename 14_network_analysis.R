##case specific interactions sparCC&Spieceasi
library(igraph)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(plyr)
library(RColorBrewer)
library(ggsci)
managedf <- function(df){
  n <- nrow(df)
  lab <- rep("p",n)
  lab[which(df[,3]<0)] <- "n"
  df <- data.frame(df,lab=lab)
  df$lab <- as.character(df$lab) 
  return(df)
}
name_rows <- function(df){
  tmp <- t(apply(df[,1:2],1,sort))
  #  a <- apply(tmp,1,function(x)paste(paste(x[1],x[2], sep = '&'),x[4], sep = '&'))
  a <- apply(tmp,1,function(x)paste(x[[1]],x[2], sep = '&'))
  b <- data.frame(a,df[,4])
  c <- apply(b,1, function(x)paste(x[1], x[2], sep = '&'))
  row.names(df) <- c
  return(df)
}
comp_d <- function(char,dfname){
  df <- get(dfname)
  tmpdf1 <- df[df[,"OTU1"] %in% char|df[,"OTU2"] %in% char, ] 
  if(is.null(tmpdf1)){
    d <- 0
  } else{
    d <- nrow(tmpdf1)
  }
  
  return(d)
}
n_link_ratio <- function(char,dfname){
  df <- get(dfname)
  tmpdf1 <- df[df[,"OTU1"] %in% char|df[,"OTU2"] %in% char, ] 
  total_link <- nrow(tmpdf1)
  n_link <- nrow(subset(tmpdf1, lab == 'n'))
  ratio <- n_link/total_link
  return(ratio)
}
count_d <- function(d,label){
  uniq_d = data.frame(table(d))
  f_d <- uniq_d$Freq/sum(uniq_d$Freq)
  uniq_d$d <- as.numeric(uniq_d$d)
  uniq_d$labels <- label
  return(uniq_d)
}

compCC <- function(char,dfname){
  df <- get(dfname)
  tmpdf1 <- df[df[,"OTU1"] %in% char|df[,"OTU2"] %in% char, ] 
  if(is.null(tmpdf1)){
    CC <- 1
  } else{
    d <- nrow(tmpdf1)
    t2 <- union(tmpdf1$OTU1,tmpdf1$OTU2)
    t2 <- t2[-which(t2==char)]
    tmpdf2 <- df[df[,"OTU1"] %in% t2 & df[,"OTU2"] %in% t2,]
    c <- nrow(tmpdf2)
    CC <- (d*(d-1)+1)/(2*c+1)
  }
  
  return(CC)
}
scale_free_index <- function(d){
  uniq_d = data.frame(table(d))
  f_d <- uniq_d$Freq/sum(uniq_d$Freq)
  uniq_d$d <- as.numeric(uniq_d$d)
  model <- summary(lm(log(f_d) ~ log(uniq_d$d)))
  return(model$adj.r.squared)
}
distribution <- function(d){
  uniq_d = data.frame(table(d))
  f_d <- uniq_d$Freq/sum(uniq_d$Freq)
  uniq_d$d <- as.numeric(uniq_d$d)
  log_d <- log(uniq_d$d)
  log_f_d <- log(f_d)
  df <- data.frame(log_d,log_f_d)
  return(df)
}

modul <- function(df){
  graph1 <- graph_from_data_frame(df,directed = F)
  wtc <- cluster_walktrap(graph1)
  mod <- modularity(wtc)
  return(mod)
}

path_length <- function(df){
  graph1 <- graph_from_data_frame(df,directed = F)
  l <- mean_distance(graph1,directed = F)
  return(l)
}

count_edges <- function(df,label1,label2) {
  df2 <- subset(df,(OTU1 %in% label1 & OTU2 %in% label2) |(OTU1 %in% label2 & OTU2 %in% label1))
  p_edge <- nrow(subset(df2,lab == 'p'))
  n_edge <- nrow(subset(df2,lab == 'n'))
  return(c(p_edge,n_edge))
}
classify_edges <- function(arow, bac_g, phage_g,eu_g) {
  if (arow[1] %in% bac_g & arow[2] %in% bac_g){
    edge_type <- 'bacteria-bacteria'
  } else if((arow[1] %in% bac_g & arow[2] %in% phage_g)|arow[1] %in% phage_g & arow[2] %in% bac_g){
    edge_type <- 'bacteria-phages'
  } else if((arow[1] %in% bac_g & arow[2] %in% eu_g)|arow[1] %in% eu_g & arow[2] %in% bac_g){
    edge_type <- 'bacteria-euk viruaes'
  } else if((arow[1] %in% eu_g & arow[2] %in% phage_g)|arow[1] %in% phage_g & arow[2] %in% eu_g){
    edge_type <- 'phages-euk viruses'
  } else if (arow[1] %in% phage_g & arow[2] %in% phage_g){
    edge_type <- 'phages-phages'
  } else if(arow[1] %in% eu_g & arow[2] %in% eu_g){
    edge_type <- 'euk viruses-euk viruses'
  } else{
    edge_type <- 'idontknow'
  }
  return(edge_type)
}
neg_bac <- function(df, bac_genus){  #count frequencies of bacteria that have negative correlations
  a <- subset(df,lab == 'n')
  tmp <- c(a$OTU1,a$OTU2)
  n_bac <- tmp[tmp %in% bac_genus]
  freq_df <- as.data.frame(table(n_bac))
  freq_df <- freq_df[order(freq_df$Freq,decreasing = T),]
  return(freq_df)
}

initial_matrix <- function(node_list, fill_value = -1){
  tmp <- as.data.frame(table(node_list))
  tmp <- tmp[order(tmp$Freq, decreasing = T),]
  
  phage_union <- as.character(tmp$node_list)
  n <- list(phage_union, disease2)
  phage_matrix <- matrix(fill_value,nrow = length(phage_union), ncol = 7, dimnames = n)
  return(phage_matrix)
}

initial_matrix2 <- function(df){
  
  tmp <- as.data.frame(table(node_list))
  tmp <- tmp[order(tmp$Freq, decreasing = T),]
  
  phage_union <- as.character(tmp$node_list)
  n <- list(phage_union, disease2)
  phage_matrix <- matrix(fill_value,nrow = length(phage_union), ncol = 7, dimnames = n)
  return(phage_matrix)
}

vir_module <- function(mem){
  mod <- unique(mem)
  mod_num <- c()
  for (num in mod){
    sub_node <- names(mem[mem==num])
    a <- sub_node[sub_node %in% virus_all]
    l1 <- length(sub_node)
    l2 <- length(a)
    if (l1>=8 & l2/l1>=0.7){  #l1:5 or 8 are the same result
      mod_num <- c(mod_num,num)
    }
  }
  return(mod_num)
}
inter_spiec <- function(file, file2){
  df0 <- read.csv(file, sep = '\t', stringsAsFactors = F)
  df0 <- df0[,c(1,2,3)]
  names(df0) <- c('OTU1','OTU2', 'cor')
  df <- managedf(df0)
  df <- name_rows(df)

  df_spiec <-  read.csv(file2, sep = '\t')
  df_spiec2 <- managedf(df_spiec)
  df_spiec2 <- name_rows(df_spiec2)
  
  common_pairs <- intersect(row.names(df), row.names(df_spiec2))
  df <- df[common_pairs,]
  #  df_mx <- df_tmp[common_pairs,]
  return(df)
}
plot_mx <- function(df){
  g <- union(df$OTU1,df$OTU2)
  n <- length(g)
  graph1 <- graph_from_data_frame(df,directed = F)
  wtc <- cluster_walktrap(graph1)
  wtc2 <- sort(membership(wtc))
  mx_name <- list(names(wtc2), names(wtc2))
  mx <- matrix(0,nrow = n, ncol = n, dimnames = mx_name)
  l = nrow(df)
  for (j in 1:l){
    otu1 = df[j,1]
    otu2 = df[j,2]
    mx[otu1,otu2] <- df[j,3]
    mx[otu2,otu1] <- df[j,3]
  }
  return(mx)
}
genus2phylum <- read.table('genus2phylum.txt', sep = '\t', header = T, stringsAsFactors = F)
genus2phylum2 <- subset(genus2phylum,!is.na(genus))
row.names(genus2phylum2) <- genus2phylum2$genus

disease <- c('../ibs','../test','../crohn','../CRC','../cirrhosis')
disease2 <- c( 'IBS', 'T2D','CD', 'CRC', 'LC')

spe2genus <- read.table('spe2genus.txt',sep = '\t', header = T, row.names = 2,stringsAsFactors = F,quote="")
phages <- read.table('phage.txt',sep = '\t')
viruses <- read.table('viruses.txt', stringsAsFactors = F)[,1]
i = 1
df_pl <- data.frame() #power-low fit
df_md <- data.frame() #mean degree
df_md_bac <- data.frame() # mean degree of bacteria edges
df_sfi <- data.frame() #scale-free index
df_edges <- data.frame()
df_nodes <- data.frame()
df_mod <- data.frame()
df_path_len <- data.frame()
bac_n_phage <- data.frame() # bacteria that nagatively correlated with phages
bac_n_eu <- data.frame()
dfnames <- c()
mxnames <- c()
g_names <- c()
g_names2 <- c()
complexity <- c()
bac_all <- c()
phage_all <- c()
eu_all <- c()

for(dir in disease){
  diff_file <- paste(disease2[i], 'specific_genus.int', sep = ‘_’)#output of boot_diff.R
  file_case2 <- paste('./SpiecEasi/', paste(disease2[i],'case_spiec_int.txt', sep = '_'), sep = '') # output of SpiecEasi
  file3 <- paste('./dfs/', paste(disease2[i],'df.txt', sep = '_'), sep = '') #output file
  file4 <- paste(dir,'/interaction/bac_vir/nodeAttr.txt', sep = ‘’)#output file

  df <- inter_spiec(diff_file, file_case2)
  #df_con <- inter_spiec(file_con, file_con2)
  
  g <- union(df$OTU1,df$OTU2)
  all <- df
  d <- mapply(comp_d,g,"all")
  d_freq <- count_d(d,'all')
  d_freq$d <- log2(d_freq$d)
  d_freq$Freq <- log2(d_freq$Freq)
  
  #############complexity  
  d2 <- d*log2(d)
  tmp <- sum(d2)
  complexity <- c(complexity,tmp)
  #############complexity
  

  vir_genus <- g[g %in% viruses]
  bac_genus <- setdiff(g,vir_genus)
  bac_all <- c(bac_all, bac_genus)
  df1 <- subset(spe2genus,genus %in% vir_genus)
  phage_genus <- unique(df1[row.names(df1) %in% phages[,2],])
  phage_all <-c(phage_all,phage_genus)
  eu_genus <- setdiff(vir_genus,phage_genus)
  eu_all <- c(eu_all,eu_genus)
  nonphage_genus <- setdiff(g,phage_genus)
  noneu_genus <- setdiff(g,eu_genus)
  
  nonvir_df <- subset(df,!(OTU1 %in% viruses|OTU2 %in% viruses))
  noneu_df <- subset(df,!(OTU1 %in% eu_genus|OTU2 %in% eu_genus))
  nonphage_df <- subset(df,!(OTU1 %in% phage_genus|OTU2 %in% phage_genus))
  df_phage_bac <- subset(df,(OTU1 %in% phage_genus & OTU2 %in% bac_genus) |(OTU1 %in% bac_genus & OTU2 %in% phage_genus))
  df_eu_bac <- subset(df,(OTU1 %in% eu_genus & OTU2 %in% bac_genus) |(OTU1 %in% bac_genus & OTU2 %in% eu_genus))
  # df_phage_vir <- subset(df,(OTU1 %in% phage_genus & OTU2 %in% vir_genus) |(OTU1 %in% vir_genus & OTU2 %in% phage_genus))
  # df_eu_vir <- subset(df,(OTU1 %in% eu_genus & OTU2 %in% vir_genus) |(OTU1 %in% eu_genus & OTU2 %in% phage_genus))
  phage_df <- subset(df,OTU1 %in% phage_genus & OTU2 %in% phage_genus)
  eu_df <- subset(df,OTU1 %in% eu_genus & OTU2 %in% eu_genus)
  ########node number
  
  # df2 <- subset(spe2genus,genus %in% viruses)
  # phage_genus0 <- unique(df2[row.names(df2) %in% phages[,2],]) #all detected phages
  tmp <- data.frame(total = length(g), bacteria = length(bac_genus), phages = length(phage_genus), 'euk-viruses' = length(eu_genus),row.names = disease2[i])
  df_nodes <- rbind(df_nodes,tmp)
  
  ########node number
  ########output node attribute
  g3_names <- c(bac_genus,phage_genus,eu_genus)
  g3 <- c(rep('bacteria', length(bac_genus)), rep('phages', length(phage_genus)), rep('euk-viruses',length(eu_genus)))
  names(g3) <- g3_names
  CC <- mapply(compCC,g,"all")
  #  labels <- unlist(lapply(g,find_class))
  labels <- g3[g]
  cyto_table <- data.frame(allgenus=g, CC, degree = d, labels, stringsAsFactors = F)
  write.table(cyto_table,file4, quote = F,row.names = F, sep = '\t')
  ########output node attribute
  
  #############degree distribution 
  
  d_bac <- mapply(comp_d,bac_genus,"nonvir_df")
  # d_bac <- d[bac_genus]
  d_bac_freq <- count_d(d_bac,'without viral')
  d_bac_freq$d <- log2(d_bac_freq$d)
  d_bac_freq$Freq <- log2(d_bac_freq$Freq)
  

  d_nophage <- mapply(comp_d,nonphage_genus,"nonphage_df")
  d_nophage_freq <- count_d(d_nophage,'without phage')
  d_nophage_freq$d <- log2(d_nophage_freq$d)
  d_nophage_freq$Freq <- log2(d_nophage_freq$Freq)
  
  d_noeuk <- mapply(comp_d,noneu_genus,"noneu_df")
  d_noeuk_freq <- count_d(d_noeuk,'without euk-v')
  d_noeuk_freq$d <- log2(d_noeuk_freq$d)
  d_noeuk_freq$Freq <- log2(d_noeuk_freq$Freq)
  
  tmp <- rbind(d_freq,d_bac_freq,d_nophage_freq,d_noeuk_freq)
  tmp$case <- disease2[i]
  df_pl <- rbind(df_pl, tmp)
########
  sfi_all <- scale_free_index(d)
  sfi_bac <- scale_free_index(d_bac)
  sfi_nophage <- scale_free_index(d_nophage)
  sfi_noeu <- scale_free_index(d_noeuk)
  tmp <- data.frame(all = sfi_all, 'without viral' = sfi_bac, 'without phage' = sfi_nophage, 'without euk' = sfi_noeu)
  row.names(tmp) = disease2[i]
  df_sfi <- rbind(df_sfi, tmp)
  ########degree distrubution
  
  ########mean degree
  tmp <- data.frame(bac_md = mean(d[bac_genus]), phage_md = mean(d[phage_genus]), euk_md = mean(d[eu_genus]), total_md = mean(d))
  row.names(tmp) = disease2[i]
  df_md <- rbind(df_md,tmp)
  
  g_tmp <- union(nonvir_df$OTU1,nonvir_df$OTU2)
  d_nonvir <- mapply(comp_d,g_tmp,'nonvir_df')
  
  # df_phages <-  subset(df,OTU1 %in% phage_genus  |OTU2 %in% phage_genus)
  # df_phage2 <- subset(df,OTU1 %in% phage_genus & OTU2 %in% phage_genus)
  # df_phage3 <- subset(df,(OTU1 %in% phage_genus & OTU2 %in% eu_genus) |(OTU1 %in% eu_genus & OTU2 %in% phage_genus))
  
  g_tmp <- phage_genus[phage_genus %in% union(df_phage_bac$OTU1,df_phage_bac$OTU2)]
  d_phage_bac <- mapply(comp_d,g_tmp,'df_phage_bac')
  g_tmp <- eu_genus[eu_genus %in% union(df_eu_bac$OTU1,df_eu_bac$OTU2)]
  d_eu_bac <- mapply(comp_d,g_tmp,'df_eu_bac')
  tmp <- data.frame(bacteria=mean(d_nonvir),'phage-bacteria'=mean(d_phage_bac),'euk viruses-bacteria'=mean(d_eu_bac), row.names = disease2[i])
  df_md_bac <- rbind(df_md_bac, tmp)
  
  #########mean degree
  ########edge table
  
  edge_num1 <- count_edges(df,bac_genus,bac_genus)
  edge_num2 <- count_edges(df,bac_genus,phage_genus)
  edge_num3 <- count_edges(df,bac_genus,eu_genus)
  edge_num4 <- count_edges(df,phage_genus,phage_genus)
  edge_num5 <- count_edges(df,phage_genus,eu_genus)
  edge_num6 <- count_edges(df,eu_genus,eu_genus)
  
  tmp <- data.frame(matrix(c(edge_num1,edge_num2,edge_num3,edge_num4,edge_num5,edge_num6), nrow = 1))
  tmp <- mutate(tmp,sum = X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12)
  row.names(tmp) = disease2[i]
  df_edges <- rbind(df_edges,tmp)
  
  ########edge table
  
  ######## negatively correlated bacteria
  n_bac_df <- neg_bac(df_phage_bac,bac_genus)
  n_bac_df <- data.frame(n_bac_df,group = disease2[i])
  bac_n_phage <- rbind(bac_n_phage,n_bac_df)
  tmp <- neg_bac(df_eu_bac,bac_genus)
  if (length(tmp) != 0){
    tmp <- data.frame(tmp,group = disease2[i])
    bac_n_eu <- rbind(bac_n_eu,tmp)
  }
  
  ######## negatively correlated bacteria
  ########modularity
  mod1 <- modul(df)
  mod2 <- modul(nonvir_df)
  mod3 <- modul(nonphage_df)
  mod4 <- modul(noneu_df)
  tmp <- data.frame(all = mod1, nonvir= mod2, 'without phage' = mod3, 'without euk-viruses' = mod4)
  row.names(tmp) <- disease2[i]
  df_mod <- rbind(df_mod,tmp)
  ########modularity
  ########path length
  length1 <- path_length(df)
  length2 <- path_length(nonvir_df)
  length3 <- path_length(phage_df)
  length4 <- path_length(eu_df)
  length5 <- path_length(nonphage_df)
  length6 <- path_length(noneu_df)
  tmp <- data.frame(all = length1, bacteria= length2, phage = length3, 'euk-viruses' = length4, 'without phage' = length5, 'without euk-viruses' = length6)
  row.names(tmp) <- disease2[i]
  df_path_len <- rbind(df_path_len,tmp)
  ########path length
  ########aggregate df
  a <- apply(df,1,classify_edges,bac_genus,phage_genus,eu_genus)
  df <- data.frame(df, edge_type = a)
  dfname <- paste(disease2[i], '_df', sep = '')
  dfnames <- c(dfnames,dfname)
  assign(dfname,df)
  write.table(df,file3, sep = '\t', quote = F)
  ########aggregate df
  #   i = i+1
  # }
  ########modularity matrix
  #  df <- cutDupPairs(df0)
  mx <- plot_mx(df)
  mxname <- paste(disease2[i], '_mx', sep = '')
  mxnames <- c(mxnames,mxname)
  assign(mxname,mx)
  g2 <- c(rep('bacteria', length(bac_genus)), rep('phage', length(phage_genus)), rep('euk viruses', length(eu_genus)))
  tmp <- c(bac_genus, phage_genus, eu_genus)
  phylum <- genus2phylum2[tmp,'phylum']
  #  df1 <- subset(genus2phylum,genus %in% bac)
  phylum[is.na(phylum)] <-'Viruses'
  phylum[!(phylum %in% c('Firmicutes','Bacteroidetes','Proteobacteria','Actinobacteria','Viruses'))] <- 'Others'
  g2 <- data.frame('node class' = g2, 'phylum or viruses' = phylum)
  row.names(g2) <- tmp
  g_name <- paste(disease2[i], '_g', sep = '')
  g_names <- c(g_names,g_name)
  assign(g_name, g2)
  
  mx <- plot_mx(nonvir_df)
  mxname <- paste(disease2[i], '_nonvir_mx', sep = '')
  assign(mxname,mx)
  phylum <- genus2phylum2[bac_genus,'phylum']

  phylum[!(phylum %in% c('Firmicutes','Bacteroidetes','Proteobacteria','Actinobacteria'))] <- 'Others'
  g2 <- data.frame('phylum or viruses' = phylum)
  row.names(g2) <- bac_genus
  g_name <- paste(disease2[i], '_nonvir_g', sep = '')
  assign(g_name, g2)

  i = i + 1
}

tmp <- as.data.frame(table(bac_all))
tmp <- tmp[order(tmp$Freq, decreasing = T),]
bac_all <- as.character(tmp$bac_all)
tmp <- as.data.frame(table(phage_all))
tmp <- tmp[order(tmp$Freq, decreasing = T),]
phage_all <- as.character(tmp$phage_all)

tmp <- as.data.frame(table(eu_all))
tmp <- tmp[order(tmp$Freq, decreasing = T),]
eu_all <- as.character(tmp$eu_all)
virus_all <- union(phage_all, eu_all)
genus_all <- c(bac_all,phage_all,eu_all)
write(genus_all, 'genus_all.txt')

names(df_edges) <- c('positive bacteria-bacteria','negative bacteria-bacteria','positive bacteria-phages','negative bacteria-phages','positive bacteria-euk viruses','negative bacteria-euk viruses','positive phages-phages','negative phages-phages','positive phages-euk viruses','negative phages-euk viruses','positive euk viruses-euk viruses','negative euk viruses-euk viruses','sum')
sums2  <- colSums(df_edges)
df_edges <- rbind(df_edges, sums2)
row.names(df_edges)[6] <- 'sum'
df_cplxt <- data.frame(complexity,row.names = disease2)
write.table(df_cplxt,'complexity.txt', quote = F, sep = '\t')
write.table(df_edges,'edge_number.txt', quote = F, sep = '\t')
write.table(df_nodes,'node_number.txt', quote = F, sep = '\t')
write.table(df_sfi,'scale_free_index.txt', quote = F, sep = '\t')
write.table(df_md,'mean_degree.txt', quote = F, sep = '\t')
write.table(df_md_bac,'mean_degree_bac.txt', quote = F, sep = '\t')
write.table(df_mod,'modularity.txt', quote = F, sep = '\t')
write.table(df_path_len,'avg_path_len.txt', quote = F, sep = '\t')

r1 <- df_edges$`negative bacteria-bacteria`/(df_edges$`positive bacteria-bacteria` + df_edges$`negative bacteria-bacteria`)
r2 <- df_edges$`negative bacteria-phages`/(df_edges$`positive bacteria-phages`+ df_edges$`negative bacteria-phages`)
r3 <- df_edges$`negative bacteria-euk viruses`/(df_edges$`positive bacteria-euk viruses` + df_edges$`negative bacteria-euk viruses`)
r4 <- df_edges$`negative phages-phages`/(df_edges$`positive phages-phages` + df_edges$`negative phages-phages`)
r5 <- df_edges$`negative phages-euk viruses`/(df_edges$`positive phages-euk viruses` + df_edges$`negative phages-euk viruses`)
r6 <- df_edges$`negative euk viruses-euk viruses`/(df_edges$`positive euk viruses-euk viruses` + df_edges$`negative euk viruses-euk viruses`)
tmp1 <- data.frame(negative_ratio = r1[-6], cases = row.names(df_edges)[-6],label = 'bacteira-bacteria')
tmp2 <- data.frame(negative_ratio = r2[-6], cases = row.names(df_edges)[-6],label = 'bacteira-phages')
tmp3 <- data.frame(negative_ratio = r3[-6], cases = row.names(df_edges)[-6],label = 'bacteira-euk viruses')
tmp4 <- data.frame(negative_ratio = r4[-6], cases = row.names(df_edges)[-6],label = 'phages-phages')
tmp5 <- data.frame(negative_ratio = r5[-6], cases = row.names(df_edges)[-6],label = 'phages-euk viruses')
tmp6 <- data.frame(negative_ratio = r6[-6], cases = row.names(df_edges)[-6],label = 'euk viruses-euk viruses')

df_neg_ratito <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)
df_neg_ratito$cases <- factor(df_neg_ratito$cases, levels = c('IBS', 'T2D','CD','CRC','LC'))
p <- ggplot(df_neg_ratito,aes(x = cases, y = negative_ratio, fill = label)) + geom_bar(stat = 'identity', position=position_dodge()) + theme_classic()
p + labs(y="Proportions of negative correlations") + theme(axis.text.x = element_text(size = 15, color = 'black'), axis.title.y = element_text(size = 16), axis.title.x = element_blank(), axis.text.y = element_text(size = 15, color = 'black')
   , legend.text = element_text(size = 15), legend.title = element_blank(),legend.key.size = unit(0.4, 'cm'), legend.position = c(0.82,0.87), legend.background = element_rect(linetype = 'solid',size = 0.5, colour ="black")) + scale_fill_tron()
# show scale_fill_tron() color:
#library("scales")
#show_col(pal_tron("legacy")(7))
#result will show in the plot panel

df_pl$case <- factor(df_pl$case, levels = c('IBS','T2D','CD','CRC','LC'))
(p <- ggplot(df_pl,aes(d,Freq,color = labels)) + geom_point(shape = 1) + stat_smooth(method="lm",se = F)
+ facet_wrap(facets= ~ case, scales = "free", nrow = 2) + theme_bw())
p + labs(x = "log2(node degree)", y="log2(node number)") + scale_color_tron() + theme(axis.text = element_text(size = 14, color = 'black'), 
          legend.position = c(0.85,0.25), legend.background = element_rect(linetype = 'solid',size = 0.5, colour ="black"), axis.title = element_text(size = 15), strip.text = element_text(size = 11),  legend.text = element_text(size = 13), legend.title = element_blank())

#######modularity heatmap
ann_colors <- list('node.class' = c(bacteria="#70A6FF", phage = "#FD7C76", 'euk viruses' = "#1BD24A"), 'phylum.or.viruses' = c('Firmicutes' = '#1258DC', 'Bacteroidetes' = '#F94D58', 'Proteobacteria' = '#FFED29', 'Actinobacteria' = '#A51D26','Viruses' = '#8FAFEA', 'Others' = '#98CA32'))
#ann_colors <- list('node.class' = c(bacteria="#70A6FF", phage = "#FD7C76", 'euk viruses' = "#1BD24A"), phylum = c('Actinobacteria' = '#DF5000', 'Bacteroidetes' = '#FAF05D', 'Firmicutes' = '#AAD0F2', 'Proteobacteria' = '#49359B', 'Others' = '#98CA32','Viruses' = '#8FAFEA'))

ann_colors2 <- list( 'phylum.or.viruses' = c('Firmicutes' = '#1258DC', 'Bacteroidetes' = '#F94D58', 'Proteobacteria' = '#FFED29', 'Actinobacteria' = '#A51D26', 'Others' = '#98CA32'))
heatmp <- function(mx, g){
  pheatmap(mx, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), border_color=NA, cluster_cols = F, fontsize = 12,  cluster_rows = F, legend = F, annotation_legend = FALSE, annotation_names_row = F, annotation_names_col = F, fontsize_row = 1, fontsize_col = 1, annotation_row = g, annotation_col = g, annotation_colors = ann_colors, width = 100, height = 100 )#legend_breaks = c(0,0.5), legend_labels = c("this is a     really long item", '0.5')
}
heatmp_nonvir <- function(mx, g){
  pheatmap(mx, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), border_color=NA, cluster_cols = F, fontsize = 12,  cluster_rows = F, legend = F, annotation_legend = FALSE, annotation_names_row = F, annotation_names_col = F, fontsize_row = 1, fontsize_col = 1,  annotation_row = g, annotation_col = g, annotation_colors = ann_colors2) #
}

mxnames
g_names
heatmp(T2D_mx, T2D_g)
heatmp_nonvir(T2D_nonvir_mx, T2D_nonvir_g)

heatmp(CD_mx, CD_g)
heatmp_nonvir(CD_nonvir_mx, CD_nonvir_g)

heatmp(CRC_mx, CRC_g)
heatmp_nonvir(CRC_nonvir_mx, CRC_nonvir_g)

heatmp(LC_mx, LC_g)
heatmp_nonvir(LC_nonvir_mx, LC_nonvir_g)

heatmp(IBS_mx, IBS_g)
heatmp_nonvir(IBS_nonvir_mx, IBS_nonvir_g)
