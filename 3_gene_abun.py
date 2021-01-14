#!usr/bin/python
from multiprocessing import Pool
import sys
import numpy as np
import pandas as pd

def sam_abun(abun_file):
    df_abun = pd.read_csv(abun_file,header=None,sep = '\t',usecols=[0,2])
    grp = df_abun.groupby([2],as_index=False).size().reset_index(name='num')
    grp.columns = ['contigs','num']
    ##adjust contig length of contig abundance
    contigs_abun = grp["contigs"].tolist()
    tmp = [x.split('_')[3] for x in contigs_abun]
    contig_len2 = map(float, tmp)
    adj_num = grp["num"]*1000/contig_len2
    grp['adj_num'] = pd.Series(adj_num, index=grp.index)
    return(grp)

def gene_count(gene_file, grp):
    genes = []
    contigs = []
    for line in open(gene_file):
        if line[0] == '>':
            gene0 = line.strip('>').split('#')[0]
            genes.append(gene0)
            contig0 = line.strip('\n').split('#')[1]
            contigs.append(contig0)

    i = 0
    gene_abun = []
    for gene in genes:
        contig = contigs[i]
        num = grp.loc[grp.loc[:,'contigs']==contig,'adj_num']
        if len(num) == 0:
            num = 0 
        num2 = float(num)
        gene_abun.append(num2)
        i = i+1

    d = {'gene': genes, 'abun': gene_abun}
    df_abun = pd.DataFrame.from_dict(d)
    return(df_abun)

if __name__=='__main__':
    filelist = open(sys.argv[1]).readlines() #arg1 is a filelist
    SRRs = [line.strip("\n") for line in filelist]
    ofile = [sys.argv[3] + '/' + x + '.abun' for x in SRRs] # arg3 is output file dir
    file_gene = [sys.argv[2] + '/' + x + ‘.gene' for x in SRRs]#arg2 is the gene sequence file dir
    file_abun = ['./sam/' + x + '.sam' for x in SRRs]
    pool = Pool(16)
    abun_grp_list = pool.map(sam_abun,file_abun)
    
    i = 0
    while i< len(SRRs):
        a = SRRs[i]
        df1 = gene_count(file_gene[i], abun_grp_list[i])
        df1.to_csv(ofile[i], sep = '\t', index = False)
        i = i+1
