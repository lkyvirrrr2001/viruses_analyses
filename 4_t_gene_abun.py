###adjust gene abundance by gene number on a contain
import sys
import numpy as np
import pandas as pd
from collections import Counter
def adjust_abun(abun_file,gene_file):
    genes_ctg = {}
    genes = []
    contigs = []
    for line in open(gene_file):
        if line[0] == '>':
            gene0 = line.strip('>').split('#')[0]
            genes.append(gene0)
            contig0 = line.strip('\n').split('#')[1]
            contigs.append(contig0)
            genes_ctg[gene0] = contig0
            
    ctg_num = dict(Counter(contigs))
    num = [ctg_num[genes_ctg[x]] for x in genes]
    gene_num = dict(zip(genes,num))

    i = 0
    genes2 = []
    abuns = []
    for line in open(abun_file):
        if i != 0:
            tmp = line.strip('\n').split('\t')
            gene = tmp[1]
            genes2.append(gene)
            abun = float(tmp[0])
            num2 = gene_num[gene]
            abun2 = abun/num2
            abuns.append(abun2)
        i = i+1
    d = {'gene': genes2, 'abun': abuns}
    df_abun = pd.DataFrame.from_dict(d)
    return(df_abun)

if __name__=='__main__':
    filelist = open(sys.argv[1]).readlines() #arg1 is the prefixes file list
    SRRs = [line.strip("\n") for line in filelist]
    ofile = [sys.argv[4] + '/' + x + '.abun' for x in SRRs]# arg4 is output file dir
    file_gene = [sys.argv[2] + '/' + x + '.gene' for x in SRRs]#arg2 is the gene sequence file dir
    file_abun = [sys.argv[3] + '/' + x + '.abun' for x in SRRs] #arg3 is the original abundance file output by gene_abun.py dir
    i = 0
    while i< len(SRRs):
        a = SRRs[i]
        df1 = adjust_abun(file_abun[i],file_gene[i])
        df1.to_csv(ofile[i], sep = '\t', index = False)
        i = i+1
