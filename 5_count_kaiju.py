#!usr/bin/python
#usage: python count_kaiju.py filelist dir_anno dir_abun outfile
import sys
import numpy as np
import pandas as pd


if __name__=='__main__':
    filelist = open(sys.argv[1]).readlines() #arg1 is a filelist
    SRRs = [line.strip("\n") for line in filelist]
    ofile = open(sys.argv[4],'w')  #arg4 is the output file

    file_anno = [sys.argv[2] + '/' + x + '.kaiju' for x in SRRs]  #arg2 is the kaiju output dir 
    file_abun = [sys.argv[3] + "/" + x + '.abun' for x in SRRs] # agr3 is the gene abun file dir
    abun_dict = {}
    i=0
    cogs = []
    while i< len(SRRs):
        a = SRRs[i]
        ofile.write('\t%s' %(a))
        gene_anno = pd.read_csv(file_anno[i],sep='\t',header=None, usecols=[1,3])
        yp = gene_anno[3].tolist()
        tmp = gene_anno[1].tolist()
        genes = [x.split('#')[0] for x in tmp]
        

        gene_abun = pd.read_csv(file_abun[i],sep='\t')
        genes2 = gene_abun['gene'].tolist()
        abun = gene_abun['abun'].tolist()
        d2 = dict(zip(genes2,abun))
        abun2 = [d2[x] for x in genes]
        d3 = {'yp':yp, 'abun':abun2}
        df3 = pd.DataFrame(data = d3)
        df4 = df3.groupby('yp',as_index=False).sum()
        virus_dict = df4.set_index('yp').to_dict()
        
        
        cogs = cogs + yp
        abun_dict[a] = virus_dict.pop('abun')
        i = i+1
    ofile.write('\n')
    
    cogs = list(set(cogs))
    for cog in cogs:
        ofile.write("%s" %(cog))
        i = 0
        for sample2 in SRRs:
            l = len(SRRs)
            dict2 = abun_dict[sample2]
            abun = dict2.get(cog)
            if not abun:
                abun = 0
            ofile.write("\t%s" %(abun))
            i = i + 1
            if i == l:
                ofile.write('\n')

    ofile.close()
