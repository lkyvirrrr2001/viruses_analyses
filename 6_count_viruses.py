#!usr/bin/python
import sys
import numpy as np
import pandas as pd

YP_virus = {}
for line in open('viral.1.protein.faa.pre'):
    if line[0] == '>':
        a = line.strip('\n').split('[')[1].strip(']')
        YP = line.strip('\n').strip('>').split(' ')[0]
        YP_virus[YP] = a

if __name__=='__main__':
    filelist = open(sys.argv[1]).readlines() #arg1 is a filelist
    SRRs = [line.strip("\n") for line in filelist]
    ofile = open(sys.argv[4],'w')  #arg4 is the output file
fasta 
    file_anno = [sys.argv[2] + '/' + x + '.blast_ncbiVir2' for x in SRRs]  #arg2 is the blast ncbi result 
    file_abun = [sys.argv[3] + "/" + x + '.abun' for x in SRRs] # agr3 is the abun(.sam) file dir
    gene_ref = YP_virus.keys()
    abun_dict = {}
    i=0
    cogs = []
    while i< len(SRRs):
        a = SRRs[i]
        ofile.write('\t%s' %(a))
        gene_anno = pd.read_csv(file_anno[i],sep='\t',header=None,usecols=[0,1])
        yp = gene_anno[1].tolist()
        tmp = gene_anno[0].tolist()
        genes = [x.split('#')[0] for x in tmp]

        gene_abun = pd.read_csv(file_abun[i],sep='\t')
        genes2 = gene_abun['gene'].tolist()
        abun = gene_abun['abun'].tolist()
        d2 = dict(zip(genes2,abun))
        abun2 = [d2[x] for x in genes]
        d3 = {'yp':yp, 'abun':abun2}
        df3 = pd.DataFrame(data = d3)

        yp_tmp = list(set(yp).intersection(set(gene_ref)))
        yp2 = [x for x in yp if x in yp_tmp]
        key_virus = [YP_virus[x] for x in yp2]
        cogs = cogs + key_virus
        df4 = df3[df3['yp'].isin(yp2)]
        tmp_d = {'virus':key_virus, 'abun':df4['abun'].tolist()}
        virus_df = pd.DataFrame.from_dict(tmp_d)
        virus_df2 = virus_df.groupby('virus', as_index=False).sum()
        virus_dict = virus_df2.set_index('virus').to_dict()
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
    
