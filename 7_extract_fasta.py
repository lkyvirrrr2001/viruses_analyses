#!usr/bin/python
from extract_fasta import extract_fasta
import sys
def extract_fasta(infile,headers, ofile):
    file3 = open(infile)
    fa_dict = {}
    i = 1
    header0 = []
    while True:
        line = file3.readline()
        if len(line) == 0:
            fa_dict[key] = seq
            break
        if i == 1:
            key = line.strip('\n').lstrip('>') 
            header0.append(key)
            seq = []
        else:
            if line[0] == '>':
                fa_dict[key] = seq
                key = line.strip('\n').lstrip('>')#also modify this line
                header0.append(key)
                seq = []
            else:
                seq.append(line)
        i = i+1  
    file3.close()

    file1 = open(ofile,'w')
    for key in header0:
        if key not in headers:
            seqs = fa_dict[key]
            key = '>' + key
            file1.write('%s\n' %(key))
            for seq in seqs:
                file1.write('%s' %(seq))
    file1.close()

filelist = sys.argv[1]
infile_dir = sys.argv[2] #dir of fasta file extracted from
header_dir = sys.argv[3] #dir of header files of genes to be extracted( header is the sequence name that followed after ‘>’ in fasta file. selected gene fasta header)
ofile_dir = sys.argv[4]  #dir of output file

for line in open(filelist):
    a = line.strip('\n')
    infile = infile_dir + '/' + a + '.gene'
    outfile = ofile_dir + '/' + a + '.seqs'
    header_file = header_dir + '/' + a + '.union.genes'

    header = []
    for line2 in open(header_file):
        header.append(line2.rstrip('\n'))

    extract_fasta(infile,header,outfile)


