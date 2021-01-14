{\rtf1\ansi\ansicpg936\cocoartf1404\cocoasubrtf470
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 for name in `cat $1`  #the first argument is a file list of prefixes of reads file\
do\
    sed -i 's/ /_/g' $name'.seqs.csqa'  #input viral translated gene sequence file\
    /lustre1/hqzhu_pkuhpc/xqwang/blast/blast-2.2.24/bin/blastall -p blastp -i $dir/$name'.seqs.csqa' -d /home/hqzhu_pkuhpc/lustre1/mli/5_proj/aclame/aclame_proteins -a 20 -e 1e-5 -v 1 -b 1 -m 8 -o $name'.ACLe5'\
\
done}