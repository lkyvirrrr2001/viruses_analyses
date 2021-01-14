#!/bin/bash
#$1 is abun matrix $2 is the output file form sparcc_rval.sh
prefix=${2%.*}
python3 /home/hqzhu_pkuhpc/lustre1/tools/SparCC3-master/MakeBootstraps.py $1 -n 100 -t $prefix'_boot_#.txt'
for ((i=0;i<=99;i++))
do
python3 /home/hqzhu_pkuhpc/lustre1/tools/SparCC3-master/SparCC.py $prefix'_boot_'$i'.txt' -i 5 --cor_file=$prefix'_cor_'$i'.txt'
done
python3 /home/hqzhu_pkuhpc/lustre1/tools/SparCC3-master/PseudoPvals.py $2 $prefix'_cor_#.txt' 100 -o $prefix'_two_sided.txt' -t 'two_sided'
