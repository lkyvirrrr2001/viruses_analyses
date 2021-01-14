for name in `cat $1`  #the first argument is a file list of prefixes of reads file
do
####quality control
    fq1=$name'_1.fastq'
    fq2=$name'_2.fastq'
    prinseq-lite -verbose -fastq $fq1 -ns_max_p 10 -min_qual_mean 25 -out_good ./tmpfiles/$name'_1_good' -out_bad ./tmpfiles/$name'_1_bad'
    prinseq-lite -verbose -fastq $fq2 -ns_max_p 10 -min_qual_mean 25 -out_good ./tmpfiles/$name'_2_good' -out_bad ./tmpfiles/$name'_2_bad'
    manageread extractpair pairfastq=./tmpfiles/$name'_1_good.fastq',./tmpfiles/$name'_2_good.fastq' outprefix=./tmpfiles/$name
    bowtie2 -x /home/hqzhu_pkuhpc/lustre1/mli/5_proj/GRCh38/G38 -1 ./tmpfiles/$name'_1.fastq' -2 ./tmpfiles/$name'_2.fastq' --very-fast --un-conc ./clean_reads/$name'.fastq' -S ./tmpfiles/$name'.sam' --no-unal --no-head -p 20

####assembly
    spades.py --meta -1 ./clean_reads/$name'.1.fastq' -2 ./clean_reads/$name'.2.fastq' -t 10  -o ./contigs/$name
    mv $3/$name/contigs.fasta $3/$name.contig.fasta
    rm -rf $3/$name

####select contigs longer than 1.5kb
    prinseq-lite -fasta ./contigs/$name'.contig.fasta' -min_len 1500 -out_good ./contigs/$name'_1.5k' -out_bad ./contigs/$name'_1.5k_short'

####gene prediction
    gmhmmp -m MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod -d -D ./genes/$name'.gene' -o ./genes/$name'.geneLoc' ../contigs/$name'_1.5k.fasta'

####map reads to contigs to calculate abundance
    bowtie2-build  ./contigs/$name'_1.5k.fasta' ./index/$name
    bowtie2 -x ./index/$name -1 ./clean_reads/$name'.1.fastq' -2 ./clean_reads/$name'.2.fastq' -S ./sam/$name'.sam' --no-unal --no-head -p 16
    
####viruses annotation
    nt2aa ./genes/$name'.gene'
    blastp -query ./genes/$name'.gene.csqa' -db ref_vir_prot -out $name'.blast_ncbiVir' -evalue 1e-5 -max_target_seqs 1 -outfmt 6 -num_threads 16
    hmmsearch --tblout $name'.hmm.result' -E 1e-5 --cpu 16 final_list.hmms ./genes/$name'.gene.csqa'
 
####bacteria annotation
    kaiju -t nodes.dmp -f kaiju/refseq/kaiju_db_refseq.fmi -i ./genes/$name.gene.csqa -E 1e-5 -s 75 -z 20 -o $name.tmp -p
    kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i $name.tmp -o $name.kaiju -r superkingdom,phylum,genus,species -u
  rm -f $name.tmp

####remove bacteria genes
    blastn -query ./genes/$name'.gene' HGG/db/HGG_bacteria.fa -perc_identity 98 -qcov_hsp_perc 90 -out $name.blast_HGG -evalue 1e-10 -max_target_seqs 1 -outfmt 6 -num_threads 20

####KEGG annotation
    sed -i 's/ /_/g' ./genes/$name'.gene.csqa'
    blastall -p blastp -i ./genes/$name'.gene.csqa' -d database/KEGG/KEGG -a 16 -e 1e-5 -v 1 -b 1 -m 8 -o $name'.KEGGe5'

done
