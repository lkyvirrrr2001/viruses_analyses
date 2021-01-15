filelist=$1
for name in `cat $filelist`
do
##keep only the top hit of each query sequence
  sort -k1,1 -k12,12nr -k11,11n $name.blast_ncbiVir | sort -u -k1,1 --merge > $name.blast_ncbiVir2 

#get gene names of aligned sequences of blast output
  awk -F '#' '{OFS="#";print $1,$2}' $name'.blast_ncbiVir2’ > $name'.blast.gene’ 

##remove the last 10 rows and first 3 rows of hmm result files (these are annotation lines).
  max=`sed -n '$=' $name'.hmm.result'`
  let sLine=max-10+1
  sed -i $sLine',$d' $name'.hmm.result'
  sed -i '1,3d' $name'.hmm.result'

##get gene names of aligned sequences of hmm output
  awk -F '#' '{OFS="#";print $1,$2}' $name'.hmm.result' > $name'.hmm.gene'

##get gene names of aligned sequences of blast gut bacteria genomes (HGG) output
  awk -F '\t' '{print $1}' $name.blast_HGG|uniq > $name.HGG.gene

##get union of blast and hmm genes
  sort $name'.hmm.gene' $name'.blast.gene' | uniq > $name'.union.genes'

##get different genes between the union genes and HGG genes
  sort $name.union.genes $name.HGG.gene $name.HGG.gene |uniq -u > $name.viral.genes

done
