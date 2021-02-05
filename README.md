# Script description

The scripts in this depository include the pipeline of sequence processing and the major analyses in the paper. Several scripts should input a file list, with each prefix of file name in a row. ‘dir’ in the following descriptions stands for directory.

1_pre-processing-pipeline.sh: This script includes commands of quality control of raw reads, assembly, selection of contigs, gene prediction, reads mapping, virus and bacteria annotation and KEGG annotation. Before running blast or bowtie2, an index should be built from the template database. Input the file list.

2_functional_annotation.sh: This script is used for functional annotation for detected viral genes.

3_gene_abun.py: This script is used for gene abundance calculation based on the sam file generated in the mapping step of script 1_pre-processing-pipeline.sh. This abundance is equal to the contig abundance where the gene locates. This abundance is used for functional abundance matrix calculation (ACLAME, VGFunC). Usage: `python 3_gene_abun.py file_list gene_sequence.fa output_dir`.

4_t_gene_abun.py: This script is used to adjust contig abundance by the gene number of each contig. The reason to adjust gene number is that species abundance should not relate to the gene number in its genome or contig. Usage: `python 4_t_gene_abun.py file_list gene_sequence.fa 3_gene_abun_output_dir output_dir`.

5_count_kaiju.py: This script is used to calculate bacteria abundance based on the output of kaiju and gene abundance output by 4_t_gene_abun.py. Usage: python 5_count_kaiju.py file_list.txt kaiju_output_dir gene_abundance_dir output_file_name
6_viral_genes.sh: This script is used to select viral genes from the blast and hmmsearch output files. Input the file list.

7_extract_fasta.py: This script is used to extract viral sequences from the gene sequence fasta files. Usage: `python extract_fasta.py file_list.txt gene_seq_file.fa viral_gene_header_file_dir output_dir`.

8_count_viruses.py: This script is used to calculate virus abundance in species level based on the RefSeq annotation and gene abundance output by 4_t_gene_abun.py. Usage: `python 8_count_viruses.py file_list.txt blast_output_dir gene_abundance_dir output_file_name`.

9_sparcc_rval.sh: Comput SparCC correlation between each pair of items in the input abundance matrix. Usage: `bash 9_sparcc_rval.sh abundance_matrix.txt output`. The input abundance matrix should contain both bacteria and viruses.

10_sparcc_pval.sh: Calculate p-values of the correlation calculated by 9_sparcc_rval.sh through bootstrap. Usage: `bash 10_sparcc_pval.sh abundance_matrix.txt 9_sparcc_rval _output`. This command will generate a output file with suffix of ‘_two_sided.txt’ when completed.

11_matrix.r: This script is used to select interaction that satisfy the threshold requirements for r and p. Usage: `Rscript matrix.r 9_sparcc_rval_output 10_sparcc_pval_output output_file`.

12_boot_diff.R: This script is used to generate disease-specific network. The network of case and control generated in the last step should input into this script.

13_spieceasi.r: This script is used to construct network by SpiecEasi. The input is the same abundance matrix file as the input of 9_sparcc_rval.sh.

14_network_analysis.R: This script includes the major code for network analyses.

# Data description
disease-specific-networks.xlsx: disease-specific Networks that can be input into Cytoscape or circos for visualization.

viral_gene_annotation.tar.gz: Annotations for the detected viral genes in this study including taxonomy, KEGG, ACLAME, and VGFunC category.

family_annotation.txt: The VGFunC functional annotation of ACLAME viral protein families.

function_category.txt: The mapping file of functional annotation to VGFunC. categories.

cate-function.txt: The functional annotation of each VGFunC category.

metabo-profiles.xlsx: The host motabolome profiles.
