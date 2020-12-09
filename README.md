# Co-evolution-TRNP1-and-GI
This repository contains data files and scripts to reproduce the analyses and results presented in our paper.

##Protein

1) Protein-coding sequence collection 
scripts/downl_genomes.sh - to download the genome sequences used
scripts/processingL_final.R - process our own sequence assemblies from targeted re-sequencing
scripts/collect_coding_seqs2 - blast wrapper to extract orthologous protein sequences from genomes 
scripts/collect_coding_seqs.R

2) [Multiple Alignments with PRANK](http://wasabiapp.org/software/prank/)
align_with_prank.sh

3) [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
codeml script.ctl; readme in the folder
select_sign_sites_PAML_M8.R

4) [COEVOL](https://github.com/bayesiancook/coevol)
scripts/run_coevol.sh - wrapper to run coevol
summarize_cor_output.R 

5) Analysis of NPC proliferation assay
proliferation_analysis.R

##Regulation

1) CRE-orthologue re-sequencing 

2) MPRA-design

3) MPRA count pre-processing

3) Motif identification (JASPAR) and expression via RNA-seq 

-zumis
motifs_JASPAR2020.R 
TF_expression_analysis.R + run_clusterbuster.sh + cbust_subset_per_species.sh
PGLS_motifs.R 
