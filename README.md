# Co-evolution-TRNP1-and-GI
This repository contains data files and scripts to reproduce the analyses and results presented in our paper.

## Protein

All the relevant scripts listed below can be found within protein/scripts folder.

### Protein-coding sequence collection 

1) `downl_genomes.sh` - to download the genome sequences from ENSEMBL and NCBI

2) `processingL_final.R` - process our own sequence assemblies from targeted re-sequencing

3) `collect_coding_seqs2/run_ccs2_function.sh` - blast wrapper to extract orthologous protein-coding sequences from genomes 

4) `collect_coding_seqs.R` - gather the orthologous TRNP1 protein-coding sequences from all included sources (genomes, resequencing). Intersect with the available trait data. Save sequences and traits for the downstream analyses

#### Protein-coding sequence alignment
align_with_prank.sh - [Multiple Alignments with PRANK](http://wasabiapp.org/software/prank/)

#### Phylogenetic analysis of TRNP1 protein

1) [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
First, run PAML site models as described in the readme in the folder PAML.
`select_sign_sites_PAML_M8.R` - pull out the identified sites under positive selection.

4) [COEVOL](https://github.com/bayesiancook/coevol)
`run_coevol.sh` - wrapper to run Coevol

`summarize_cor_output1.R` - access the estimated correlations and posterior probabilities

`summarize_cor_output2.R` - access the estimated omega of the protein

5) Analysis of NPC proliferation assay
protein/scripts/proliferation_analysis.R - gather proliferation assay data, estimate proliferation rates using logistic regression, infer association with GI using PGLS

## Regulation

1) CRE-orthologue re-sequencing 

2) MPRA-design

3) MPRA count pre-processing
regulation/scripts/MPRA/collect_MPRA_fastas.R - extract sequences from each of the 7 CRE regions, align using [mafft](https://mafft.cbrc.jp/alignment/software/)

regulation/scripts/MPRA/MPRA_analysis.R - filter and summarize CRE activities. Plug into PGLS and compare to brain mass and gyrification

3) Motif identification (JASPAR) and expression via RNA-seq 

regulation/scripts/TFs/motifs_JASPAR2020.R - download PWMs and motif clustering from [JASPAR 2020](http://jaspar.genereg.net/downloads/), transform PWMs for Cluster-Buster

regulation/scripts/TFs/TF_expression_analysis.R - find the expressed transcription factors in our NPCs (from bulk RNA-seq data). Run [Cluster-Buster](http://cagt.bu.edu/page/ClusterBuster_download) on the intron sequences including only the PWMs of the expressed TFs to identify overrepresented motifs

regulation/scripts/TFs/PGLS_motifs.R - investigate binding score assocation with intron CRE activity and GI among the 22 most abundant motifs on the intron sequence using PGLS



Generate figures using figures/figures_forPaper.R 
