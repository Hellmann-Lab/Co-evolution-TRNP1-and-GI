# Co-evolution-TRNP1-and-GI
This repository contains data files and scripts to reproduce the analyses and results presented in our paper.

## Protein

All relevant scripts listed below can be found within the `protein/scripts` folder.

### Protein-coding sequence collection 

1) `collect_coding_seqs2/run_ccs2_function.sh` - blast wrapper to extract orthologous protein-coding sequences from genomes. If you want to do the blast search using our scripts, you will need to download genomes listed in table `protein/data/TRNP1_source_genomes.csv`

2) `processingL_final.R` - process our own sequence assemblies from targeted re-sequencing

3) `collect_coding_seqs.R` - gather the orthologous TRNP1 protein-coding sequences from all included sources. Intersect with the available trait data. Save sequences and traits for the downstream analyses


### Evolutionary analysis of TRNP1 coding-sequence

#### [Multiple Alignments with PRANK](http://wasabiapp.org/software/prank/) (v150803)
`align_with_prank.sh` -  protein-coding sequence alignment

#### [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (v4.8)
First, run PAML site models as described in the readme from folder `PAML`.
`select_sign_sites_PAML_M8.R` - pull out the identified sites under positive selection

#### [COEVOL](https://github.com/bayesiancook/coevol) (v1.4)

1) `run_coevol.sh` - wrapper to run Coevol

2) `summarize_coevol_output1.R` - access the estimated correlations and posterior probabilities

3) `summarize_coevol_output2.R` - access the estimated omega of the protein

### Analysis of NPC proliferation assay
`proliferation_analysis.R` - gather proliferation assay data, estimate proliferation rates using logistic regression, infer association with GI using PGLS


## Regulation

All relevant scripts listed below can be found within the `regulation/scripts` folder.


### MPRA assay 

1) `MPRA/MPRA_sequences.R` - identify and collect orthologous TRNP1 CRE sequences across mammals from our sequenced data as well as published genomes

2) `MPRA/MPRA_oligolib_construction.R` - MPRA design -  using a sliding window, construct enhancer tiles based on the orthologous CRE sequences from the previous script to test within the MPRA assay

3) `MPRA/preprocessing`  MPRA count pre-processing. Extract reporter gene expression counts for each included enhancer tile. This folder contains a README with further details 

4) `MPRA/collect_MPRA_fastas.R` - separate and save the relevant sequences from each of the 7 CRE regions, align using [MAFFT](https://mafft.cbrc.jp/alignment/software/) (v7.407)

5) `MPRA/MPRA_analysis.R` - filter and summarize CRE activities. Plug into PGLS and compare to brain mass and gyrification

6) `MPRA/combine_dnds_intron.R` - combine TRNP1 protein evolution rates inferred using Coevol with the intron activity across catharrines within the same model

### Transcription factor analysis

1) `TFs/download_motifs_JASPAR2020.R` - download PWMs and motif clustering from [JASPAR 2020](http://jaspar.genereg.net/downloads/), transform PWMs for Cluster-Buster

2) `TFs/MPRNAseq_NPC.yaml` - [zUMIs](https://github.com/sdparekh/zUMIs) (v2.5.4) yaml file for mapping RNA-seq reads from NPCs. Input raw data for this processing can be accessed under E-MTAB-9951

3) `TFs/TF_expression_analysis.R` - find the expressed transcription factors in our NPCs (from bulk RNA-seq data). Run [Cluster-Buster](http://cagt.bu.edu/page/ClusterBuster_download) (Jun 13 2019) on the intron sequences including only the PWMs of the expressed TFs to identify overrepresented motifs

4) `TFs/PGLS_motifs.R` - investigate binding score assocation with intron CRE activity and GI among the 22 most abundant motifs on the intron sequence using PGLS


Tree construction: `regulation/scripts/MPRA/tree_construction.R`

Throughout the workflow, we are using job scheduling system [slurm](https://github.com/mattthias/slurm) (v0.4.3).
