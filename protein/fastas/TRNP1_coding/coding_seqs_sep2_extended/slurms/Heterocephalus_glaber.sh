#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Heterocephalus_glaber.%J.err
#SBATCH --output=Heterocephalus_glaber.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Heterocephalus_glaber /data/share/ngs/genomes/UCSC/Heterocephalus_glaber hetGla2.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
