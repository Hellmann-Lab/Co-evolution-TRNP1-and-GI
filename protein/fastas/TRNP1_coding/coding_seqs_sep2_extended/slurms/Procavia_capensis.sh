#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Procavia_capensis.%J.err
#SBATCH --output=Procavia_capensis.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Procavia_capensis /data/share/ngs/genomes/UCSC/Procavia_capensis proCap1.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
