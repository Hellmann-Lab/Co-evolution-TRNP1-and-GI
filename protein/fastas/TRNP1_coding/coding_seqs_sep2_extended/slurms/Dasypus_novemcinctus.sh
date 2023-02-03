#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Dasypus_novemcinctus.%J.err
#SBATCH --output=Dasypus_novemcinctus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Dasypus_novemcinctus /data/share/ngs/genomes/UCSC/Dasypus_novemcinctus dasNov3.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
