#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Papio_hamadryas.%J.err
#SBATCH --output=Papio_hamadryas.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Papio_hamadryas /data/share/ngs/genomes/UCSC/Papio_hamadryas papHam1.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
