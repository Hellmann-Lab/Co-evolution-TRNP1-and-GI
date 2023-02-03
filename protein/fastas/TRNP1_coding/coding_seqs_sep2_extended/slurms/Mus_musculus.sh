#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Mus_musculus.%J.err
#SBATCH --output=Mus_musculus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Mus_musculus /data/share/ngs/genomes/UCSC/Mus_musculus mm10.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
