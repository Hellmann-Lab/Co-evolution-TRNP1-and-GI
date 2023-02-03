#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Chlorocebus_sabeus.%J.err
#SBATCH --output=Chlorocebus_sabeus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Chlorocebus_sabeus /data/share/ngs/genomes/UCSC/Chlorocebus_sabeus chlSab2.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
