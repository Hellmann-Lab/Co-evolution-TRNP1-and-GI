#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Microcebus_murinus.%J.err
#SBATCH --output=Microcebus_murinus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Microcebus_murinus /data/share/ngs/genomes/UCSC/Microcebus_murinus/micMur2 micMur2.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
