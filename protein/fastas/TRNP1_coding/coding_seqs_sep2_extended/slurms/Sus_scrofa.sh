#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Sus_scrofa.%J.err
#SBATCH --output=Sus_scrofa.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Sus_scrofa /data/share/ngs/genomes/UCSC/Sus_scrofa susScr3.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
