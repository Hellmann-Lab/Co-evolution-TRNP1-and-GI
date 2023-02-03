#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Rattus_norvegicus.%J.err
#SBATCH --output=Rattus_norvegicus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Rattus_norvegicus /data/share/ngs/genomes/UCSC/Rattus_norvegicus rn6.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
