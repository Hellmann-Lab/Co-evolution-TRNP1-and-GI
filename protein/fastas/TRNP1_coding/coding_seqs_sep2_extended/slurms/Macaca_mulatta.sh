#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Macaca_mulatta.%J.err
#SBATCH --output=Macaca_mulatta.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Macaca_mulatta /data/share/ngs/genomes/UCSC/Macaca_mulatta/RheMac8 rheMac8.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
