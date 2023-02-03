#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Otolemur_garnettii.%J.err
#SBATCH --output=Otolemur_garnettii.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Otolemur_garnettii /data/share/ngs/genomes/UCSC/Otolemur_garnettii otoGar3.fa /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
