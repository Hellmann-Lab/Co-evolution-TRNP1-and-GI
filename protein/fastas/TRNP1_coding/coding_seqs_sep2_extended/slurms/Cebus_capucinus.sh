#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Cebus_capucinus.%J.err
#SBATCH --output=Cebus_capucinus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Cebus_capucinus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Cebus_capucinus Cebus_capucinus.Cebus_imitator-1.0.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
