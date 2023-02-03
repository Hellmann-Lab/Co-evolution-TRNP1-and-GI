#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Mus_caroli.%J.err
#SBATCH --output=Mus_caroli.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Mus_caroli /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Mus_caroli Mus_caroli.CAROLI_EIJ_v1.1.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
