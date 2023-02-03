#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Callithrix_jacchus.%J.err
#SBATCH --output=Callithrix_jacchus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Callithrix_jacchus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Callithrix_jacchus Callithrix_jacchus.ASM275486v1.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
