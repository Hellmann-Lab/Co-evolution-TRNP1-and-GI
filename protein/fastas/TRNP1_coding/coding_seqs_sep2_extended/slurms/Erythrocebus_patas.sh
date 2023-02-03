#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Erythrocebus_patas.%J.err
#SBATCH --output=Erythrocebus_patas.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Erythrocebus_patas /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Erythrocebus_patas GCA_004027335.1_EryPat_v1_BIUU_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
