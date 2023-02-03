#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Tursiops_truncatus.%J.err
#SBATCH --output=Tursiops_truncatus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Tursiops_truncatus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Tursiops_truncatus Tur_tru_Illumina_dip_v1.fasta /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
