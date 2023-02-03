#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Pan_troglodytes.%J.err
#SBATCH --output=Pan_troglodytes.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Pan_troglodytes /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Pan_troglodytes GCF_002880755.1_Clint_PTRv2_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
