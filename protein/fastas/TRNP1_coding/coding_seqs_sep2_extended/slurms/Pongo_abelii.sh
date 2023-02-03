#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Pongo_abelii.%J.err
#SBATCH --output=Pongo_abelii.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Pongo_abelii /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Pongo_abelii GCF_002880775.1_Susie_PABv2_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
