#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Lipotes_vexillifer.%J.err
#SBATCH --output=Lipotes_vexillifer.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Lipotes_vexillifer /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Lipotes_vexillifer GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
