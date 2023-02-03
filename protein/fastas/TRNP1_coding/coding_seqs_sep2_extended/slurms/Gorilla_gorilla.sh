#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Gorilla_gorilla.%J.err
#SBATCH --output=Gorilla_gorilla.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Gorilla_gorilla /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Gorilla_gorilla GCA_900006655.3_Susie3_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
