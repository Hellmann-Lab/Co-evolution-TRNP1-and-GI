#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Delphinapterus_leucas.%J.err
#SBATCH --output=Delphinapterus_leucas.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Delphinapterus_leucas /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Delphinapterus_leucas GCF_002288925.2_ASM228892v3_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
