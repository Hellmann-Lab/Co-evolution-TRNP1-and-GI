#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Phoca_vitulina.%J.err
#SBATCH --output=Phoca_vitulina.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Phoca_vitulina /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Phoca_vitulina GCF_004348235.1_GSC_HSeal_1.0_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
