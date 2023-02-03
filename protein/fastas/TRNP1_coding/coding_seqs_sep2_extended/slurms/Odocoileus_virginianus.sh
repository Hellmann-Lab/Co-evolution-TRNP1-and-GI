#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Odocoileus_virginianus.%J.err
#SBATCH --output=Odocoileus_virginianus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Odocoileus_virginianus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Odocoileus_virginianus GCF_002102435.1_Ovir.te_1.0_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
