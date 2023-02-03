#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Orcinus_orca.%J.err
#SBATCH --output=Orcinus_orca.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Orcinus_orca /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Orcinus_orca GCF_000331955.2_Oorc_1.1_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
