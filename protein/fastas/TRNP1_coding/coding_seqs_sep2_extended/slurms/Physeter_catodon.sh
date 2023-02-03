#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Physeter_catodon.%J.err
#SBATCH --output=Physeter_catodon.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Physeter_catodon /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Physeter_catodon GCF_002837175.2_ASM283717v2_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
