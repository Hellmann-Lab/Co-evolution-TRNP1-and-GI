#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Piliocolobus_tephrosceles.%J.err
#SBATCH --output=Piliocolobus_tephrosceles.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Piliocolobus_tephrosceles /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Piliocolobus_tephrosceles GCF_002776525.3_ASM277652v3_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
