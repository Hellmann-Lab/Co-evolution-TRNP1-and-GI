#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Felis_catus.%J.err
#SBATCH --output=Felis_catus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Felis_catus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Felis_catus Felis_catus.Felis_catus_9.0.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
