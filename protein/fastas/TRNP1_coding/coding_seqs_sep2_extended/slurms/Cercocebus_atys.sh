#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Cercocebus_atys.%J.err
#SBATCH --output=Cercocebus_atys.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Cercocebus_atys /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Cercocebus_atys Cercocebus_atys.Caty_1.0.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
