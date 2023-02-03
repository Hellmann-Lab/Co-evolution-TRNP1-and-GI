#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Macaca_fascicularis.%J.err
#SBATCH --output=Macaca_fascicularis.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Macaca_fascicularis /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Macaca_fascicularis Macaca_fascicularis.Macaca_fascicularis_5.0.dna.toplevel.fa.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
