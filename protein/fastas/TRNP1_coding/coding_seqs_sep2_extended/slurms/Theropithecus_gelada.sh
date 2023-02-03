#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Theropithecus_gelada.%J.err
#SBATCH --output=Theropithecus_gelada.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Theropithecus_gelada /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Theropithecus_gelada GCF_003255815.1_Tgel_1.0_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
