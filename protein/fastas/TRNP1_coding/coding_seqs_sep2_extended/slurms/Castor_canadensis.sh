#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Castor_canadensis.%J.err
#SBATCH --output=Castor_canadensis.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Castor_canadensis /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Castor_canadensis GCF_001984765.1_C.can_genome_v1.0_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
