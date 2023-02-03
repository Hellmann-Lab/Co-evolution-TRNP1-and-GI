#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Elephantulus_edwardii.%J.err
#SBATCH --output=Elephantulus_edwardii.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Elephantulus_edwardii /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Elephantulus_edwardii GCF_000299155.1_EleEdw1.0_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
