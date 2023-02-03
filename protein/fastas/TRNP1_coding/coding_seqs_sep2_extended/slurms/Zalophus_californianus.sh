#!/bin/bash
#SBATCH -n 1
#SBATCH --error=Zalophus_californianus.%J.err
#SBATCH --output=Zalophus_californianus.%J.out
#SBATCH --nodelist=gorilla1
#SBATCH -p high
Rscript /data/share/htp/TRNP1/paper_data/protein/scripts/collect_coding_seqs2_extended/collect_coding_seqs2_function.R Zalophus_californianus /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Zalophus_californianus GCF_009762305.2_mZalCal1.pri.v2_genomic.fna.gz /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/coding_seqs_sep2_extended
