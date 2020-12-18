#!/bin/bash
#SBATCH --workdir=/data/share/htp/TRNP1/MPRA_full_300919/MPRA_data/deML/demult/files
#SBATCH --error=MPRA_II_filt2_%J.err
#SBATCH --output=MPRA_II_filt2_%J.out
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Variables #
maindir="/data/share/htp/TRNP1/MPRA_full_300919/MPRA_data/"

# Matthias Pipeline for barcode filtering based on Phred quality and the constant region (beginning of GFP) #
###runs as bash pipeline.sh <list of input files to process> <Phred quality to filter> <study name> 
bash ${maindir}/scripts/III_quality.filtering_counting_pipeline_GFP2.sh \
demult*r1.fq.gz \
10 \
full_MPRA

# Reads through filtering
bash ${maindir}/scripts/IV_reads.through.filtering_all.samples2.sh \
filtered.demult*/report.demult*fq.gz2.txt

cp ${maindir}/deML/demult/files/filtered.demult*/full_MPRA.demult*fq.gz2.txt ${maindir}/deML/demult/count_summaries/counts

