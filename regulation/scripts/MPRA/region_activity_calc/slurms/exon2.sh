#!/bin/bash
#SBATCH -n 1
#SBATCH --error=exon2.%J.err
#SBATCH --output=exon2.%J.out
Rscript /data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc/MPRA_summarization_function.R exon2
