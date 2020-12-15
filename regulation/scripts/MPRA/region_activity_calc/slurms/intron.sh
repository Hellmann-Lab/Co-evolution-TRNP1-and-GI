#!/bin/bash
#SBATCH -n 1
#SBATCH --error=intron.%J.err
#SBATCH --output=intron.%J.out
Rscript /data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc/MPRA_summarization_function.R intron
