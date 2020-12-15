#!/bin/bash
#SBATCH -n 1
#SBATCH --error=upstream1.%J.err
#SBATCH --output=upstream1.%J.out
Rscript /data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc/MPRA_summarization_function.R upstream1
