#!/bin/bash
#SBATCH -n 1
#SBATCH --error=downstream.%J.err
#SBATCH --output=downstream.%J.out
Rscript /data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc/MPRA_summarization_function.R downstream
