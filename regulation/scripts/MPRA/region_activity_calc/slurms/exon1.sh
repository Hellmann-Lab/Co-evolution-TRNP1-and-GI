#!/bin/bash
#SBATCH -n 1
#SBATCH --error=exon1.%J.err
#SBATCH --output=exon1.%J.out
Rscript /data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc/MPRA_summarization_function.R exon1
