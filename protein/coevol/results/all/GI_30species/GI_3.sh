#!/bin/bash
#SBATCH -n 1
#SBATCH --error=GI_3.%J.err
#SBATCH --output=GI_3.%J.out
#SBATCH --time=4-00:00:00
srun coevol -d /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol_longer.best.nuc.phy         -t /data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_30sp_forCoevol.txt         -c /data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data//GI_exon1_30species.lht -dsom exon1_GI_dsomggc3 
