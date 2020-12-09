#!/bin/bash
#SBATCH -n 1
#SBATCH --error=brain_mass_2.%J.err
#SBATCH --output=brain_mass_2.%J.out
srun coevol -d /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol.best.nuc.phy         -t /data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt         -c /data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data//brain_mass_exon1_31species.lht -dsom exon1_brain_mass_dsomggc2 
