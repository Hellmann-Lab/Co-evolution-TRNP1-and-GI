#!/bin/bash
#SBATCH -n 1
#SBATCH --error=body_mass_brain_mass_GI_5.%J.err
#SBATCH --output=body_mass_brain_mass_GI_5.%J.out
srun coevol -d /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol.best.nuc.phy         -t /data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt         -c /data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data//body_mass_brain_mass_GI_exon1_31species.lht -dsom exon1_body_mass_brain_mass_GI_dsomggc5 
