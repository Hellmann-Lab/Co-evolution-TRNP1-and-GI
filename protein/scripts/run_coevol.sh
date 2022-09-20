#!/bin/bash

declare -a arr=("brain_mass" "body_mass" "body_mass_brain_mass_GI" "GI") 
declare -a arr_run=(1 2 3)
declare -a arr_sp=(30)

for sp in "${arr_sp[@]}"
do


#define file locations
prank_alignment_loc=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_"$sp"species_forCoevol_longer.best.nuc.phy
tree_loc=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_"$sp"sp_forCoevol.txt
pheno_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data/
coevol_output_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/results/all



for i in "${arr[@]}"
do
mkdir $coevol_output_loc/"$i"_"$sp"species

for run in "${arr_run[@]}"

   do
        # HEADER
        echo '#!/bin/bash' > $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh  
        echo '#SBATCH -n 1' >> $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh
        echo '#SBATCH --error='$i'_'$run'.%J.err' >> $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh 
        echo '#SBATCH --output='$i'_'$run'.%J.out' >> $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh
        echo '#SBATCH --time=4-00:00:00' >> $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh
        
        # COEVOL COMMAND
        echo "srun coevol -d $prank_alignment_loc \
        -t $tree_loc \
        -c $pheno_loc/"$i"_exon1_"$sp"species.lht -dsom exon1_"$i"_dsomggc"$run" " >> $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh 

        # SUBMIT THE TASK 
        cd $coevol_output_loc/"$i"_"$sp"species/
        sbatch $coevol_output_loc/"$i"_"$sp"species/"$i"_"$run".sh

done
done
done


