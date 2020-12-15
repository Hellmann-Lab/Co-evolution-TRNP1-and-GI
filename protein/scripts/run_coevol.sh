#!/bin/bash

#define file locations-- NEED TO ADJUST (INES?)
#coevol_input_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/scripts
prank_alignment_loc=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol.best.nuc.phy
tree_loc=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt
pheno_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data/
coevol_output_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/results/all


declare -a arr=("body_mass_brain_mass_GI" "GI" "brain_mass" "body_mass")
declare -a arr_run=(1 2 3)



for i in "${arr[@]}"
do
mkdir $coevol_output_loc/"$i"_31species
done


for i in "${arr[@]}"
do

for run in "${arr_run[@]}"

   do
        # HEADER
        echo '#!/bin/bash' > $coevol_output_loc/"$i"_31species/"$i"_"$run".sh  
        echo '#SBATCH -n 1' >> $coevol_output_loc/"$i"_31species/"$i"_"$run".sh
        echo '#SBATCH --error='$i'_'$run'.%J.err' >> $coevol_output_loc/"$i"_31species/"$i"_"$run".sh 
        echo '#SBATCH --output='$i'_'$run'.%J.out' >> $coevol_output_loc/"$i"_31species/"$i"_"$run".sh

        # COEVOL COMMAND
        echo "srun coevol -d $prank_alignment_loc \
        -t $tree_loc \
        -c $pheno_loc/"$i"_exon1_31species.lht -dsom exon1_"$i"_dsomggc"$run" " >> $coevol_output_loc/"$i"_31species/"$i"_"$run".sh 

        # SUBMIT THE TASK 
        cd $coevol_output_loc/"$i"_31species/
        sbatch $coevol_output_loc/"$i"_31species/"$i"_"$run".sh

done
done



