#!/bin/bash

#define file locations
coevol_input_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/scripts
prank_alignment_loc=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol.best.nuc.phy
tree_loc=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_30sp_forCoevol.txt
pheno_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/pheno_data/
coevol_output_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/results/no_dolphin


declare -a arr=("GI" "GI_EQ_body_mass" "GI_brain_mass" "brain_mass")
declare -a arr_run=(1 2 3)



for i in "${arr[@]}"
do
mkdir $coevol_output_loc/"$i"_30species
done


for i in "${arr[@]}"
do

for run in "${arr_run[@]}"

   do
        # HEADER
        echo '#!/bin/bash' > $coevol_output_loc/"$i"_30species/"$i"_"$run".sh  
        echo '#SBATCH -n 1' >> $coevol_output_loc/"$i"_30species/"$i"_"$run".sh
        echo '#SBATCH --error='$i'_'$run'.%J.err' >> $coevol_output_loc/"$i"_30species/"$i"_"$run".sh 
        echo '#SBATCH --output='$i'_'$run'.%J.out' >> $coevol_output_loc/"$i"_30species/"$i"_"$run".sh

        # COEVOL COMMAND
        echo "srun coevol -d $prank_alignment_loc \
        -t $tree_loc \
        -c $pheno_loc/"$i"_exon1_31species.lht -dsom exon1_"$i"_dsomggc"$run" " >> $coevol_output_loc/"$i"_30species/"$i"_"$run".sh 

        # SUBMIT THE TASK 
        cd $coevol_output_loc/"$i"_30species/
        sbatch $coevol_output_loc/"$i"_30species/"$i"_"$run".sh

done
done



