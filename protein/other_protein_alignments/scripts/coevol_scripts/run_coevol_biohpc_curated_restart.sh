#!/bin/bash

#define file locations
workdir=~/COEVOL/later_batch
tree_loc=~/COEVOL/trees/tree_30sp_for_coevol.txt
pheno_loc=~/COEVOL/pheno_coevol
coevol_output_loc=$workdir/results_biohpc_good_curated


#list of leftover proteins

declare -a arr=("GI" "body_mass_brain_mass_GI")
declare -a arr_run=(1 2 3)

cd $workdir/ccds_alignments/
ls good_curated | while read prot;

do

prot=$(echo $prot | sed -r 's/.best.nuc.phy//g')

for i in "${arr[@]}"
    do
    cd $coevol_output_loc/$prot/"$i"_30species
    
    for run in "${arr_run[@]}"
    do
   
      run_file=$coevol_output_loc/$prot/"$i"_30species/"$run"_"$prot"_"$i"_restart.sh 

        # HEADER
        echo '#!/bin/bash' > $run_file  
        echo '#SBATCH --error='$run'_'$prot'_'$i'.%J.err' >> $run_file
        echo '#SBATCH --output='$run'_'$prot'_'$i'.%J.out' >> $run_file
        echo '#SBATCH --partition=biohpc_gen_normal' >> $run_file
        echo '#SBATCH --clusters=biohpc_gen' >> $run_file
        echo '#SBATCH --time=2-00:00:00' >> $run_file
        echo '#SBATCH --get-user-env' >> $run_file
        echo '#SBATCH --job-name='$run'_'$prot'_'$i'' >> $run_file

        # COEVOL COMMAND
        echo '1' > $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc$run.run
        echo "srun coevol $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc$run" >> $run_file
        
        # SUBMIT THE TASK 
        sbatch $run_file
        done


done
done
