#!/bin/bash

#define file locations
workdir=~/COEVOL/later_batch
tree_loc=~/COEVOL/trees/tree_30sp_for_coevol.txt
pheno_loc=~/COEVOL/pheno_coevol
coevol_output_loc=$workdir/results_biohpc_good_noCuration


#list of 113 proteins


declare -a arr=("GI" "body_mass" "brain_mass" "body_mass_brain_mass_GI")
declare -a arr_run=(1 2 3)

cat $workdir/ccds_alignments/good_noCuration.tsv | while read prot;
do

prank_alignment_loc=$workdir/ccds_alignments/good_noCuration/"$prot".best.nuc.phy

if ! [ -f $prank_alignment_loc ]; then echo "alignment_missing in $prank_alignment_loc" > $coevol_output_loc/"$prot"_error_alignment.txt; else

mkdir -p $coevol_output_loc/"$prot"

for i in "${arr[@]}"
    do
    mkdir -p $coevol_output_loc/$prot/"$i"_30species
    tracecomp_summary=$coevol_output_loc/$prot/"$i"_30species/tracecomp_summary
    mkdir -p $tracecomp_summary
    cd $coevol_output_loc/$prot/"$i"_30species
    
    for run in "${arr_run[@]}"
    
    do
   
      run_file=$coevol_output_loc/$prot/"$i"_30species/"$run"_"$prot"_"$i".sh 

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
        echo "srun coevol -d $prank_alignment_loc \
        -t $tree_loc \
        -c $pheno_loc/"$i"_30sp.lht -dsom "$prot"_"$i"_dsomggc"$run" " >> $run_file 

        # SUBMIT THE TASK 
        sbatch $run_file
        
        done

done
fi
done


