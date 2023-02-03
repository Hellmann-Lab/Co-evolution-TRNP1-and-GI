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


for i in "${arr[@]}"
    do
    tracecomp_summary=$coevol_output_loc/$prot/"$i"_30species/tracecomp_summary
    cd $coevol_output_loc/$prot/"$i"_30species
    rm -r tracecomp_summary
    mkdir -p $tracecomp_summary
     
    # BUILD IN A COMMAND THAT CHECKS CONVERGENCE
    tracecomp -x 1000 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc1 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc2 > $tracecomp_summary/tracecomp_"$i"_run1vs2.txt
    tracecomp -x 1000 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc1 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc3  > $tracecomp_summary/tracecomp_"$i"_run1vs3.txt
    tracecomp -x 1000 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc2 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc3  > $tracecomp_summary/tracecomp_"$i"_run2vs3.txt
         
    wait
     
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs2.txt | awk '$2<300 || $3>0.1' > $tracecomp_summary/tracecomp_"$i".txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run2vs3.txt | awk '$2<300 || $3>0.1' >> $tracecomp_summary/tracecomp_"$i".txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs3.txt | awk '$2<300 || $3>0.1' >> $tracecomp_summary/tracecomp_"$i".txt
    

        for run in "${arr_run[@]}"
        do
        echo '0' > $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc$run.run
        /dss/dsshome1/lxc03/di39qiz/readcoevol -x 1000 $coevol_output_loc/$prot/"$i"_30species/"$prot"_"$i"_dsomggc$run

        done

done
done
