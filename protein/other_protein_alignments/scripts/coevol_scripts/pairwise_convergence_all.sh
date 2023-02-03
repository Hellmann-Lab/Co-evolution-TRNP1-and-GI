#!/bin/bash

#define file locations
workdir=~/COEVOL/later_batch
tree_loc=~/COEVOL/trees/tree_30sp_for_coevol.txt
pheno_loc=~/COEVOL/pheno_coevol


#list of 113 proteins

declare -a arr=("GI" "body_mass_brain_mass_GI" "brain_mass" "body_mass")
declare -a arr_out=("good_curated" "good_noCuration")

for type in "${arr_out[@]}"
do

coevol_output_loc=$workdir/results_biohpc_"$type"

cd $workdir/ccds_alignments/
ls $type | while read prot;

do

prot=$(echo $prot | sed -r 's/.best.nuc.phy//g')

for i in "${arr[@]}"
    do
    tracecomp_summary=$coevol_output_loc/$prot/"$i"_30species/tracecomp_summary
    cd $coevol_output_loc/$prot/"$i"_30species
     
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs2.txt | awk '$2<300 || $3>0.1' > $tracecomp_summary/tracecomp_"$i"_run1vs2_unconv.txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run2vs3.txt | awk '$2<300 || $3>0.1' > $tracecomp_summary/tracecomp_"$i"_run2vs3_unconv.txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs3.txt | awk '$2<300 || $3>0.1' >> $tracecomp_summary/tracecomp_"$i"_run1vs3_unconv.txt


done
done
done
