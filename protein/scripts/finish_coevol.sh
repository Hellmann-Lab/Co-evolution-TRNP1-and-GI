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
prot="exon1"


for i in "${arr[@]}"
    do
    tracecomp_summary=$coevol_output_loc/"$i"_"$sp"species/tracecomp_summary
    cd $coevol_output_loc/"$i"_"$sp"species
    rm -r tracecomp_summary
    mkdir -p tracecomp_summary
     
    # BUILD IN A COMMAND THAT CHECKS CONVERGENCE
    tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc1 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc2 > $tracecomp_summary/tracecomp_"$i"_run1vs2.txt
    tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc1 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc3  > $tracecomp_summary/tracecomp_"$i"_run1vs3.txt
    tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc2 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc3  > $tracecomp_summary/tracecomp_"$i"_run2vs3.txt
         
    wait
     
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs2.txt | awk '$2<300 || $3>0.1' > $tracecomp_summary/tracecomp_"$i".txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run2vs3.txt | awk '$2<300 || $3>0.1' >> $tracecomp_summary/tracecomp_"$i".txt
    sed '1,2d' $tracecomp_summary/tracecomp_"$i"_run1vs3.txt | awk '$2<300 || $3>0.1' >> $tracecomp_summary/tracecomp_"$i".txt
    

      for run in "${arr_run[@]}"
      do
      echo '0' > $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc$run.run
      readcoevol -x 1000 $coevol_output_loc/"$i"_"$sp"species/"$prot"_"$i"_dsomggc$run

      done

done
done
