#!/bin/bash

#define file locations
coevol_output_loc=/data/share/htp/TRNP1/paper_data/protein/coevol/results/all
mkdir $coevol_output_loc/tracecomp_summary
mkdir $coevol_output_loc/correlation_output

tracecomp_summary=$coevol_output_loc/tracecomp_summary
cor_summary=$coevol_output_loc/correlation_output

sp=31


#declare the test groups
declare -a arr=("GI" "GI_EQ_body_mass" "brain_mass" "GI_body_mass" "GI_brain_mass" "GI_EQ" "body_mass" "EQ" "EQ_body_mass")
#"GI_brain_mass"
declare -a arr_run=(1 2 3)


for i in "${arr[@]}"
do

#first, calculate convergence of the runs

mkdir $tracecomp_summary/"$i"_"$sp"sp

tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc1 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc2 > $tracecomp_summary/"$i"_"$sp"sp/tracecomp_"$i"_run1vsrun2.txt
tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc1 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc3 > $tracecomp_summary/"$i"_"$sp"sp/tracecomp_"$i"_run1vsrun3.txt
tracecomp -x 1000 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc2 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc3 > $tracecomp_summary/"$i"_"$sp"sp/tracecomp_"$i"_run2vsrun3.txt



#second, calculate and export the correlation coefficients, posterior probabilities and the variables included
for run in "${arr_run[@]}"
do

#readcoevol -x 1000 $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run"

#export parts of the .cov table 
mkdir $cor_summary/"$i"_"$sp"sp
# #from https://stackoverflow.com/questions/45314145/print-lines-after-a-pattern-until-second-occurrence-of-a-different-pattern
# sed -n '/^correlation coefficients/,/^ *$/{/^correlation coefficients/N;p}' $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run".cov > $cor_summary/"$i"_"$sp"sp/cors_"$i"_run"$run".txt
# sed -n '/^posterior probs/,/^ *$/{/^posterior probs/N;p}' $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run".cov | awk -v N=2 '{print}/^ *$/&&--N<=0{exit}' > $cor_summary/"$i"_"$sp"sp/pp_"$i"_run"$run".txt
# awk -v N=1 '{print}/^ *$/&&--N<=0{exit}' $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run".cov | sed '1d' > $cor_summary/"$i"_"$sp"sp/names_"$i"_run"$run".txt
# 
# sed -n '/^partial correlation coefficients/,/^ *$/{/^ *$/N;p}' $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run".cov > $cor_summary/"$i"_"$sp"sp/partial_cors_"$i"_run"$run".txt
#sed -n '/^posterior probs/,/^ *$/{/^posterior probs/N;p}' $coevol_output_loc/"$i"_"$sp"species/exon1_"$i"_dsomggc"$run".cov | awk -v N=2 '{print}/^ *$/&&--N<=0{exit}' > $cor_summary/"$i"_"$sp"sp/pp_"$i"_run"$run".txt



#extract the ln probabilities of the model 
files=($coevol_output_loc/"$i"_"$sp"species/*.err)
awk '/exit/{getline; print}' "${files[$run-1]}" > $cor_summary/"$i"_"$sp"sp/ln_prob_"$i"_run"$run".txt


done
done