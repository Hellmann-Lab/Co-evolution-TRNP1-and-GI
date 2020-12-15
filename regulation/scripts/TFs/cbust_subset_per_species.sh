#!/bin/bash

JASPAR_folder_name=$1
weighted=$2
cluster_cutoff=$3
motif_cutoff=$4
gap=$5

#specify file locations
workdir=/data/share/htp/TRNP1/paper_data/regulation/data/TFs/ClusterBuster
res=$workdir/results/$JASPAR_folder_name #where jaspar_output files are to be found
subs=$workdir/results/$JASPAR_folder_name/subset_cbust #where to put all of this (keep this script in a subfolder "scripts" of this folder)


if [[ "$weighted" == "TRUE" ]]
then
jaspar_output=cbust_weighted_c"$cluster_cutoff"_m"$motif_cutoff".txt
jaspar_output_f1=cbust_weighted_c"$cluster_cutoff"_m"$motif_cutoff"_f1.txt #the same as jaspar_output, just in a different format
fi

if [[ "$weighted" == "FALSE" ]]
then
jaspar_output=cbust_g"$gap"_c"$cluster_cutoff"_m"$motif_cutoff".txt
jaspar_output_f1=cbust_g"$gap"_c"$cluster_cutoff"_m"$motif_cutoff"_f1.txt #the same as jaspar_output, just in a different format
fi


#remove old scripts in case I haven't deleted them
mkdir $subs
mkdir $subs/scripts

rm $subs/scripts/subsetting.sh
rm $subs/scripts/subsetting2.sh
rm -r $subs/motifs
rm -r $subs/scores
rm -r $subs/sequences
rm -r $subs/species

#selecting species names and dropping everything but the name
mkdir $subs/species
sed -n '/^>/p' $res/$jaspar_output_f1 | sed 's/>//g;s/[[:blank:]]//g;s/[[:space:]]//g;s/[0-9]*//g;s/bp//g;s/)//g;s/(//g'  > $subs/species/species_output.txt

#generate bash commands to subset 
mkdir $subs/scores
while IFS= read -r line; do
   echo "sed -n '/"$line"/,/^$/p' $res/$jaspar_output_f1 > $subs/scores/"$line".txt" >> $subs/scripts/subsetting.sh
done < $subs/species/species_output.txt

#run subsetting
bash $subs/scripts/subsetting.sh

#generate a list of TF cluster sequences
mkdir $subs/sequences
awk '/^>|^TT|^TA|^TC|^TG|^AA|^AT|^AC|^AG|^CC|^CA|^CT|^CG|^GG|^GC|^GT|^GA/' $res/$jaspar_output > $subs/sequences/cluster_sequences.txt

#select the specific motifs per species
mkdir $subs/motifs
mkdir $subs/motifs/tmp
while IFS= read -r line; do
   echo "sed -n '/"$line"/,/^>/p' $res/$jaspar_output > $subs/motifs/tmp/"$line".txt" >> $subs/scripts/subsetting2.sh
done < $subs/species/species_output.txt

bash $subs/scripts/subsetting2.sh

mkdir $subs/motifs/motif_tables
while IFS= read -r line; do
   awk '/[[:blank:]]-[[:blank:]]|+/' $subs/motifs/tmp/"$line".txt > $subs/motifs/motif_tables/"$line".txt
done < $subs/species/species_output.txt



