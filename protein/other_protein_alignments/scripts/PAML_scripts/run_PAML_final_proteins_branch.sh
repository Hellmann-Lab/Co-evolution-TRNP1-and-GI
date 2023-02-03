#!/bin/bash

workdir=protein/other_protein_alignments/PAML/PAML_out_branch_model
alns=protein/other_protein_alignments/prank_out

declare -a farr=("good_noCuration" "good_curated")

declare -a arr=("1")


for f in "${farr[@]}"
do

cd $alns/$f
for aln in $(ls *best.nuc.phy ); do echo $aln | cut -d'_' -f 2; done > $workdir/protein_list_"$f".txt
#sort $workdir/protein_list_"$f".txt | uniq > $workdir/protein_list_"$f".txt
cat $workdir/protein_list_"$f".txt | while read prot; #can use the following part before while arg: | tail -n 75

do

mkdir -p $workdir/$prot

for i in "${arr[@]}"
do

mkdir -p $workdir/$prot/$i


echo "seqfile = $alns/$f/$prot" > $workdir/$prot/$i/codonml_branch_model_"$i".ctl
echo "treefile = protein/other_protein_alignments/tree_30sp_for_coevol.txt" >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl
echo "outfile = $workdir/$prot/$i/"$i"_model_output.txt" >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl
#echo \n >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl
cat protein/other_protein_alignments/PAML/constant_codonml_branch_model_"$i".txt >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl


cd $workdir/$prot/$i
sbatch --wrap="codeml codonml_branch_model_"$i".ctl"

done
done
done