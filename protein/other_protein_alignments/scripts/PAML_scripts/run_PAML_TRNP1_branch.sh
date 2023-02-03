#!/bin/bash

workdir=protein/other_protein_alignments/PAML/PAML_out_branch_model
alns=protein/fastas/prank_output


declare -a arr=("1")

prot=TRNP1

mkdir -p $workdir/$prot

for i in "${arr[@]}"
do

mkdir -p $workdir/$prot/$i


echo "seqfile = $alns/prank_"$prot"_coding_30species_forCoevol_longer.best.nuc.phy"  > $workdir/$prot/$i/codonml_branch_model_"$i".ctl
echo "treefile = protein/other_protein_alignments/tree_30sp_for_coevol.txt" >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl
echo "outfile = $workdir/$prot/$i/"$i"_model_output.txt" >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl
cat /data/share/htp/TRNP1/other_protein_alignments/blat_322cand/PAML/constant_codonml_branch_model_"$i".txt >> $workdir/$prot/$i/codonml_branch_model_"$i".ctl


cd $workdir/$prot/$i
sbatch --wrap="codeml codonml_branch_model_"$i".ctl"

done
