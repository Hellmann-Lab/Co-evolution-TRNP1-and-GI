#!/bin/bash

#where my PAML outputs are located with 2 folders called NS7 and NS8
workdir=protein/other_protein_alignments/PAML/PAML_out_branch_model

rm $workdir/../branch_model_out.txt

declare -a arr=("1")
declare -a farr=("good_noCuration" "good_curated")


for f in "${farr[@]}"
do

cp $workdir/protein_list_"$f".txt $workdir/protein_list_"$f"_withTRNP1.txt
echo "TRNP1" >> $workdir/protein_list_"$f"_withTRNP1.txt
cat $workdir/protein_list_"$f"_withTRNP1.txt | while read prot;

do

for i in "${arr[@]}"
do

dN=$(cat $workdir/$prot/$i/"$i"_model_output.txt | grep '^tree length for dN:' | sed 's/tree length for dN://' | sed 's/ //g')
dS=$(cat $workdir/$prot/$i/"$i"_model_output.txt | grep '^tree length for dS:' | sed 's/tree length for dS://' | sed 's/ //g')


echo -e "$prot\t$dN\t$dS" >> $workdir/../branch_model_out.txt

done 
done
done

sed -i '1 i\Protein\tdN\tdS' $workdir/../branch_model_out.txt 

