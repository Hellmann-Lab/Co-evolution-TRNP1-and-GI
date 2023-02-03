#!/bin/bash

#where my PAML outputs are located with 2 folders called NS7 and NS8
workdir=/data/share/htp/TRNP1/paper_data/protein/PAML

declare -a arr=("NS7" "NS8")
prot=TRNP1

for i in "${arr[@]}"
do

lnL=$(cat $workdir/$i/"$i"_model_output.txt | grep '^lnL')

if [ $i == "NS8" ]; then

  p0=$(cat $workdir/$i/"$i"_model_output.txt | grep 'p0 =') #| grep -o -E '[0-9]+' ) #| sed -e 's/\*//')

else

  p0="  p0 =   1  p =   1 q =   1"

fi

echo -e "$prot\t$i\t$lnL\t$p0" >> $workdir/M8_vs_M7.txt

done 
