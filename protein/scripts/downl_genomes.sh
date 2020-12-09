#!/bin/bash

workdir=/data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/

declare -a arr=("NCBI" "ENSEMBL")


for i in "${arr[@]}"
do

  cat /data/share/htp/TRNP1/paper_data/protein/data/"$i"_genomes.txt | while read i j;

    do
        mkdir -p $workdir/$i
        cd $workdir/$i
        sbatch -J $i --cpus-per-task=10 --wrap="wget $j"
        
  done
done