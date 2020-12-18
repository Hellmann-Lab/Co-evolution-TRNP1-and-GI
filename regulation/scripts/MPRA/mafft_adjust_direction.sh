#!/bin/bash
#SBATCH -n 1
#SBATCH --error=mafft_adjust_direction.%J.err
#SBATCH --output=mafft_adjust_direction.%J.out
#SBATCH --workdir=regulation/scripts/MPRA/slurms


declare -a arr=("intron" "upstream3" "upstream2" "upstream1" "exon2" "downstream" "exon1")


## now loop through the above array
for i in "${arr[@]}"
do

/home/zane/TRNP1/mafft7/mafft-linux64/mafft.bat --adjustdirection --reorder --genafpair regulation/fastas/separated_regions/TRNP1_"$i"_seqs.fa \
> regulation/fastas/mafft_regions/mafft_"$i"_adj.fa 

done


