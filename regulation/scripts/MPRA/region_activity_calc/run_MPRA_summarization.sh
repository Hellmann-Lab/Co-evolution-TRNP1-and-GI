#!/bin/bash

workdir=/data/share/htp/TRNP1/paper_data/regulation/scripts/MPRA/region_activity_calc
slurm=$workdir/slurms

mkdir -p $slurm
mkdir -p /data/share/htp/TRNP1/paper_data/regulation/data/MPRA/output/positional_activity
mkdir -p /data/share/htp/TRNP1/paper_data/regulation/data/MPRA/output/summarized_activity



declare -a arr=("intron" "upstream3" "upstream2" "upstream1" "exon2" "downstream" "exon1")


## now loop through the above array
for i in "${arr[@]}"
do

         # MAKING THE HEADER
        echo '#!/bin/bash' >$slurm/$i.sh  
        echo '#SBATCH -n 1' >>$slurm/$i.sh 
        echo '#SBATCH --error='$i'.%J.err' >>$slurm/$i.sh 
        echo '#SBATCH --output='$i'.%J.out' >>$slurm/$i.sh

        # RSCRIPT COMMAND
        # blast target (human protein) to query genomes
         
echo "Rscript $workdir/MPRA_summarization_function.R $i" >> $slurm/$i.sh


# SUBMIT THE TASK 

cd $slurm	
sbatch $i.sh
done