#!/bin/bash

workdir=/data/share/htp/TRNP1/paper_data/protein
output=$workdir/fastas/TRNP1_coding/coding_seqs_sep2_full
slurm=$output/slurms



cat /data/share/htp/TRNP1/paper_data/protein/data/additional_genomes_full.txt | while read i j k;

    do
         # MAKING THE HEADER
        echo '#!/bin/bash' >$slurm/$i.sh  
        echo '#SBATCH -n 1' >>$slurm/$i.sh 
        echo '#SBATCH --error='$i'.%J.err' >>$slurm/$i.sh 
        echo '#SBATCH --output='$i'.%J.out' >>$slurm/$i.sh

        # THE ACTUAL COMMANDS
        # blast target (human protein) to query genomes
         
echo "Rscript $workdir/scripts/collect_coding_seqs2/collect_coding_seqs2_function.R $i $j $k $output" >> $slurm/$i.sh
#cleanup
#echo "rm $j/$i.seqs.fasta*" >> $slurm/$i.sh

# SUBMIT THE TASK 
cd $slurm	
sbatch $slurm/$i.sh
done