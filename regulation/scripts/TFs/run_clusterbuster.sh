#!/bin/bash
#SBATCH --error=regulation/data/TFs/ClusterBuster/results/slurms/slurm-%J.err
#SBATCH --output=regulation/data/TFs/ClusterBuster/results/slurms/slurm-%J.out

workingdir=regulation/data/TFs

JASPAR_folder_name=$1
fasta=$2  
weighted=$3
cluster_cutoff=$4
motif_cutoff=$5
gap=$6


cd $workingdir/JASPAR_2020/JASPAR_collection_transposed/$JASPAR_folder_name
for i in * 
do 
echo ">$i";cat "$i"; done >  $workingdir/JASPAR_2020/ready_for_ClusterBuster/TF_psfm_jaspar2020_"$JASPAR_folder_name".txt


/opt/bin/ctrain -x $workingdir/JASPAR_2020/ready_for_ClusterBuster/TF_psfm_jaspar2020_"$JASPAR_folder_name".txt \
$workingdir/ClusterBuster/input_fastas/$fasta \
-c $workingdir/ClusterBuster/results/$JASPAR_folder_name/ctrain_matrix.txt  > \
$workingdir/ClusterBuster/results/$JASPAR_folder_name/ctrain_output.txt


#gap parameter and weights inferred from ctrain

if [[ "$weighted" == "TRUE" ]]
then
/opt/bin/cbust -c$cluster_cutoff -m$motif_cutoff $workingdir/ClusterBuster/results/$JASPAR_folder_name/ctrain_matrix.txt \
$workingdir/ClusterBuster/input_fastas/$fasta > \
$workingdir/ClusterBuster/results/$JASPAR_folder_name/cbust_weighted_c"$cluster_cutoff"_m"$motif_cutoff".txt


/opt/bin/cbust -c$cluster_cutoff -m$motif_cutoff -f1 $workingdir/ClusterBuster/results/$JASPAR_folder_name/ctrain_matrix.txt \
$workingdir/ClusterBuster/input_fastas/$fasta > \
$workingdir/ClusterBuster/results/$JASPAR_folder_name/cbust_weighted_c"$cluster_cutoff"_m"$motif_cutoff"_f1.txt

fi


#gap parameter set by the user, no weights

if [[ "$weighted" == "FALSE" ]]
then
/opt/bin/cbust -g$gap -c$cluster_cutoff -m$motif_cutoff $workingdir/JASPAR_2020/ready_for_ClusterBuster/TF_psfm_jaspar2020_"$JASPAR_folder_name".txt \
$workingdir/ClusterBuster/input_fastas/$fasta > \
$workingdir/ClusterBuster/results/$JASPAR_folder_name/cbust_g"$gap"_c"$cluster_cutoff"_m"$motif_cutoff".txt


/opt/bin/cbust -g$gap -c$cluster_cutoff -m$motif_cutoff -f1 $workingdir/JASPAR_2020/ready_for_ClusterBuster/TF_psfm_jaspar2020_"$JASPAR_folder_name".txt \
$workingdir/ClusterBuster/input_fastas/$fasta > \
$workingdir/ClusterBuster/results/$JASPAR_folder_name/cbust_g"$gap"_c"$cluster_cutoff"_m"$motif_cutoff"_f1.txt

fi



#subset cbust output
bash regulation/scripts/TFs/cbust_subset_per_species.sh $JASPAR_folder_name $weighted $cluster_cutoff $motif_cutoff $gap 


