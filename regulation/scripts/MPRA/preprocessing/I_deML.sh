#!/bin/bash
#SBATCH --workdir=/data/share/htp/TRNP1/MPRA_full_300919/MPRA_data/slurms
#SBATCH --error=MPRA_I_deML_%J.err
#SBATCH --output=MPRA_I_deML_%J.out
#SBATCH --cpus-per-task=5
#SBATCH --mem=32G

# Variables #
maindir="/data/share/htp/TRNP1/MPRA_full_300919/MPRA_data/"
fqdir="fastq"
demldir="deML"
declare -a lane=("lane6" "lane7" "lane8")


mkdir -p ${maindir}/fastqc
mkdir -p ${maindir}/${demldir}/demult
mkdir -p ${maindir}/${demldir}/demult/files
mkdir -p ${maindir}/${demldir}/demult/deML_summaries


# deML # 
for i in "${lane[@]}"
do
deML -i ${maindir}/${demldir}/MPRA_indices.txt \
-if1 ${maindir}/${fqdir}/"$i"_R2.fastq.gz \
-if2 ${maindir}/${fqdir}/"$i"_R3.fastq.gz \
-f ${maindir}/${fqdir}/"$i"_R1.fastq.gz \
-o ${maindir}/${demldir}/demult/files/demult_"$i" \
-s ${maindir}/${demldir}/demult/deML_summaries/summary_"$i".txt \
-e ${maindir}/${demldir}/demult/deML_summaries/error_"$i".txt


# fastqc #
fastqc ${maindir}/${fqdir}/"$i"* \
-o ${maindir}/fastqc/ \
-t ${SLURM_CPUS_PER_TASK}

done
