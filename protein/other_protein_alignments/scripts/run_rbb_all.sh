#!/bin/bash

#make a bash script that for each entrance generates a 2bit file and does reciprocal best blat of human protein sequences and back
TRNP1=protein

mkdir -p $TRNP1/other_protein_alignments/blat_output_full

cat data/coevol_species_genomes.txt | while read i j k;
do

cd $j

if [ ! -f $k.2bit ]; then
    echo "Generating 2bit file"
    /opt/bin/faToTwoBit $k $k.2bit
fi

mkdir -p $TRNP1/other_protein_alignments/blat_output_full/$i
slurmpfx=$TRNP1/other_protein_alignments/blat_output_full/$i/slurm_$i.%J


sbatch -J $i --mem=15G --error=$slurmpfx.err --output=$slurmpfx.out --wrap="$TRNP1/other_protein_alignments/reciprocal_best_blat.R --ref2bit=$TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/Homo_sapiens/hg38.p12.2bit --other2bit=$k.2bit --refprotein=$TRNP1/other_protein_alignments/CCDS/CCDS_protein.current.faa.gz --outdir=$TRNP1/other_protein_alignments/blat_output_full/$i --outbase=$i --prop_prot_length=0.3"

done
