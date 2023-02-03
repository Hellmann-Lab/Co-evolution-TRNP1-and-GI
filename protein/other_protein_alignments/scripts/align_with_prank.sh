#!/bin/bash

base=protein/other_protein_alignments
prank_out=$base/prank_out
#prank_input=$base/prank_input

mkdir -p $prank_out/aln_forCoevol
mkdir -p $prank_out/aln_forChecking

cd $base
ls prank_input | while read i;

#for i in $fasta_list

do

slurmpfx=$base/slurms/slurm_$i.%J


sbatch --error=$slurmpfx.err --output=$slurmpfx.out --wrap="prank -d=$base/prank_input/$i \
-t=$base/tree_30sp_for_coevol.txt \
-o=$prank_out/aln_forChecking/$i -f='fasta' -translate"


sbatch --error=$slurmpfx.err --output=$slurmpfx.out --wrap="prank -d=$base/prank_input/$i \
-t=$base/tree_30sp_for_coevol.txt \
-o=$prank_out/aln_forCoevol/$i -f='phylips' -translate"

done

