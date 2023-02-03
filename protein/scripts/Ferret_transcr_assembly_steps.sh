#!/bin/bash

#1
######### map with bwa2 ##########
r1=/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trimmed/R1.fastq
r2=/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trimmed/R2.fastq
bam=/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/ferret_trnp1_bwa.bam
genomeFa=/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/genome/musFur1.fa

sbatch --mem=100G --wrap="bwa-mem2 index -p musFur1 $genomeFa"
wait 
sbatch -c 20  --mem=100G --wrap="bwa-mem2 mem -t 20 -M musFur1 $r1 $r2 | \
samtools fixmate -m - - | samtools sort -o $bam"

## consider cleaning bam only mate pairs ...

#2
######### Trinity ################
#Trinity Reference guided non referencguided problematic, maybe try oases
trinityDir=/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trinity_kmer_bwa/

sbatch -c 5  --mem=15G --wrap="Trinity --seqType fq --left $r1 --right $r2 \
--output $trinityDir \
--min_contig_length 300 \
--normalize_reads \
--trimmomatic \
--CPU 5 \
--max_memory 10G \
--min_kmer_cov 2 \
--genome_guided_bam $bam \
--genome_guided_max_intron 10000"

#basically include all ends of the pcr sequences, and a reverseComplement version of the starts. then cut AFTER it's occurence 
sed -E 's/(CACGGCTCTCCAGCTCCGAC|AAATGGGGTGGGGTG|CGGGGTCGGAGTCAAGGTCG|AGGCTGAGGTGCGGTCTG|GCCACAAGGGACCCAGGA|TCTCCTTCCTGCCCAGG|AGCCAATCAGAGACTGGC|CCGACCTTGACTCCGA)([^.])/\1\n>new\n\2/g' $trinityDir/Trinity-GG.fasta > $trinityDir/Trinity-GG-split.fasta

#cut BEFORE starts and reverseComplement ends
sed -E 's/(GTCGGAGCTGGAGAGCCGTG|CAGCACCCCACCCCATT|CGACCTTGACTCCGACCCC|CAGACCGCACCTCAGC|TCCTGGGTCCCTTGTG|CCTGGGCAGGAAGGAGA|GCCAGTCTCTGATTGG|TCGGAGTCAAGGTCGG)/\n>new\n\1/g' $trinityDir/Trinity-GG-split.fasta > $trinityDir/Trinity-GG-split2.fasta

sbatch -c 5  --mem=10G --wrap="Trinity --seqType fa --single $trinityDir/Trinity-GG-split2.fasta \
--output $trinityDir/trinity-final \
--min_contig_length 300 \
--CPU 5 \
--max_memory 10G"

#make a blast library
makeblastdb  -input_type fasta -dbtype nucl -in $trinityDir/trinity-final.Trinity.fasta -out $trinityDir/finalBlastDB

#make PCR library blast db
makeblastdb  -input_type fasta -dbtype nucl -in /data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/MusPut_Trnp1_UCSC_in_silico_PCR.fa -out /data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/PCRblastDB


#switch to Ferret_cutContig_makeFa.R
