#!/bin/bash
# 45 species with full names
prank -d=protein/fastas/TRNP1_coding_seqs_45sp.fa \
-t=protein/trees/tree_TRNP1_coding_45sp.txt \
-o=protein/fastas/prank_output/prank_TRNP1_coding_45species -f='fasta' -translate


#the 45 species for PAML (adjusted species names)
prank -d=protein/fastas/TRNP1_coding_seqs_45species_forPAML.fa \
-t=protein/trees/tree_TRNP1_coding_45sp_forPAML.txt \
-o=protein/fastas/prank_output/prank_TRNP1_coding_45species_forPAML -f='phylips' -translate -once


#the 31 species with phenotype data for Coevol
prank -d=protein/fastas/TRNP1_coding_seqs_31species_forCoevol.fa \
-t=protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt \
-o=protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol -f='phylips' -translate 

#fasta for inspecting the alignment in was.bi
prank -d=protein/fastas/TRNP1_coding_seqs_31species_forCoevol.fa \
-t=protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt \
-o=protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol -f='fasta' -translate 
