#!/bin/bash
# actually do not need the -showall flag that returns reconstructed ancestral sequences and the tree
# 45 species with full names
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_45sp.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_45sp.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_45species -f='fasta' -translate


#the 45 species for PAML (adjusted species names)
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_45species_forPAML.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_45sp_forPAML.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_45species_forPAML -f='phylips' -translate -once


#the 31 species with phenotype data for Coevol
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_31species_forCoevol.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol -f='phylips' -translate 

#fasta for inspecting the alignment in was.bi
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_31species_forCoevol.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_31species_forCoevol -f='fasta' -translate 
