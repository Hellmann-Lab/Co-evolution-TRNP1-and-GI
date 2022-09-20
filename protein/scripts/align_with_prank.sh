#!/bin/bash

#45 species for PAML
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_45species_forPAML_longer.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_45sp_forPAML.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_45species_forPAML_longer -f='phylips' -translate -once

prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_45species_forPAML_longer.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_45sp_forPAML.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_45species_forPAML_longer -f='fasta' -translate -once


# now align for Coevol!

#better ferret
prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_30species_forCoevol_longer.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_30sp_forCoevol.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol_longer -f='phylips' -translate -once


prank -d=/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_30species_forCoevol_longer.fa \
-t=/data/share/htp/TRNP1/paper_data/protein/trees/tree_TRNP1_coding_30sp_forCoevol.txt \
-o=/data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol_longer -f='fasta' -translate -once



