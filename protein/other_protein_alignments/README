README on evolutionary analysis of control protein sequences:
From reciprocal best blat to multiple sequence alignment and running Coevol and PAML on all good control proteins.

## IMPORTANT OUTPUT FILES

The good alignments used for evolutionary analyses are under `prank_out/good_noCuration` and `prank_out/good_curated` as .phy files (necessary for PAML and Coevol; nucleotide level), and all together under `prank_out/aln_forChecking/` in fasta format for checking (nuc and prot alignments.)

All the important output files can be found in `summaries` folder; they are all loaded and summarized as in `scripts/put_all_results_together.R`.


## WORKFLOW

1) Download the necessary data 

Download human CCDS using the following command:

`for i in BuildInfo.current.txt CCDS.current.txt CCDS2Sequence.current.txt CCDS2UniProtKB.current.txt CCDS_attributes.current.txt CCDS_exons.current.txt CCDS_nucleotide.current.fna.gz CCDS_protein.current.faa.gz CCDS_protein_exons.current.faa.gz; 
do sbatch -J $i --cpus-per-task=10 --wrap="wget https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/$i"; done`

And download the genomes listed under `../data/coevol_species_genomes.txt`


2) Run RBB using bash `scripts/run_rbb_all.sh` which calls `scripts/reciprocal_best_blat.R` script and generates an output directory `blat_output_full` for all human CCDS RBB matches per target genome.


3) Select similar genes to TRNP1 in `scripts/select_candidates.R`. 


4) I rerun RBB on these promising 322 proteins to be able to upload as much as possible of this smaller data on github, script under `scripts/run_rbb_all_322.sh`

5) Sequence extraction `scripts/extract_seqs.R`, further filtering `scripts/compare_matching_seqs.R` and generation of the final set for MSAs `scripts/final_controls.R`

6) Multiple sequence alignment using PRANK `scripts/align_with_prank.sh`
After inspecting the alignments, i got 112 good ones. Their subsetting is in the next script: `scripts/subset_alignments_forCoevol.R`

The ones that I decided are worth redoing with other genomes / curating, were redone when using alternative genome versions for 5 species gorilla gorGor5.fa, dolphin GCF_011762595.1_mTurTru1, wild boar GCF_000003025.6_Sscrofa11.1, rhesus macaque GCF_003339765.1_Mmul_10, olive baboon GCA_000264685.2_Panu_3.0) and redoing the alignment. I got out 22 additional good quality protein sequence alignments.


7) Running Coevol
See the README under `scripts/coevol_scripts`

8) Checking out additional properties:
alignment quality, total tree length and variation: `scripts/summarise_tree_lengths.R`

9) PAML site model: folder PAML. The order goes as follows (all within `scripts/PAML_scripts` folder):

`run_PAML_final_proteins_branch.sh` 
`run_PAML_TRNP1_branch.sh`
`summarise_branch_out.sh`


10) Finally, read and summarize all results under `scripts/put_all_results_together.R`

