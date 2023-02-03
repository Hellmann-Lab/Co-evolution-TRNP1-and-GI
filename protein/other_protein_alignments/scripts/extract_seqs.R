
library(Biostrings)
library(stringr)
library(tidyverse)
library(plyranges)
library(foreach)
library(doParallel)

#what I will need to do:
# 1) cut out blat matches with the same overlap to the human sequence and same blosum scores
#quantify their differences - my guess is that they might be nearly identical and not well assambled due to rather bad and short scaffolds where the same seq cames up multiple times 


#species
species_df<-read_delim("protein/data/coevol_species_genomes.txt", delim =" ", col_names = F) %>%
  filter(X1!="Homo_sapiens")

#proteins
CCDS_322<-readAAStringSet("protein/other_protein_alignments/CCDS_protein.current.322.fa")


blat_overlap_list <- list.files("protein/other_protein_alignments/results", pattern = "_blat_overlap", recursive = T, full.names = T) 

blat_toOther_list <- list.files("protein/other_protein_alignments/results", pattern = "_blat_to_other_filt", recursive = T, full.names = T)



# 1) compare sequences of blat matches ####

registerDoParallel(30)  # use multicore, set to the number of cores

foreach (sp=1:length(blat_overlap_list)) %dopar% {

  #preparation
  blat_sp_overlap<-readRDS(blat_overlap_list[sp])
  blat_sp_toOther<-data.table::fread(blat_toOther_list[sp]) %>%
    filter(index_other %in% blat_sp_overlap$index_other) %>%
    mutate(T_strand=substr(strand,2,2),
           new_index=1:nrow(.)) 
  
  species<-word(blat_overlap_list[sp],9,9,"/")
  species_base<-species_df$X2[species_df$X1==species]
  species_gen<-species_df$X3[species_df$X1==species]
  
  input_seq <-Biostrings::readDNAStringSet(paste0(species_base,"/",species_gen))
  
  #cut out the matching positions
  seq_list<-list()
  for (Q in unique(blat_sp_toOther$Q_name)){

    seq_out<-DNAStringSet()
    
    for (i in blat_sp_toOther$new_index[blat_sp_toOther$Q_name==Q]){
      
      subject_seq<-input_seq[grepl(paste0("^",blat_sp_toOther$T_name[i]," |^",blat_sp_toOther$T_name[i],"$"),names(input_seq))]
      
      subject_cut<-subseq(subject_seq, start=blat_sp_toOther$T_start[i]+1, end=blat_sp_toOther$T_end[i])
      
      if (blat_sp_toOther$T_strand[i]=="-"){
        subject_cut<-Biostrings::reverseComplement(subject_cut)
      }
      
      names(subject_cut)<-paste0(blat_sp_toOther$T_name[i])
      seq_out<-c(seq_out,subject_cut)
    }
  seq_list[[Q]]<-seq_out
  }
  saveRDS(seq_list, paste0("protein/other_protein_alignments/self_matches_fa/", species,".rds"))
}