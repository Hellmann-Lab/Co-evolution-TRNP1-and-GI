library(Biostrings)
library(stringr)
library(tidyverse)
library(plyranges)
library(foreach)
library(doParallel)


prot_seqs_sp <- list.files("protein/other_protein_alignments/self_matches_fa", recursive = T, full.names = T)



#for the cases where multiple sequences were matching equally good, count different positions
diff_per_sp<-list()
length_diff_per_sp<-list()


for (sp in 1:length(prot_seqs_sp)) {

  seq_list<-readRDS(prot_seqs_sp[sp])
  species<-word(prot_seqs_sp[sp],9,9,"/") %>% gsub(".rds","",.)
  seq_differences_list<-list()
  
  for (Q in 1:length(seq_list)){
    
  #in case the lengths are different, keep the shorter one; write down the longer one
   if (length(unique(width(seq_list[[Q]])))!=1){
     
     seq_list[[Q]]<-seq_list[[Q]][width(seq_list[[Q]])==min(width(seq_list[[Q]]))]
     length_diff_per_sp[[paste0(species,names(seq_list[Q]))]]<-max(width(seq_list[[Q]]))
     
   } else {
     
     pos_split_diff<-as.data.frame(strsplit(x=as.character(seq_list[[Q]]), split=character(0))) %>%
       rownames_to_column("pos_seq") %>%
        pivot_longer(2:ncol(.)) %>% 
        group_by(pos_seq) %>%
        mutate(n=length(unique(value))) %>%
        filter(n>1) 
   }
  
  if (nrow(pos_split_diff)!=0){
    
    seq_differences_list[[names(seq_list[Q])]]<-pos_split_diff
  }
    }

  if (length(seq_differences_list)==0){
    diff_per_sp[[species]]<-NULL
  }
  else {
    seq_differences_df<-bind_rows(seq_differences_list, .id="protein") %>%
      group_by(protein) %>%
      summarise(n_positions_nucDiff=length(unique(pos_seq)))
  saveRDS(seq_differences_df, paste0("protein/other_protein_alignments/self_matches_diff/", species,".rds"))
  
  diff_per_sp[[species]]<-seq_differences_df
}
}

saveRDS(diff_per_sp,"protein/other_protein_alignments/seq_diffs_matchingSeqs.rds")

saveRDS(length_diff_per_sp,"protein/other_protein_alignments/seq_length_diffs_matchingSeqs.rds")
