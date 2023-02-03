library(tidyverse)
library(ggtree)
library(Biostrings)

# try to calculate the average branch length across the sites (per protein)
#functions

calculate_tree_length<-function(df, pos_aln, tree){
  vec_sp<-df %>% filter(pos_alignment==pos_aln) %>% pull(., species)
  reduced_tree<-ape::drop.tip(tree, tree$tip.label[tree$tip.label %in% vec_sp])
  return(sum(reduced_tree$edge.length))
}


summarize_inform_sites<-function(tree_loc, aln_loc){
  
  #branch length sum
  tr<-read.tree(tree_loc)
  tr$tip.label<-word(tr$tip.label,1,1,"_")
  
  aln<-readAAMultipleAlignment(aln_loc)
  
  alignment_df<-strsplit(x=as.character(aln), split=character(0)) %>%
    as.data.frame() %>%
    rownames_to_column("pos_alignment")
  
  alignment_df_summarized<-alignment_df %>%
    tidyr::pivot_longer(cols = 2:ncol(.), names_to="species", values_to="AA") %>%
    dplyr::filter(AA=="-") 
  
  alignment_df_more<-alignment_df %>%
    rowwise() %>%
    mutate(tree_length=calculate_tree_length(df=alignment_df_summarized, pos_aln=pos_alignment, tree=tr)) %>%
    ungroup() %>%
    mutate(max_tree_length=max(tree_length)) %>%
    rowwise() %>%
    mutate(rel_tree_length=tree_length/max_tree_length) %>%
    ungroup()
  
  #variable sites ####
  pos_split_diff<-alignment_df %>%
    pivot_longer(2:ncol(.)) %>% 
    group_by(pos_alignment) %>%
    mutate(n=length(unique(value))) %>%
    dplyr::filter(n>1) %>%
    pivot_wider()
  
  #summarize it all
  res1<-alignment_df_more %>%
    summarise(median_tree_length=median(rel_tree_length),
              mean_tree_length=mean(rel_tree_length),
              sd_tree_length=sd(rel_tree_length)) %>% 
    mutate(variable_sites=dim(pos_split_diff)[1],
           total_sites=dim(alignment_df)[1],
           prop_variable=variable_sites/total_sites)
  return(res1)
  
}





# bininda tree ####
bininda_TRNP1<-summarize_inform_sites(tree_loc="protein/other_protein_alignments/tree_30sp_for_coevol.txt", aln_loc="protein/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol_longer.best.pep.fas")

#where I'll collect all values
bininda_prot_list<-list()

prot_list<-list.files("protein/other_protein_alignments/prank_out/aln_forChecking", pattern = "best.pep.fas") 

for (i in prot_list){
  
  bininda_prot_list[[i]]<-summarize_inform_sites(tree_loc="protein/other_protein_alignments/tree_30sp_for_coevol.txt", aln_loc=paste0("protein/other_protein_alignments/prank_out/aln_forChecking/",i))
  
}
  

bininda_prot_df<-bind_rows(bininda_prot_list, .id="protein") %>%
  bind_rows(bininda_TRNP1 %>% mutate(protein="TRNP1")) %>%
  mutate(protein=gsub(".fa|.fa.best.pep.fas","",protein))

saveRDS(bininda_prot_df,"protein/other_protein_alignments/summaries/bininda_prot_alignment_info.rds")




# plot the alignments

pdf("protein/other_protein_alignments/alns.pdf", onefile = TRUE)
for (i in prot_list){
  
  i=gsub(".fa|.fa.best.pep.fas","",i)
  loc = paste0("protein/other_protein_alignments/prank_out/aln_forChecking/",i,".fa.best.pep.fas")

  p<-msaplot(basic.tree, loc, width=7, bg_line=FALSE) +
    theme(legend.position = "none", 
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent"), 
          legend.box.background = element_rect(fill = "transparent"))+
    ggtitle(paste(i))
  #do.call("grid.arrange", msa)  
  print(p)
}
dev.off()






