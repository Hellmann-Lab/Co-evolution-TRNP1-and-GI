
library(tidyverse)
library(ggtree)
library(Biostrings)
library(readxl)
library(nlme)
library(ape)
library(geiger)
library(stringr)


# put all the results together
# CCDS IDs
#CCDS<-read_delim("protein/other_protein_alignments/CCDS/CCDS.current.txt", delim="\t")
trnp1_ccds<-"CCDS41289.1"


source("protein/scripts/summarize_cor_output_short.R")


#get convergence ####
list_convergence<-list.files("protein/other_protein_alignments/blat_322cand/tracecomp_combined_all/")

converg_out<-list()
for (i in list_convergence){
  for (run in c("1vs2","1vs3","2vs3")){
    converg_out[[paste0(i,"_",run)]]<-read.table(paste0("protein/other_protein_alignments/blat_322cand/tracecomp_combined_all/",i,"/tracecomp_body_mass_brain_mass_GI_run",run,".txt"), comment.char = "", skip = 2)
  }}

#add trnp1
for (run in c("1vs2","1vs3","2vs3")){
  converg_out[[paste0(trnp1_ccds,"_",run)]]<-read.table(paste0("protein/coevol/results/all/body_mass_brain_mass_GI_30species/tracecomp_summary/tracecomp_body_mass_brain_mass_GI_run",run,".txt"), skip=2)
}

converg_out_df<-bind_rows(converg_out, .id="protein") %>%
  tidyr::separate(protein, into = c("protein","run"), sep = "_")


#these should be excluded - did not converge
comb_unconv<-converg_out_df %>% filter(V3>0.3 | V2<50) %>%
  group_by(protein) %>%
  dplyr::summarise(n_unconverged_params=length(V1), unconverged_params=paste(V1, collapse=",")) 

#save convergence (y/n)
converg_out_df %>%
  distinct(protein) %>%
  mutate(converged=ifelse(protein %in% comb_unconv$protein, F,T)) %>%
  saveRDS(file = "protein/other_protein_alignments/summaries/convergence.rds")





# alignment info in AAs, and variation ####

bininda_prot_df<-readRDS("protein/other_protein_alignments/summaries/bininda_prot_alignment_info.rds") %>%
  mutate(protein=paste0(protein,".fa")) %>%
  mutate(protein=ifelse(protein=="TRNP1.fa","CCDS41289.1",protein)) 






# plot correlation coefficients #### 

all_corfiles<-list.files("protein/other_protein_alignments/coevol_cor_tables/", full.names = T)

outcors<-list()
for (i in all_corfiles){
  outcors[[i]]<-summarize_cor_output(i)
}


trnp1_corfiles<-list.files("protein/coevol/results/all/body_mass_brain_mass_GI_30species/", full.names = T, pattern = ".cov")
outcors_trnp1<-list()
for (i in trnp1_corfiles){
  outcors_trnp1[[i]]<-summarize_cor_output(i)
}

#combine and calculate averages 
outcors_df<-bind_rows(outcors, .id="name") %>%
  mutate(name2=word(name,10,10,"/"),
         protein=gsub(".fa.*","",name2),
         run=gsub(".*dsomggc|[.]cov","",name2)) %>%
  bind_rows(bind_rows(outcors_trnp1, .id="name") %>% 
              mutate(name2=word(name,13,13,"/"),
                     protein=trnp1_ccds,
                     run=gsub(".*dsomggc|[.]cov","",name2))) %>%
  filter(combi=="body_mass,omega" | combi=="brain_mass,omega" | combi=="GI,omega") %>%
  #exclude unconverged
  filter(! protein %in% gsub("[.]fa","",comb_unconv$protein)) %>%
  group_by(combi,protein) %>%
  dplyr::summarise(marginal_cor=mean(marginal_cor),
                   marginal_pp=mean(marginal_pp),
                   partial_cor=mean(partial_cor),
                   partial_pp=mean(partial_pp)) %>%
  group_by(combi) %>%
  mutate(mean_mar_cor=mean(marginal_cor),
         mean_par_cor=mean(partial_cor)) %>%
  ungroup() %>%
  mutate(combi=gsub(",omega","",combi)) %>%
  left_join(CCDS %>% dplyr::select(ccds_id, gene), by=c("protein"="ccds_id")) %>%
  distinct(gene, combi, .keep_all = T)



outcors_long<-bind_rows(outcors_df[,c(1:4,7)] %>% 
                          mutate(type="Marginal") %>% 
                          dplyr::rename(correlation=marginal_cor, pp = marginal_pp, mean=mean_mar_cor),
                        outcors_df[,c(1:2,5:6,8)] %>% 
                          mutate(type="Partial") %>%
                          dplyr::rename(correlation=partial_cor, pp = partial_pp, mean=mean_par_cor)) %>%
  mutate(combi2=case_when(combi=="brain_mass" ~ "Brain size",
                          combi=="body_mass" ~ "Body mass",
                          T ~ combi),
         TRNP1=ifelse(protein==trnp1_ccds,1,0)) %>%
  arrange(combi2,type, -correlation, -TRNP1) %>%
  group_by(combi2, type) %>%
  mutate(rank=1:length(type), total=length(combi)) %>%
  ungroup()

saveRDS(outcors_long,"protein/other_protein_alignments/summaries/cor_coefficients_allProt.rds")




#branch model output of PAML ####

br<-read.table("protein/other_protein_alignments/summaries/branch_model_out.txt", sep="\t", fill = T)
conv<-readRDS("protein/other_protein_alignments/summaries/convergence.rds") %>% mutate(protein=gsub(".fa","",protein))
colnames(br)<-br[1,]
br<-br[-1,]

br<-br %>% distinct() %>% mutate(dN=as.numeric(dN),
                                 dS=as.numeric(dS),
                                 `dN/dS`=dN/dS,
                                 Protein=gsub(".fa.best.nuc.phy","",Protein)) %>%
  #this one is duplicated in terms of 2 CCDS for 1 gene.. remove one 
  filter(Protein!="CCDS59451.1") %>%
  left_join(conv, by=c("Protein"="protein")) %>% 
  mutate(converged=ifelse(Protein=="TRNP1",T,converged),
         dS=as.numeric(dS))



# IDR content downloaded from D2P2 database
p.fa<-readDNAStringSet("protein/other_protein_alignments/CCDS_protein.current.322.fa")
names(p.fa)<-word(names(p.fa),1,1,"[|]")

ctrl_IDRs<-read_excel("protein/other_protein_alignments/summaries/d2p2_ctrlProt.xlsx") %>%
  fill(Protein) %>%
  left_join(data.frame(Protein = names(p.fa), width = width(p.fa))) %>%
  mutate(width_IDR = End-Start) %>%
  group_by(Protein) %>%
  summarise(width = width[1],
            width_IDR = sum(width_IDR),
            prop_IDR = width_IDR / width) %>%
  left_join(CCDS %>% dplyr::select(gene, ccds_id), by=c("Protein" = "ccds_id"))
