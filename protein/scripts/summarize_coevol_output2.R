
#load libraries 
libs<-c("tidyverse","cowplot","GenomicRanges","data.table","reshape2", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","grid","gtable","xtable","rr2")
sapply(libs, require, character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/")

#function
gather_coevol_trees<-function(omega_tree, brain_tree, GI_tree, BW_tree){
  
  omega_tree2<-read.tree(paste0(omega_tree))
  omega_estimates<-data.frame(tip.label=omega_tree2$tip.label)
  omega_estimates<-data.frame(tip.label=omega_tree2$tip.label)
  omega_estimates<-separate(omega_estimates, col="tip.label", into=c("species", "CI.low", "CI.high"), sep="_")
  omega_estimates$CI.low<-as.numeric(omega_estimates$CI.low)
  omega_estimates$CI.high<-as.numeric(omega_estimates$CI.high)
  omega_estimates<-omega_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(omega=mean(c(CI.low, CI.high)))
  
  brain_tree2<-read.tree(paste0(brain_tree))
  brain_estimates<-data.frame(tip.label=brain_tree2$tip.label)
  brain_estimates<-separate(brain_estimates, col="tip.label", into=c("species", "brain.CI.low", "brain.CI.high"), sep="_")
  brain_estimates$brain.CI.low<-as.numeric(brain_estimates$brain.CI.low)
  brain_estimates$brain.CI.high<-as.numeric(brain_estimates$brain.CI.high)
  brain_estimates<-brain_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(brain=mean(c(brain.CI.low, brain.CI.high)))
  
  GI_tree2<-read.tree(paste0(GI_tree))
  GI_estimates<-data.frame(tip.label=GI_tree2$tip.label)
  GI_estimates<-separate(GI_estimates, col="tip.label", into=c("species", "GI.CI.low", "GI.CI.high"), sep="_")
  GI_estimates$GI.CI.low<-as.numeric(GI_estimates$GI.CI.low)
  GI_estimates$GI.CI.high<-as.numeric(GI_estimates$GI.CI.high)
  GI_estimates<-GI_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(GI=mean(c(GI.CI.low, GI.CI.high)))
  
  BW_tree2<-read.tree(paste0(BW_tree))
  BW_estimates<-data.frame(tip.label=BW_tree2$tip.label)
  BW_estimates<-separate(BW_estimates, col="tip.label", into=c("species", "BW.CI.low", "BW.CI.high"), sep="_")
  BW_estimates$BW.CI.low<-as.numeric(BW_estimates$BW.CI.low)
  BW_estimates$BW.CI.high<-as.numeric(BW_estimates$BW.CI.high)
  BW_estimates<-BW_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(BodyM=mean(c(BW.CI.low, BW.CI.high)))
  
  omega_estimates<-omega_estimates %>% 
    inner_join(brain_estimates, by="species") %>%
    inner_join(GI_estimates, by="species") %>%
    inner_join(BW_estimates, by="species")
  
  return(omega_estimates)
}




coevol_estimates<-gather_coevol_trees(omega_tree = "protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.postmeanomega.tre",
                                      brain_tree = "protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean2.tre",
                                      GI_tree = "protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean3.tre",
                                      BW_tree = "protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean1.tre")



#annotate
tree.exon1.coding.coevol<-read.tree("protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt") 
ggtree(tree.exon1.coding.coevol)+geom_tiplab()+xlim(0,200) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

primates_coevol_res<-extract.clade(tree.exon1.coding.coevol, node=38)
rodents_coevol_res<-extract.clade(tree.exon1.coding.coevol, node=36)
carnivores_coevol_res<-extract.clade(tree.exon1.coding.coevol, node=58)

coevol_estimates$clade<-"Other"
coevol_estimates$clade[coevol_estimates$species %in% primates_coevol_res$tip.label]<-"Primate"
coevol_estimates$clade[coevol_estimates$species %in% rodents_coevol_res$tip.label]<-"Rodent"
coevol_estimates$clade[coevol_estimates$species %in% carnivores_coevol_res$tip.label]<-"Carnivore"

coevol_estimates$interesting_species<-NA
coevol_estimates$interesting_species[coevol_estimates$species=="Homosapien"]<-"human"
coevol_estimates$interesting_species[coevol_estimates$species=="Tursiopstr"]<-"dolphin"
coevol_estimates$interesting_species[coevol_estimates$species=="Musmusculu"]<-"mouse"
coevol_estimates$interesting_species[coevol_estimates$species=="Mustelaput"]<-"ferret"
coevol_estimates$interesting_species[coevol_estimates$species=="Macacamula"]<-"macaque"
coevol_estimates$interesting_species[coevol_estimates$species=="Otolemurga"]<-"galago"

coevol_estimates$experiment<-"no"
coevol_estimates$experiment[is.na(coevol_estimates$interesting_species)==F]<-"yes"

saveRDS(coevol_estimates, "protein/coevol/results/for_figures/coevol_3phenos_31sp_summarized.rds")
