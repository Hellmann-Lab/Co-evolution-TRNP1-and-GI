library(ggtree)
library(tidyverse)


gather_coevol_trees<-function(omega_tree, ds_tree, brain_tree=NULL, GI_tree, BW_tree=NULL){
  
  omega_tree2<-read.tree(paste0(omega_tree))
  omega_estimates<-data.frame(tip.label=omega_tree2$tip.label)
  omega_estimates<-data.frame(tip.label=omega_tree2$tip.label)
  omega_estimates<-separate(omega_estimates, col="tip.label", into=c("species", "CI.low", "CI.high"), sep="_")
  omega_estimates$CI.low<-as.numeric(omega_estimates$CI.low)
  omega_estimates$CI.high<-as.numeric(omega_estimates$CI.high)
  omega_estimates<-omega_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(omega=mean(c(CI.low, CI.high)))
  
  ds_tree2<-read.tree(paste0(ds_tree))
  ds_estimates<-data.frame(tip.label=ds_tree2$tip.label)
  ds_estimates<-data.frame(tip.label=ds_tree2$tip.label)
  ds_estimates<-separate(ds_estimates, col="tip.label", into=c("species", "ds.CI.low", "ds.CI.high"), sep="_")
  ds_estimates$ds.CI.low<-as.numeric(ds_estimates$ds.CI.low)
  ds_estimates$ds.CI.high<-as.numeric(ds_estimates$ds.CI.high)
  ds_estimates<-ds_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(ds=mean(c(ds.CI.low, ds.CI.high)))
  
  GI_tree2<-read.tree(paste0(GI_tree))
  GI_estimates<-data.frame(tip.label=GI_tree2$tip.label)
  GI_estimates<-separate(GI_estimates, col="tip.label", into=c("species", "GI.CI.low", "GI.CI.high"), sep="_")
  GI_estimates$GI.CI.low<-as.numeric(GI_estimates$GI.CI.low)
  GI_estimates$GI.CI.high<-as.numeric(GI_estimates$GI.CI.high)
  GI_estimates<-GI_estimates %>% dplyr::rowwise() %>%
    dplyr::mutate(GI=mean(c(GI.CI.low, GI.CI.high)))
  
  omega_estimates<-omega_estimates %>% 
    inner_join(ds_estimates, by="species")  %>%
    inner_join(GI_estimates, by="species") 
  
  if (!is.null(brain_tree)& !is.null(BW_tree)){
    brain_tree2<-read.tree(paste0(brain_tree))
    brain_estimates<-data.frame(tip.label=brain_tree2$tip.label)
    brain_estimates<-separate(brain_estimates, col="tip.label", into=c("species", "brain.CI.low", "brain.CI.high"), sep="_")
    brain_estimates$brain.CI.low<-as.numeric(brain_estimates$brain.CI.low)
    brain_estimates$brain.CI.high<-as.numeric(brain_estimates$brain.CI.high)
    brain_estimates<-brain_estimates %>% dplyr::rowwise() %>%
      dplyr::mutate(brain=mean(c(brain.CI.low, brain.CI.high)))
    
    BW_tree2<-read.tree(paste0(BW_tree))
    BW_estimates<-data.frame(tip.label=BW_tree2$tip.label)
    BW_estimates<-separate(BW_estimates, col="tip.label", into=c("species", "BW.CI.low", "BW.CI.high"), sep="_")
    BW_estimates$BW.CI.low<-as.numeric(BW_estimates$BW.CI.low)
    BW_estimates$BW.CI.high<-as.numeric(BW_estimates$BW.CI.high)
    BW_estimates<-BW_estimates %>% dplyr::rowwise() %>%
      dplyr::mutate(BodyM=mean(c(BW.CI.low, BW.CI.high)))
    
    omega_estimates<-omega_estimates %>%
      inner_join(brain_estimates, by="species") %>%
      inner_join(BW_estimates, by="species")
  }
  
  return(omega_estimates)
}








coevol_estimates<-gather_coevol_trees(omega_tree = "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/exon1_body_mass_brain_mass_GI_dsomggc1.postmeanomega.tre",
                                      ds_tree= "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/exon1_body_mass_brain_mass_GI_dsomggc1.postmeanbranchsynrate.tre",
                                      brain_tree = "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean2.tre",
                                      GI_tree = "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean3.tre",
                                      BW_tree = "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/exon1_body_mass_brain_mass_GI_dsomggc1.postmean1.tre")

saveRDS(coevol_estimates, "results_TRNP1_betterFer/body_mass_brain_mass_GI_30species/coevol_estimates.rds")



#also collect the GI only tree
coevol_estimates_GI<-gather_coevol_trees(omega_tree = "results_TRNP1_betterFer/GI_30species/exon1_GI_dsomggc1.postmeanomega.tre",
                                         ds_tree= "results_TRNP1_betterFer/GI_30species/exon1_GI_dsomggc1.postmeanbranchsynrate.tre",
                                         GI_tree = "results_TRNP1_betterFer/GI_30species/exon1_GI_dsomggc1.postmean1.tre")

saveRDS(coevol_estimates_GI, "results_TRNP1_betterFer/GI_30species/coevol_estimates.rds")




#also collect the brain only tree
coevol_estimates_brain<-gather_coevol_trees(omega_tree = "results_TRNP1_betterFer/brain_mass_30species/exon1_GI_dsomggc1.postmeanomega.tre",
                                            ds_tree= "results_TRNP1_betterFer/brain_mass_30species/exon1_GI_dsomggc1.postmeanbranchsynrate.tre",
                                            GI_tree = "results_TRNP1_betterFer/brain_mass_30species/exon1_GI_dsomggc1.postmean1.tre") %>%
  dplyr::rename(brain.CI.low=GI.CI.low, brain.CI.high=GI.CI.high)

saveRDS(coevol_estimates_brain, "results_TRNP1_betterFer/brain_mass_30species/coevol_estimates.rds")




#also collect the body only tree
coevol_estimates_body<-gather_coevol_trees(omega_tree = "results_TRNP1_betterFer/body_mass_30species/exon1_GI_dsomggc1.postmeanomega.tre",
                                            ds_tree= "results_TRNP1_betterFer/body_mass_30species/exon1_GI_dsomggc1.postmeanbranchsynrate.tre",
                                            GI_tree = "results_TRNP1_betterFer/body_mass_30species/exon1_GI_dsomggc1.postmean1.tre") %>%
  dplyr::rename(body.CI.low=GI.CI.low, body.CI.high=GI.CI.high)


saveRDS(coevol_estimates_body, "results_TRNP1_betterFer/body_mass_30species/coevol_estimates.rds")

