library(Biostrings)
library(tidyverse)
library(nlme)
library(ape)
library(geiger)


# PERMUTATIONS
# SHUFFLE TILES TO TEST FOR INTRON SIGNIFICANCE ####

#plot only the ones with GI and brain weight present +human1 activity --> 45 species btw
pheno_data<-read_rds("pheno_data/pheno_data.rds")
mammaltree<-read.tree("protein/trees/mammaltree.txt") 


#data
activity_summary_hum1_withPheno<-readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds") %>%
  filter(cell_line=="human1") %>%
  left_join(pheno_data[,c("brain_mass","GI","species")],
            by=c("species")) %>%
  filter(!is.na(GI) & !is.na(brain_mass))
                                         
                                    
cath_intron<-activity_summary_hum1_withPheno %>%
  filter(region == "intron", primate_clade %in% c("Old World monkey", "Great ape")) 


#select shuffles
shuffles_for_intron<-activity_summary_hum1_withPheno %>%
  filter(! (region == "intron" & primate_clade %in% c("Old World monkey", "Great ape"))) 

#prepare for PGLS
intron_tree_apes_OWMs<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% cath_intron$species[cath_intron$primate_clade %in% c("Great ape","Old World monkey")]])



#shuffle 10 sequences 1000 times, do PGLS
N<-1000
permutation_out<-data.frame(n=integer(), L.Ratio=numeric(), p.value=numeric())
set.seed(123)

for (i in 1:N){
  
  shuffled<-data.frame(log2_total_activity = shuffles_for_intron %>% 
                         sample_n(10) %>% 
                         pull(log2_total_activity),
                       GI=cath_intron$GI, species=cath_intron$species)
  rownames(shuffled)<-shuffled$species
  
  mod1<-gls(log2(GI)~log2_total_activity,
            data=shuffled,
            correlation=corBrownian(value=1,phy=intron_tree_apes_OWMs, form=~species), 
            method="ML")
  mod0<-gls(log2(GI)~1, 
            data=shuffled,
            correlation=corBrownian(value=1,phy=intron_tree_apes_OWMs, form=~species), 
            method="ML")
  
  permutation_out<-permutation_out %>%
    bind_rows(data.frame(n=i,
                         L.Ratio=anova(mod1,mod0)[2,8],
                         p.value=anova(mod1,mod0)[2,9]))
}

table(permutation_out$p.value>=0.003)
8/1000