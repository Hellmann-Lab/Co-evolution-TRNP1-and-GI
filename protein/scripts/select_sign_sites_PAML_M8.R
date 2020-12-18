#selecting the sites under positive selection according to Naive Empirical Bayes analysis (M8)

#Positively selected sites (*: P>95%; **: P>99%)
#(amino acids refer to 1st sequence: Mus_pahari)
library(tidyverse)
library(xtable)

sites_under_sel<-read.table("protein/PAML/NS8_NEB_sites")
colnames(sites_under_sel)<-c("Alignment_position","AA_in_Mpahari", "Pr_w>1", "post_mean_w")

sites_under_sel_pp95<-sites_under_sel %>%
  mutate(`Pr_w_2`=gsub("[*]|[**]","", `Pr_w>1`)) %>%
  filter(Pr_w_2>0.95) %>%
  dplyr::select(Alignment_position, `Pr_w>1`, post_mean_w) %>%
  dplyr::rename("Alignment position"=Alignment_position,
                "Post mean omega"=post_mean_w)

print(xtable(sites_under_sel_pp95), include.rownames=FALSE, file="protein/PAML/PAML_M8_NEB_sites_xtable.txt")
saveRDS(sites_under_sel_pp95,"protein/PAML/PAML_M8_NEB_sites_pp95.rds")
