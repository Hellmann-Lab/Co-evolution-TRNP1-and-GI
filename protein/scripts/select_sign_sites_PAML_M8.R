# first, compare site model with a fraction of proteins with omega>1 with a null model
setwd("/data/share/htp/TRNP1/paper_data/")

system("bash protein/scripts/summarise_PAML_out.sh")

m8_vs_m7_trnp1<-read.table("protein/PAML/M8_vs_M7.txt", sep="\t") %>%
  distinct() %>%
  mutate(lnL=word(V3,7,7, sep=" "),
         p0=word(V4,7,7, sep=" ")) %>%
  dplyr::rename(protein=V1, model=V2) %>%
  dplyr::select(protein, model, lnL, p0) %>%
  tidyr::pivot_wider(names_from=model, values_from=c(lnL, p0)) %>%
  rowwise() %>%
  ##compare models using likelihood ratio test (chi-sq distr, 2 dfs of freedom)
  mutate(LRT=as.numeric(lnL_NS8)-as.numeric(lnL_NS7),
         pval=1-pchisq(2*LRT, df=2))

saveRDS(m8_vs_m7_trnp1,"protein/PAML/M8_vs_M7_LTR_pvals.rds")


#compare models using likelihood ratio test (chi-sq distr, 2 dfs of freedom)
# LRT= lnL(M8)-lnL(M7). 
# Compare 2*LRT to Χ2 to calculate p-value (df=2).
# 1 - pchisq(2*2.08554, df = 2)  
# 1 - pchisq(2*(-5438.53-(-5447.20)), df = 2)  





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



