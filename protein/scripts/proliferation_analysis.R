#analyse the proliferation rates induced by the different TRNP1 orthologues in vitro (mouse E14 neural stem cells)

#load libraries  
libs<-c("tidyverse", "broom","multcomp","ape", "readr","geiger","nlme","phytools","xtable")
sapply(libs, require, character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/Co-evolution-TRNP1-and-GI/")



################################################################
# Logistic regression to estimate proliferation rates ####
################################################################

combined_prolif_all<-read.csv2("protein/data/proliferation/input_data/prolif_data_combined.csv")

# generate a set excluding GFP
combined_prolif_Trnp1<- combined_prolif_all %>%
  filter(species!="GFP")


#compare model fits
#orthologue effects####
combined_prolif_Trnp1$species<-factor(combined_prolif_Trnp1$species, levels=c("macaque", "galago","ferret","mouse","human","dolphin"))
combined_prolif_Trnp1$n<-as.factor(combined_prolif_Trnp1$n)


mod_null_orth<-glm(perc_prolif ~ 1, weights = GFP_pos, 
                          data = combined_prolif_Trnp1, family = "binomial")

mod_species_orth<-glm(perc_prolif ~ species, weights = GFP_pos, 
                          data = combined_prolif_Trnp1, family = "binomial")

mod_species_n_orth<-glm(perc_prolif ~ species+n, weights = GFP_pos, 
                             data = combined_prolif_Trnp1, family = "binomial")

anova_3mod_prolif_orth<-anova(mod_null_orth, mod_species_orth, mod_species_n_orth, test="LRT") # the best model is the full one
print(xtable(anova_3mod_prolif_orth,digits=c(1,0,2,0,2,-1)), include.rownames=FALSE, file="protein/data/proliferation/output_xtables/prolif_mod_sel1.txt")





#get the estimates from the best -full- model (set intercept to zero to get the absolute proliferation rates)
mod_combined_glm_all<-glm(perc_prolif ~ 0+species+n, weights = GFP_pos, 
                           data = combined_prolif_Trnp1, family = "binomial")

combined_glm_all<-tidy(glm(perc_prolif ~ 0+species+n, weights = GFP_pos, 
                           data = combined_prolif_Trnp1, family = "binomial"))

combined_glm_all$term<-gsub("species","",combined_glm_all$term)
combined_glm_all<-combined_glm_all %>% filter(!grepl("n[1-9|10-12]", term))
combined_glm_all$term<-factor(combined_glm_all$term, levels=c("macaque", "galago","ferret","mouse","human","dolphin"))

combined_glm_all<-combined_glm_all %>%
  mutate(clade=case_when(term %in% c("human","galago","macaque") ~ "Primate",
                         term=="mouse" ~ "Rodent",
                         term=="dolphin" ~ "Cetacean",
                         term=="ferret" ~ "Carnivore")) %>%
  rowwise() %>%
  mutate(prolif_prob=exp(estimate)/(1+exp(estimate)),
         prolif_stderr_max=exp(estimate+std.error)/(1+exp(estimate+std.error)),
         prolif_stderr_min=exp(estimate-std.error)/(1+exp(estimate-std.error)),
         prolif_stderr=prolif_stderr_max-prolif_prob)


saveRDS(combined_glm_all, "protein/data/proliferation/proliferation_LR_res_orthologues.rds")

#prepare the table for the supplementary material
combined_glm_pretty<-combined_glm_all %>%
  dplyr::mutate(prolif_prob=round(prolif_prob,2), prolif_stderr=round(prolif_stderr,3)) %>%
  dplyr::select(term, prolif_prob, prolif_stderr) %>%
  dplyr::rename(Species=term, `Proliferation rate` = prolif_prob, `Proliferation SE`=prolif_stderr)

print(xtable(combined_glm_pretty, display=rep("s",ncol(combined_glm_pretty)+1)), include.rownames=FALSE, file="protein/data/proliferation/output_xtables/prolif_prob1.txt")



#get the pairwise comparisons
#info for general linear hypothesis https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
#it's a comparison of means: kind of like Tukey test (but assuming z-distribution), where YA is the larger of the two means being compared, YB is the smaller of the two means being compared, and SE is the standard error of the sum of the means (summed variances / sqrt sum of ns?)

wrap_glht_output<-function(glht_output){
  
  pq<-summary(glht_output)$test
  mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
  error <- attr(pq$pvalues, "error")
  pname <- switch(glht_output$alternativ, 
                  less = paste("Pr(<", ifelse(glht_output$df ==0, "z", "t"), ")", sep = ""), 
                  greater = paste("Pr(>", ifelse(glht_output$df == 0, "z", "t"), ")", sep = ""), 
                  two.sided = paste("Pr(>|", ifelse(glht_output$df == 0, "z", "t"), "|)", sep = ""))     
  colnames(mtests) <- c("Estimate", "Std. Error", ifelse(glht_output$df ==0, "z value", "t value"), pname)
  return(mtests)
}


x<-glht(mod_combined_glm_all,  linfct=mcp(species=c("human - mouse = 0",
                                                    "dolphin - human = 0",
                                                    "human - macaque = 0",
                                                    "human - galago = 0")))   

mtests1<-wrap_glht_output(x)
xtable(mtests1, digits=4)
print(xtable(mtests1,digits=c(1,2,3,3,4)), file="protein/data/proliferation/output_xtables/prolif_comparison1.txt")






#general TRNP1 effect####
#general effect of TRNP1 presence on the proliferation (irrespectively of the species)
combined_prolif_all$Trnp1<-as.factor(combined_prolif_all$Trnp1)
combined_prolif_all$n<-as.factor(combined_prolif_all$n)

mod_null_Trnp1<-glm(perc_prolif ~ 1, weights = GFP_pos, 
                   data = combined_prolif_all, family = "binomial")

mod_Trnp1_pres<-glm(perc_prolif ~ Trnp1, weights = GFP_pos, 
                      data = combined_prolif_all, family = "binomial")

mod_Trnp1_pres_n<-glm(perc_prolif ~ Trnp1+n, weights = GFP_pos, 
                        data = combined_prolif_all, family = "binomial")

anova_3mod_prolif_Trnp1<-anova(mod_null_Trnp1, mod_Trnp1_pres, mod_Trnp1_pres_n, test="LRT") # the best model is the full one
print(xtable(anova_3mod_prolif_Trnp1,digits=c(1,0,2,0,2,-1)), include.rownames=FALSE, file="protein/data/proliferation/output_xtables/prolif_mod_sel2.txt")


#best model: Trnp1 +n; set intercept to 0 to backcalculate the absolute prolif rates
combined_glm_Trnp1<-tidy(glm(perc_prolif ~ 0+Trnp1+n, weights = GFP_pos, 
                             data = combined_prolif_all, family = "binomial"))

combined_glm_Trnp1$term<-gsub("Trnp1","",combined_glm_Trnp1$term)
combined_glm_Trnp1<-combined_glm_Trnp1 %>% filter(!grepl("n[1-9|10-12]", term))


#backcalculate proliferation probabilities
combined_glm_Trnp1<-combined_glm_Trnp1 %>%
  rowwise() %>%
  mutate(prolif_prob=exp(estimate)/(1+exp(estimate)),
         prolif_stderr_max=exp(estimate+std.error)/(1+exp(estimate+std.error)),
         prolif_stderr_min=exp(estimate-std.error)/(1+exp(estimate-std.error)),
         prolif_stderr=prolif_stderr_max-prolif_prob)

saveRDS(combined_glm_Trnp1, "protein/data/proliferation/proliferation_LR_res_TRNP1.rds")



#prepare the table for the supplementary material
combined_glm_Trnp1_pretty<-combined_glm_Trnp1 %>%
  dplyr::mutate(prolif_prob=round(prolif_prob,2), prolif_stderr=round(prolif_stderr,3)) %>%
  dplyr::select(term, prolif_prob, prolif_stderr) %>%
  dplyr::rename("TRNP1 present"=term, `Proliferation rate` = prolif_prob, `Proliferation SE`=prolif_stderr)

print(xtable(combined_glm_Trnp1_pretty, display=rep("s",ncol(combined_glm_pretty)+1)), include.rownames=FALSE, file="protein/data/proliferation/output_xtables/prolif_prob2.txt")





# do multiple comparisons of means to get test statistics on the differences between conditions
comb_trnp1<-glm(perc_prolif ~ 0+Trnp1+n, weights = GFP_pos, 
                data = combined_prolif_all, family = "binomial")
summary(comb_trnp1)

x2<-glht(comb_trnp1, linfct=mcp(Trnp1=c("yes - no = 0")))
mtests2<-wrap_glht_output(x2)
xtable(mtests2, digits=3)

#p-value rounds up to a zero -- > use the one from summary output
summary(x2)
#2e-16
mtests2[4]<-2E-16
print(xtable(mtests2,digits=c(1,2,3,3,-1)), file="protein/data/proliferation/output_xtables/prolif_comparison2.txt")


 





################################################################
# PGLS: GI vs proliferation ####
################################################################

wrap_summary_table_BM<-function(model_output){
  repl_summary<-summary(model_output)
  repl_summary_table<-as.data.frame(repl_summary$tTable) %>%
    rownames_to_column(var="Predictor")%>%
    dplyr::rename(p.value="p-value", t.value="t-value") %>%
    mutate(Value=round(Value,3),
           p.value=round(p.value, 4),
           t.value=round(t.value, 4),
           logLik=model_output$logLik)
  return(repl_summary_table)
}

#compare to GI
pheno_data<-readRDS("pheno_data/pheno_data.rds")

prolif_vs_GI<-combined_glm_all %>% 
  mutate(species=case_when(term=="macaque" ~ "Macaca_mulatta",
                           term=="galago" ~ "Otolemur_garnettii",
                           term=="ferret" ~ "Mustela_putorius",
                           term=="mouse" ~ "Mus_musculus",
                           term=="human" ~ "Homo_sapiens",
                           term=="dolphin" ~ "Tursiops_truncatus")) %>%
  left_join(pheno_data)
saveRDS(prolif_vs_GI, "protein/data/proliferation/prolif_vs_GI.rds")



#tree
tree.coding31<-read.tree("protein/trees/tree_TRNP1_coding_31sp.txt") 

tree.prolif<-drop.tip(tree.coding31, tree.coding31$tip.label[!tree.coding31$tip.label %in% prolif_vs_GI$species ])
rownames(prolif_vs_GI)<-prolif_vs_GI$species

mod_GI<-gls(log2(GI) ~ log2(prolif_prob), data=prolif_vs_GI, correlation=corBrownian(1, tree.prolif,form=~species), method="ML")
summary(mod_GI)
#qq plots
qqnorm(mod_GI$residuals, pch = 1, frame = FALSE)
qqline(mod_GI$residuals, col = "steelblue", lwd = 2)
plot(mod_GI$fitted, mod_GI$residuals)

mod_GI_ctrl<-gls(log2(GI) ~ 1, data=prolif_vs_GI, correlation=corBrownian(1, tree.prolif,form=~species), method="ML")
anova(mod_GI,mod_GI_ctrl)


anova_GI_reg<- as.data.frame(anova(mod_GI_ctrl,mod_GI)) 
anova_GI_reg$call<-c("log2(GI) ~ 1","log2(GI) ~ log2(Proliferation rate)")

print(xtable(anova_GI_reg,digits=c(1,1,1,1,2,2,2,1,2,3)),
      include.rownames=FALSE, file="protein/data/proliferation/output_xtables/GI_prolif_mod_sel_BM.tex")


mod_GI_prolif<-wrap_summary_table_BM(mod_GI) %>%
  mutate(Predictor=gsub("prolif_prob","Proliferation rate",Predictor),
         R2=round(R2.resid(mod_GI,mod_GI_ctrl),digits=3))
mod_GI_prolif$logLik[1]<-NA
mod_GI_prolif$R2[1]<-NA

print(xtable(mod_GI_prolif,digits=c(1,1,2,2,2,3,2,3)),
      include.rownames=FALSE, file="protein/data/proliferation/output_xtables/GI_prolif_mod_BM.tex")




#exclude mouse####
prolif_vs_GI_nomouse<-prolif_vs_GI %>%
  filter(term!="mouse")

tree.prolif.nomouse<-drop.tip(tree.coding31, tree.coding31$tip.label[!tree.coding31$tip.label %in% prolif_vs_GI_nomouse$species ])
rownames(prolif_vs_GI_nomouse)<-prolif_vs_GI_nomouse$species


mod_GI_nomouse<-gls(log2(GI) ~ log2(prolif_prob), data=prolif_vs_GI_nomouse, correlation=corBrownian(1, tree.prolif.nomouse,form=~species), method="ML")
summary(mod_GI_nomouse)

mod_GI_ctrl_nomouse<-gls(log2(GI) ~ 1, data=prolif_vs_GI_nomouse, correlation=corBrownian(1, tree.prolif.nomouse,form=~species), method="ML")
anova(mod_GI_nomouse,mod_GI_ctrl_nomouse)



