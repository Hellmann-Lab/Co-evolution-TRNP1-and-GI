#analyse the proliferation rates induced by the different TRNP1 orthologues in vitro (mouse E14 neural stem cells)

#load libraries  
libs<-c("tidyverse", "broom","multcomp","ape", "readr","geiger","nlme","phytools","xtable")
sapply(libs, require, character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/protein/data/proliferation/")



################################################################
# Collect the proliferation data from different experiments ####
################################################################


#data from the first experiment
prolif_rate_Trnp1_exp1 <-read_delim("input_data/goetz_prolif_cells_full.csv", 
                                    ";", escape_double = FALSE, 
                                    locale = locale(decimal_mark = ","), 
                                    trim_ws = TRUE)

prolif_rate_Trnp1_exp1<-prolif_rate_Trnp1_exp1 %>%
  mutate(species=case_when(plasmid=="msTRNP" ~ "mouse",
                           plasmid=="hmTRNP" ~ "human",
                           plasmid=="ferTRNP" ~ "ferret",
                           plasmid=="GFP" ~ "GFP"))


#sum up for each biological replicate & orthologue since these were not separate wells according to miriam
prolif_rate_Trnp1_exp1<-prolif_rate_Trnp1_exp1 %>%
  dplyr::group_by(n, species) %>%
  dplyr::summarise(GFP_pos=sum(GFP), GFP_Ki67_pos=sum(GFP_Ki67)) %>%
  mutate(perc_prolif=GFP_Ki67_pos/GFP_pos)




#data from the second experiment
prolif_rate_Trnp1_exp2 <-read_delim("input_data/prolif_cells_full_repeated.csv", 
                                             ";", escape_double = FALSE, 
                                             locale = locale(decimal_mark = ","), 
                                             trim_ws = TRUE)

#exclude one batch because it was reported to behave weirdly during cell culturing (a lot of dead, unhealthy cells)
prolif_rate_Trnp1_exp2<-prolif_rate_Trnp1_exp2 %>%
  dplyr::rename(n=X1) %>%
  #exclude batch 3 because it was reported to behave weirdly during cell culturing (a lot of dead, unhealthy cells)
  filter(!n=="N3") %>%
  #rename the rest of the batches to be able to combine results with the first experiment
  mutate(n=case_when(n=="N1" ~ 6,
                      n=="N2" ~ 7,
                      n=="N4" ~ 8,
                      n=="N5" ~ 9,
                      n=="N6" ~ 10))


#data from the third experiment
prolif_rate_Trnp1_exp3 <-read_delim("input_data/repl_6_7.csv", 
                                    ",", escape_double = FALSE, 
                                    locale = locale(decimal_mark = ","), 
                                    trim_ws = TRUE)

prolif_rate_Trnp1_exp3<-prolif_rate_Trnp1_exp3 %>%
  dplyr::rename(n=X1) %>%
  mutate(n=case_when(n=="N6" ~ 11,
                     n=="N7" ~ 12))


prolif_rate_Trnp1_full_repeated<-prolif_rate_Trnp1_exp2 %>%
  bind_rows(prolif_rate_Trnp1_exp3) %>%
  tidyr::gather(species, no_of_cells, 3:9) %>%
  #exclude BrdU measurements since we do not use them in the downstream analysis
  filter(X2!="number of GFP+ BrdU+ cells") %>%
  #shorten the other names
  mutate(X2=case_when(X2=="number GFP + cells" ~ "GFP_pos",
                      X2=="number of GFP+ ki67+ cells" ~ "GFP_Ki67_pos")) %>%
  tidyr::spread(X2, no_of_cells) %>%
  mutate(perc_prolif=GFP_Ki67_pos/GFP_pos) %>%
  rowwise() %>%
  mutate_at("species", .funs=tolower) %>%
  mutate(species=gsub("gfp","GFP",species),
         species=gsub("macaca","macaque",species))



#combine all data
combined_prolif_all<-bind_rows(prolif_rate_Trnp1_exp1, prolif_rate_Trnp1_full_repeated) %>%
  mutate(Trnp1=case_when(species=="GFP"~ "no",
                         T ~ "yes"))
write.csv2(combined_prolif_all, "input_data/prolif_data_combined.csv")







################################################################
# Logistic regression to estimate proliferation rates ####
################################################################

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
print(xtable(anova_3mod_prolif_orth,digits=c(1,0,2,0,2,-1)), include.rownames=FALSE, file="output_xtables/prolif_mod_sel1.txt")





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


saveRDS(combined_glm_all, "proliferation_LR_res_orthologues.rds")

#prepare the table for the supplementary material
combined_glm_pretty<-combined_glm_all %>%
  dplyr::mutate(prolif_prob=round(prolif_prob,2), prolif_stderr=round(prolif_stderr,3)) %>%
  dplyr::select(term, prolif_prob, prolif_stderr) %>%
  dplyr::rename(Species=term, `Proliferation rate` = prolif_prob, `Proliferation SE`=prolif_stderr)

print(xtable(combined_glm_pretty, display=rep("s",ncol(combined_glm_pretty)+1)), include.rownames=FALSE, file="output_xtables/prolif_prob1.txt")



#get the pairwise comparisons
#info for general linear hypothesis https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
#it's a comparison of means: kind of like Tukey test (but assuming z-distribution), where YA is the larger of the two means being compared, YB is the smaller of the two means being compared, and SE is the standard error of the sum of the means (summed variances / sqrt sum of ns?)

x<-glht(mod_combined_glm_all,  linfct=mcp(species=c("human - mouse = 0",
                                                    "dolphin - human = 0",
                                                    "human - macaque = 0",
                                                    "human - galago = 0")))   

pq<-summary(x)$test

mtests <- cbind(pq$coefficients, pq$sigma, pq$tstat, pq$pvalues)
error <- attr(pq$pvalues, "error")
pname <- switch(x$alternativ, 
                less = paste("Pr(<", ifelse(x$df ==0, "z", "t"), ")", sep = ""), 
                greater = paste("Pr(>", ifelse(x$df == 0, "z", "t"), ")", sep = ""), 
                two.sided = paste("Pr(>|", ifelse(x$df == 0, "z", "t"), "|)", sep = ""))                                                                   
colnames(mtests) <- c("Estimate", "Std. Error", ifelse(x$df ==0, "z value", "t value"), pname)
xtable(mtests, digits=4)
print(xtable(mtests,digits=c(1,2,3,3,4)), file="output_xtables/prolif_comparison1.txt")






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
print(xtable(anova_3mod_prolif_Trnp1,digits=c(1,0,2,0,2,-1)), include.rownames=FALSE, file="output_xtables/prolif_mod_sel2.txt")



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

saveRDS(combined_glm_Trnp1, "proliferation_LR_res_TRNP1.rds")



#prepare the table for the supplementary material
combined_glm_Trnp1_pretty<-combined_glm_Trnp1 %>%
  dplyr::mutate(prolif_prob=round(prolif_prob,2), prolif_stderr=round(prolif_stderr,3)) %>%
  dplyr::select(term, prolif_prob, prolif_stderr) %>%
  dplyr::rename("TRNP1 present"=term, `Proliferation rate` = prolif_prob, `Proliferation SE`=prolif_stderr)

print(xtable(combined_glm_Trnp1_pretty, display=rep("s",ncol(combined_glm_pretty)+1)), include.rownames=FALSE, file="output_xtables/prolif_prob2.txt")





# do multiple comparisons of means to get test statistics
comb_trnp1<-glm(perc_prolif ~ 0+Trnp1+n, weights = GFP_pos, 
                data = combined_prolif_all, family = "binomial")
summary(comb_trnp1)

x2<-glht(comb_trnp1, linfct=mcp(Trnp1=c("yes - no = 0")))
pq2<-summary(x2)$test

mtests2 <- cbind(pq2$coefficients, pq2$sigma, pq2$tstat, pq2$pvalues)
error2 <- attr(pq2$pvalues, "error")
pname2 <- switch(x2$alternativ, 
                less = paste("Pr(<", ifelse(x2$df ==0, "z", "t"), ")", sep = ""), 
                greater = paste("Pr(>", ifelse(x2$df == 0, "z", "t"), ")", sep = ""), 
                two.sided = paste("Pr(>|", ifelse(x2$df == 0, "z", "t"), "|)", sep = ""))                                                                   
colnames(mtests2) <- c("Estimate", "Std. Error", ifelse(x2$df ==0, "z value", "t value"), pname)
xtable(mtests2, digits=3)

#p-value rounds up to a zero -- > use the one from summary output
summary(x2)
#2e-16
mtests2[4]<-2E-16
print(xtable(mtests2,digits=c(1,2,3,3,-1)), file="output_xtables/prolif_comparison2.txt")






################################################################
# PGLS: GI vs proliferation ####
################################################################

#compare to GI
pheno_data_new<-readRDS("/data/share/htp/TRNP1/paper_data/data_tables/pheno_data/pheno_data_new.rds")

combined_glm_all<-combined_glm_all %>%
  mutate(species=case_when(term=="macaque" ~ "Macaca_mulatta",
                           term=="galago" ~ "Otolemur_garnettii",
                           term=="ferret" ~ "Mustela_putorius",
                           term=="mouse" ~ "Mus_musculus",
                           term=="human" ~ "Homo_sapiens",
                           term=="dolphin" ~ "Tursiops_truncatus"))

prolif_vs_dnds<-combined_glm_all %>% left_join(pheno_data_new)
saveRDS(prolif_vs_dnds, "prolif_vs_dnds.rds")



#tree
tree.exon1.coding.29sp.full<-read.tree("../../trees/tree_TRNP1_coding_31sp.txt") 

tree.prolif<-drop.tip(tree.exon1.coding.29sp.full, tree.exon1.coding.29sp.full$tip.label[!tree.exon1.coding.29sp.full$tip.label %in% prolif_vs_dnds$species ])
rownames(prolif_vs_dnds)<-prolif_vs_dnds$species

mod_GI<-gls(log2(GI) ~ log2(prolif_prob), data=prolif_vs_dnds, correlation=corBrownian(1, tree.prolif,form=~species), method="ML")
summary(mod_GI)
qqnorm(mod_GI$residuals, pch = 1, frame = FALSE)
qqline(mod_GI$residuals, col = "steelblue", lwd = 2)
plot(mod_GI$fitted, mod_GI$residuals)

mod_GI_ctrl<-gls(log2(GI) ~ 1, data=prolif_vs_dnds, correlation=corBrownian(1, tree.prolif,form=~species), method="ML")
anova(mod_GI,mod_GI_ctrl)


anova_GI_reg<- as.data.frame(anova(mod_GI_ctrl,mod_GI)) 
anova_GI_reg$call<-c("log2(GI) ~ 1","log2(GI) ~ log2(Proliferation rate)")

print(xtable(anova_GI_reg,digits=c(1,1,1,1,2,2,2,1,2,3)),
      include.rownames=FALSE, file="output_xtables/GI_prolif_mod_sel_BM.tex")


mod_GI_prolif<-wrap_summary_table_BM(mod_GI) %>%
  mutate(Predictor=gsub("prolif_prob","Proliferation rate",Predictor),
         R2=round(R2.resid(mod_GI,mod_GI_ctrl),digits=3))
mod_GI_prolif$logLik[1]<-NA
mod_GI_prolif$R2[1]<-NA

print(xtable(mod_GI_prolif,digits=c(1,1,2,2,2,3,2,3)),
      include.rownames=FALSE, file="output_xtables/GI_prolif_mod_BM.tex")
