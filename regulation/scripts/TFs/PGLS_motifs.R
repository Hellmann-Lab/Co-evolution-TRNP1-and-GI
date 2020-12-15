#load libraries 
libs<-c("tidyverse","cowplot","data.table","reshape2", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","GGally","RColorBrewer","rr2","purrr")
sapply(libs, require, character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/")
source("regulation/scripts/MPRA/MPRA_helper_functions.R")


#############################
# LOAD THE RELEVANT DATA ####
#############################

# MPRA data 
activity_overlap_summary<-readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds")

#load the combined phenotype data
pheno_data<-readRDS("pheno_data/pheno_data.rds") %>%
  dplyr::select("species","EQ","GI") 

#load the upgraded mammal tree from Bininda-Emonds 2007
mammaltree<-read.tree("protein/trees/mammaltree.txt") 

#expressed TFs
NPC_expr_tf_df<-readRDS("regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/NPC_expressed_TF_IDs_jaspar2020.rds")
unique_motifs_vec<-unique(NPC_expr_tf_df$motif)

#ctrain weights
ctrain_weights<-read.table("regulation/data/TFs/ClusterBuster/results/OWMs_apes_intron_weighted/ctrain_output.txt", skip = 19, nrows = length(unique_motifs_vec), sep=",") %>%
  tidyr::separate(V1, into=c("score","motif"), sep="MA",extra="merge") %>%
  tidyr::separate(motif, into=c("motif","SYMBOL"), sep="\\s") %>%
  mutate(motif=paste0("MA",motif),
         score=as.numeric(score),
         SYMBOL_clean=toupper(SYMBOL),
         SYMBOL_clean=gsub("::|-","_",SYMBOL_clean),
         SYMBOL_clean=gsub("[(]VAR.2[)]","_VAR2",SYMBOL_clean))

intron_motifs_1_OWMs<-ctrain_weights %>%
  filter(score>=1)
saveRDS(intron_motifs_1_OWMs, "regulation/data/TFs/ClusterBuster/results/intron_motifs_1_OWMs.rds")


# cbust scores
scores_NPC_expressed_df<-gather_scores("OWMs_apes_intron_weighted")
scores_NPC_expressed_df$End-scores_NPC_expressed_df$Start #reasonable length
#all species have only one cluster


#####################
# PREPARE MATRIX ####
#####################

scores_top22TFs_long<-scores_NPC_expressed_df %>%
  gather(motif, score, grep("MA",names(scores_NPC_expressed_df))) %>%
  inner_join(intron_motifs_1_OWMs[,c("SYMBOL_clean","SYMBOL","motif")])

ggplot(scores_top22TFs_long, aes(x=score))+geom_density()+facet_wrap(~SYMBOL_clean,ncol=5, nrow=6, scales = "free")
ggplot(scores_top22TFs_long, aes(x=species,y=score))+geom_point()+coord_flip()+facet_wrap(~SYMBOL_clean,ncol=5, nrow=5, scales = "free")




###########
# PGLS ####
###########

#prepare the data and tree
scores_top22TFs_mat<-scores_top22TFs_long %>%
  dplyr::select(species, score, SYMBOL_clean) %>%
  spread(SYMBOL_clean, score) %>%
  left_join(scores_NPC_expressed_df[,c("species","Score")]) %>%
  left_join(activity_overlap_summary) %>%
  left_join(pheno_data[,c("GI","species")]) %>%
  filter(!is.na(GI) & region=="intron" & cell_line=="human1") %>%
  left_join(intron_activity_GI_hum1)

ggcorr(scores_top22TFs_mat[, colnames(scores_top22TFs_mat) %in% intron_motifs_1_OWMs$SYMBOL_clean],label=TRUE, label_round=2, label_size = 2, size=2, hjust=0.7)

tree.intron.tfs22<- drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% scores_top22TFs_mat$species])
rownames(scores_top22TFs_mat)<-scores_top22TFs_mat$species
name.check(tree.intron.tfs22,scores_top22TFs_mat)



#intron activity ####
#test the individual effects
activity_indSel_owms<-select_indiv_motifs("scores_top22TFs_mat",pheno="log2_total_activity", intron_motifs_1_OWMs$SYMBOL_clean, "tree.intron.tfs22", cutoff=0.05)
#ZBTB26,SOX8,CTCF
activity_all_lrts<-activity_indSel_owms[[3]]

TFs<-bind_rows(activity_all_lrts, .id="TF") %>% filter(`p-value`<0.05)
activity_sign_lrts_df<-bind_rows(activity_all_lrts, .id="TF") %>% filter(TF %in% TFs$TF)



#compare the combined model to the null model
mod_activity_3tf<-gls(log2_total_activity~CTCF+SOX8+ZBTB26, data=scores_top22TFs_mat, 
                             correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
scores_top22TFs_mat$fitted_intron_3TFs<-mod_activity_3tf$fitted
drop_combined<-as.data.frame(summary(mod_activity_3tf)$tTable) %>% rownames_to_column("motif")

mod_activity_null<-gls(log2_total_activity~1, data=scores_top22TFs_mat, 
                       correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
act_drop0<-anova(mod_activity_null,mod_activity_3tf) # full better
R2.resid(mod_activity_3tf,mod_activity_null) #0.78

#diagnostics
qqnorm(mod_activity_3tf$residuals, pch = 1, frame = FALSE)
qqline(mod_activity_3tf$residuals, col = "steelblue", lwd = 2)
plot(mod_activity_3tf$fitted, mod_activity_3tf$residuals)




# drop either TF
mod_activity_2tf<-gls(log2_total_activity~CTCF+ZBTB26, data=scores_top22TFs_mat, correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
act_drop1<-anova(mod_activity_2tf,mod_activity_3tf) # full better

mod_activity_2tf_2<-gls(log2_total_activity~SOX8+ZBTB26, data=scores_top22TFs_mat, correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
act_drop2<-anova(mod_activity_2tf_2,mod_activity_3tf)

mod_activity_2tf_3<-gls(log2_total_activity~CTCF+SOX8, data=scores_top22TFs_mat, correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
act_drop3<-anova(mod_activity_2tf_3,mod_activity_3tf) 
#the combination of the 3 is the best

lrt_activity_drop<-bind_rows(act_drop0,act_drop1,act_drop2,act_drop3) %>%
  mutate(call=gsub("(gls[(]model = )","",call),
         call=gsub("\\,.*","",call),
         call=gsub("log2_total_activity","log2(intron)",call))

print(xtable(lrt_activity_drop, digits=c(1,1,1,1,2,2,2,1,2,3)),include.rownames=FALSE,
      hline.after = c(-1,0,2,4,6,8),file="regulation/data/TFs/xtables/LRT_activity_droppingTFs.tex")




# GI ####
#individual effects
GI_indSel_owms<-select_indiv_motifs("scores_top22TFs_mat",pheno="log2(GI)", intron_motifs_1_OWMs$SYMBOL_clean, "tree.intron.tfs22", cutoff=0.05) #CTCF
GI_all_lrts<-GI_indSel_owms[[3]]

TFs2<-bind_rows(GI_all_lrts, .id="TF") %>% filter(`p-value`<0.05)
GI_sign_lrts_df<-bind_rows(GI_all_lrts, .id="TF") %>% filter(TF %in% TFs2$TF)


mod_GI_CTCF<-gls(log2(GI)~CTCF, data=scores_top22TFs_mat, correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
mod_GI_null<-gls(log2(GI)~1, data=scores_top22TFs_mat, correlation=corBrownian(1, tree.intron.tfs22, form=~species), method="ML")
R2.resid(mod_GI_CTCF,mod_GI_null) #0.41

#diagnostics
qqnorm(mod_GI_CTCF$residuals, pch = 1, frame = FALSE)
qqline(mod_GI_CTCF$residuals, col = "steelblue", lwd = 2)
plot(mod_GI_CTCF$fitted, mod_GI_CTCF$residuals)


#save the fitted values
scores_top22TFs_mat$fitted_GI_CTCF<-mod_GI_CTCF$fitted
saveRDS(scores_top22TFs_mat, "regulation/data/TFs/ClusterBuster/results/scores_top22TFs_mat.rds")

drop_indiv_GI<-as.data.frame(summary(mod_GI_CTCF)$tTable) %>% rownames_to_column("motif")

#save the coefficients from both models for plotting
combined_coeffs<-bind_rows(drop_indiv_GI %>% mutate(y="GI"),
                           drop_combined %>% mutate(y="intron activity")) %>%
  mutate(y=factor(y, levels=c("intron activity","GI")))
saveRDS(combined_coeffs,"regulation/data/TFs/ClusterBuster/results/PGLS_3TF_coeffs.rds")





################################################
# VISUALIZE BINDING SCORES - prepare tables ####
################################################

scores_top22TFs_long2<-scores_top22TFs_long %>%
  group_by(SYMBOL_clean) %>%
  mutate(mean=mean(score), sd=sd(score), max=max(score)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(stand_score=(score-mean)/sd,
         rel_score=score/max,
         centered_score=score-mean) 

#order the data in the same order as the tree tips
tf_tree<-ggtree(tree.intron.tfs22)
tf_tree_ord<-tf_tree$data[which(tf_tree$data$isTip),]
desiredOrder_tf  <- rev(tf_tree_ord[order(tf_tree_ord$y),]$label)

scores_top22TFs_long2$order<-match(scores_top22TFs_long2$species, desiredOrder_tf)
scores_top22TFs_long2$species2<-gsub("_"," ",scores_top22TFs_long2$species)
saveRDS(scores_top22TFs_long2,"regulation/data/TFs/ClusterBuster/results/scores_22TFs_long.rds")

tree.intron.tfs22.2<-tree.intron.tfs22
tree.intron.tfs22.2$tip.label<-gsub("_"," ",tree.intron.tfs22.2$tip.label)
saveRDS(tree.intron.tfs22.2,"regulation/data/TFs/ClusterBuster/results/tree_22TFs.rds")


