#combining dN/dS estimates of the protein with the intronic activity and running PGLS on this
#load libraries 
libs<-c("tidyverse","cowplot","GenomicRanges","data.table", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","grid","gtable","xtable","rr2")
sapply(libs, require, character.only=T)


################
# load data ####
################

#get the species names with protein seq available
trnp1.coding.all.fa<-readDNAStringSet('protein/fastas/TRNP1_coding_seqs_45sp_full.fa')
trnp1.coding.all.df<-data.frame(species=names(trnp1.coding.all.fa)) %>%
  mutate(species2=gsub("_","",species),
         species2=substr(species2,1,10))

#intersect with coevol results to get the full species name
coevol_estimates_31sp<-readRDS("protein/coevol/results/for_figures/coevol_3phenos_31sp_summarized.rds") %>%
  dplyr::rename(species2=species) %>%
  right_join(trnp1.coding.all.df) %>%
  dplyr::select(species,GI,omega) %>%
  filter(!is.na(GI))

#add intron activity
coevol_estimates_9sp_PGLS<-coevol_estimates_31sp %>%
  inner_join(readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds") %>%  
  filter(region=="intron", primate_clade %in% c("Old World monkey", "Great ape"), cell_line=="human1")) 

#adjust tree
mammaltree<-read.tree("protein/trees/mammaltree.txt")
tree.intron.dnds.GI<-drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% coevol_estimates_9sp_PGLS$species])
ggtree(tree.intron.dnds.GI)+geom_tiplab()+xlim(0,200)



##########################
# standardized values ####
##########################

#standardize activity and dNdS (subtract the mean and divide by the standard deviation)
values_to_standardize<-coevol_estimates_9sp_PGLS %>%
  dplyr::select(species,omega,log2_total_activity) %>%
  mutate(log2_omega=log2(omega)) %>%
  gather(group,value, "log2_total_activity", "log2_omega") %>%
  group_by(group) %>%
  mutate(mean=mean(value), sd=sd(value)) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(stand=(value-mean)/sd,
         group_stand=paste0(group,"_stand")) %>%
  dplyr::select(species, group_stand, stand) %>%
  spread(group_stand,stand)

coevol_estimates_9sp_PGLS_stand<-left_join(coevol_estimates_9sp_PGLS,values_to_standardize)
rownames(coevol_estimates_9sp_PGLS_stand)<-coevol_estimates_9sp_PGLS_stand$species


###########
# PGLS ####
###########

#null model
mod.null.stand<-gls(log2(GI)~1, data=coevol_estimates_9sp_PGLS_stand, 
                    correlation=corBrownian(value=1,phy=tree.intron.dnds.GI,form=~species), 
                    method="ML")

#individual - only omega as a predictor
mod.dnds.stand<-gls(log2(GI)~log2_omega_stand, data=coevol_estimates_9sp_PGLS_stand, 
                    correlation=corBrownian(value=1,phy=tree.intron.dnds.GI,form=~species), 
                    method="ML")
anova(mod.dnds.stand,mod.null.stand)
mod_GI_dnds_stand<-wrap_summary_table_BM(mod.dnds.stand) %>%
  mutate(R2=round(R2.resid(mod.dnds.stand,mod.null.stand),digits=3))


#combined - omega and intron activity as predictors
mod.dnds.intron.stand<-gls(log2(GI)~log2_omega_stand+log2_total_activity_stand, 
                           data=coevol_estimates_9sp_PGLS_stand, 
                           correlation=corBrownian(value=1,phy=tree.intron.dnds.GI,form=~species),
                           method="ML")
summary(mod.dnds.intron.stand)
anova(mod.dnds.intron.stand,mod.dnds.stand)
mod_GI_dnds_intron_stand<-wrap_summary_table_BM(mod.dnds.intron.stand)%>%
  mutate(R2=round(R2.resid(mod.dnds.intron.stand,mod.null.stand),digits=3))

#LRT table
pgls_combined_mod_LRT_stand<-anova(mod.dnds.stand, mod.dnds.intron.stand) %>%
  mutate(call=gsub("(gls[(]model = )","",call),
         call=gsub("\\,.*","",call),
         call=gsub("log2_total_activity_stand","log2(stand.intron)",call),
         call=gsub("log2_omega_stand","log2(stand.omega)",call),
         df=as.factor(df)) 

#combine
comb_stand<-bind_rows(mod_GI_dnds_stand,mod_GI_dnds_intron_stand) %>%
  filter(Predictor!="(Intercept)") %>%
  mutate(Predictor=gsub("log2_total_activity_stand","log2(stand.intron)",Predictor),
         Predictor=gsub("log2_omega_stand","log2(stand.omega)",Predictor),
         mod=c(1,2,2)) %>%
  left_join(pgls_combined_mod_LRT_stand[,c("call","df","L.Ratio","p-value")] %>% 
              mutate(mod=c(1,2))) %>%
  dplyr::select(call,Predictor,Value,Std.Error,df,logLik,L.Ratio,`p-value`) %>%
  dplyr::rename(Model=call, `LRT p-value`=`p-value`)
comb_stand[2,c("df","logLik","L.Ratio","LRT p-value")]<-NA

print(xtable(comb_stand,digits=c(1,1,1,2,3,1,2,2,3)),include.rownames=FALSE, file="regulation/data/MPRA/xtables/LRT_coeff_dnds_intron_stand_OWMs_apes.tex")


