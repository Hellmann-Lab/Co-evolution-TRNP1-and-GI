#function

summarize_cor_output<-function(path_to_cov_file){
  L <- readLines(paste0(path_to_cov_file))

  #get the included parameters
  params1<-which(grepl("entries are", L)==T)
  params2<-which(grepl("covariances", L)==T)
  params<-L[(params1+1):(params2-2)]

  #correlations & their posterior probabilities
  cor1<-which(grepl("^correlation coefficients$", L)==T)
  cor2<-which(grepl("posterior", L)==T)[1]
  cors<-data.frame(L[(cor1+2):(cor2-2)]) %>% 
    tidyr::separate(col = 1, into=paste0(params), sep="\t")
  rownames(cors)<-params
  cors_long<-cors %>% 
    tibble::rownames_to_column("param1") %>% 
    tidyr::gather("param2","marginal_cor",2:ncol(.)) %>%
    dplyr::mutate(marginal_cor=round(as.numeric(marginal_cor),3))


  pp1<-which(grepl("^posterior", L)==T)[1]
  pp2<-which(grepl("precisions", L)==T)
  pps<-data.frame(L[(pp1+2):(pp2-2)]) %>%
    tidyr::separate(col = 1, into=paste0(params), sep="\t") #it's the last \t that is being excluded
  rownames(pps)<-params
  pps_long<-pps %>% 
    tibble::rownames_to_column("param1") %>% 
    tidyr::gather("param2","marginal_pp",2:ncol(.)) %>%
    dplyr::mutate(marginal_pp=round(as.numeric(marginal_pp),2))


  #now, partial cors and their pps
  partial_cor1<-which(grepl("^partial correlation coefficients$", L)==T)
  partial_cor2<-which(grepl("posterior", L)==T)[2]
  partial_cors<-data.frame(L[(partial_cor1+2):(partial_cor2-2)]) %>% 
    tidyr::separate(col = 1, into=paste0(params), sep="\t")
  rownames(partial_cors)<-params
  partial_cors_long<-partial_cors %>% 
    tibble::rownames_to_column("param1") %>% 
    tidyr::gather("param2","partial_cor",2:ncol(.)) %>%
    dplyr::mutate(partial_cor=round(as.numeric(partial_cor),3))


  partial_pp1<-which(grepl("^posterior", L)==T)[2]
  partial_pps<-data.frame(L[(partial_pp1+2):(length(L)-1)]) %>% 
    tidyr::separate(col = 1, into=paste0(params), sep="\t")
  rownames(partial_pps)<-params
  partial_pps_long<-partial_pps %>% 
    tibble::rownames_to_column("param1") %>% 
    tidyr::gather("param2","partial_pp",2:ncol(.)) %>%
    dplyr::mutate(partial_pp=round(as.numeric(partial_pp),2))



  #get them all into the same table

  full_cor_df<-cors_long %>%
    dplyr::left_join(pps_long) %>%
    dplyr::left_join(partial_cors_long) %>%
    dplyr::left_join(partial_pps_long) %>%
    dplyr::filter(!(is.na(marginal_pp)& grepl(1, marginal_cor))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(combi=paste(sort(c(param1,param2)), collapse=",")) %>%
    dplyr::distinct(combi, .keep_all=T) %>%
    #for the negative cors, pp is given by 1-pp. Turn it around again
    dplyr::mutate(marginal_pp=dplyr::case_when(marginal_cor<0 ~ (1-marginal_pp),
                                              T ~ marginal_pp)) %>%
    dplyr::mutate(partial_pp=dplyr::case_when(partial_cor<0 ~ (1-partial_pp),
                                              T ~ partial_pp))
}



#use the function (warning messages saying:
# 1: Expected 4 pieces. Additional pieces discarded in 4 rows [1, 2, 3, 4]. 
# 2: Expected 4 pieces. Additional pieces discarded in 4 rows [1, 2, 3, 4]. 
# 3: NAs introduced by coercion 
# 4: Expected 4 pieces. Additional pieces discarded in 4 rows [1, 2, 3, 4]. 
# 5: Expected 4 pieces. Additional pieces discarded in 4 rows [1, 2, 3, 4]. 
# 6: NAs introduced by coercion 
#are okay, as long as the pieces are the expected number (i.e., here 4). These warnings are due to the additional "\t" at the very end of each row that makes the function "separate" think that there should be one more column. and the NAs are because of the rounding up of pps --> it turns non-numerics (-) into NAs. no difference for us.

setwd("/data/share/htp/TRNP1/paper_data/")

full_cor_df_31_3<-summarize_cor_output("protein/coevol/results/all/GI_EQ_body_mass_31species/exon1_GI_EQ_body_mass_dsomggc1.cov")

full_cor_df_31_3$param2<-factor(full_cor_df_31_3$param2, levels=c("omega","dS","EQ","body_mass"))

full_cor_df_31_3<-full_cor_df_31_3 %>%
  arrange(param2, desc(marginal_cor)) %>%
  mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
         partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
  dplyr::select(-combi, -marginal_pp, -partial_pp)


#prettify
full_cor_df_pretty<-full_cor_df_31_3
full_cor_df_pretty$param1<-gsub("_"," ",full_cor_df_pretty$param1)
full_cor_df_pretty$param2<-gsub("_"," ",full_cor_df_pretty$param2)
names(full_cor_df_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")

print(xtable(full_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_GI_EQ_body_mass_31sp.txt")





# THE FINAL ONE: brain, body mass, GI together
full_cor_df_31_3<-summarize_cor_output("protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.cov")

full_cor_df_31_3$param2<-factor(full_cor_df_31_3$param2, levels=c("omega","dS","brain_mass","body_mass"))

full_cor_df_31_3<-full_cor_df_31_3 %>%
  arrange(param2, desc(marginal_cor)) %>%
  mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
         partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
  dplyr::select(-combi, -marginal_pp, -partial_pp)


#prettify
full_cor_df_pretty<-full_cor_df_31_3
full_cor_df_pretty$param1<-gsub("_"," ",full_cor_df_pretty$param1)
full_cor_df_pretty$param2<-gsub("_"," ",full_cor_df_pretty$param2)
names(full_cor_df_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")

print(xtable(full_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_GI_brain_mass_body_mass_31sp.txt")




#the individual cors
#GI alone
full_cor_df_31_GI<-summarize_cor_output("protein/coevol/results/all/GI_31species/exon1_GI_dsomggc1.cov")
full_cor_df_31_GI$param2<-factor(full_cor_df_31_GI$param2, levels=c("omega","dS"))

full_cor_df_31_GI<-full_cor_df_31_GI %>%
  arrange(param2, desc(marginal_cor)) %>%
  mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
         partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
  dplyr::select(-combi, -marginal_pp, -partial_pp)

#prettify
GI_cor_df_pretty<-full_cor_df_31_GI
GI_cor_df_pretty$param1<-gsub("_"," ",GI_cor_df_pretty$param1)
GI_cor_df_pretty$param2<-gsub("_"," ",GI_cor_df_pretty$param2)
names(GI_cor_df_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")

print(xtable(GI_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_GI_31sp.txt")




#brain mass alone
full_cor_df_31_brain<-summarize_cor_output("protein/coevol/results/all/brain_mass_31species/exon1_brain_mass_dsomggc1.cov")
full_cor_df_31_brain$param2<-factor(full_cor_df_31_brain$param2, levels=c("omega","dS"))

full_cor_df_31_brain<-full_cor_df_31_brain %>%
  arrange(param2, desc(marginal_cor)) %>%
  mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
         partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
  dplyr::select(-combi, -marginal_pp, -partial_pp)

#prettify
brain_cor_df_pretty<-full_cor_df_31_brain
brain_cor_df_pretty$param1<-gsub("_"," ",brain_cor_df_pretty$param1)
brain_cor_df_pretty$param2<-gsub("_"," ",brain_cor_df_pretty$param2)
names(brain_cor_df_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")

print(xtable(brain_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_brain_mass_31sp.txt")




#brain mass alone
full_cor_df_31_body<-summarize_cor_output("protein/coevol/results/all/body_mass_31species/exon1_body_mass_dsomggc1.cov")
full_cor_df_31_body$param2<-factor(full_cor_df_31_body$param2, levels=c("omega","dS"))

full_cor_df_31_body<-full_cor_df_31_body %>%
  arrange(param2, desc(marginal_cor)) %>%
  mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
         partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
  dplyr::select(-combi, -marginal_pp, -partial_pp)

#prettify
body_cor_df_pretty<-full_cor_df_31_body
body_cor_df_pretty$param1<-gsub("_"," ",body_cor_df_pretty$param1)
body_cor_df_pretty$param2<-gsub("_"," ",body_cor_df_pretty$param2)
names(body_cor_df_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")

print(xtable(body_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_body_mass_31sp.txt")







#compare the different models ####
full_cor_df_31_1_EQ<-summarize_cor_output("protein/coevol/results/all/GI_EQ_body_mass_31species/exon1_GI_EQ_body_mass_dsomggc1.cov") %>%
  mutate(model="GI_EQ_body_mass")

#now do GI, body mass, brain mass
full_cor_df_31_1_BM<-summarize_cor_output("protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.cov") %>%
  mutate(model="GI_brain_mass_body_mass")


#now do all 4 separately
GI_cor_df_31_1<-summarize_cor_output("protein/coevol/results/all/GI_31species/exon1_GI_dsomggc1.cov") %>%
  mutate(model="GI_only")

EQ_cor_df_31_1<-summarize_cor_output("protein/coevol/results/all/EQ_31species/exon1_EQ_dsomggc1.cov") %>%
  mutate(model="EQ_only")

body_cor_df_31_1<-summarize_cor_output("protein/coevol/results/all/body_mass_31species/exon1_body_mass_dsomggc1.cov") %>%
  mutate(model="body_only")

brain_cor_df_31_1<-summarize_cor_output("protein/coevol/results/all/brain_mass_31species/exon1_brain_mass_dsomggc1.cov") %>%
  mutate(model="brain_only")


all_coevol_models<-bind_rows(full_cor_df_31_1_EQ,
                             full_cor_df_31_1_BM,
                             GI_cor_df_31_1, EQ_cor_df_31_1, body_cor_df_31_1, brain_cor_df_31_1) %>%
  #filter(grepl("omega",combi)) %>%
  filter(!grepl("dS",combi)) %>%
  mutate(combi=case_when(model=="GI_EQ_body_mass" ~ "GI_EQ_body_mass",
                         model=="GI_brain_mass_body_mass" ~ "GI_brain_mass_body_mass",
                         T ~ "single"))


ggplot(all_coevol_models %>% filter(param2=="omega"), aes(x=param1, y=marginal_cor, color=combi, size=marginal_pp))+ geom_point() + scale_color_manual(values=c("#646E68","#9E768F","grey"))+theme(panel.background = element_rect(fill="grey95"), axis.title.x = element_blank())+ylab("Marginal correlation with omega")+
  scale_size_continuous(range=c(3,7))

ggsave("protein/coevol/results/all/figures/compare_models/marginal_cor_models.pdf", height=4, width=7)




all_coevol_models_long<-all_coevol_models %>%
  filter(param2=="omega") %>%
  #mutate(param_combi=paste0(param1,"_",param2)) %>%
  reshape(direction="long",
          varying = c("marginal_cor","marginal_pp","partial_cor","partial_pp"),
          #timevar="model",
          times=c("marginal","partial"),
          v.names=c("pp","cor")) %>%
  mutate(line_grouping=paste0(param1,"_",model))


ggplot(all_coevol_models_long, aes(x=param1, y=cor, color=combi, size=pp,shape=time))+ 
  geom_point(aes(group=time), position=position_dodge(width=0.40)) + 
  #geom_line(aes(group=time), position=position_dodge(width=0.4)) +
  scale_color_manual(values=c("#646E68","#9E768F","grey"))+
  theme(panel.background = element_rect(fill="grey95"), axis.title.x = element_blank())+
  ylab("Correlation with omega")+
  scale_size_continuous(range=c(3,7))


ggplot(all_coevol_models_long, aes(x=time, y=cor, color=combi, group=model))+ 
  geom_point(aes(size=pp,shape=time)) + 
  geom_line() +
  scale_color_manual(values=c("#646E68","#9E768F","grey"))+
  theme(panel.background = element_rect(fill="grey95"), axis.title.x = element_blank())+
  facet_grid(.~param1)+ylab("Correlation with omega")+
  scale_size_continuous(range=c(2,7))

ggsave("protein/coevol/results/all/figures/compare_models/partial_cor_models.pdf", height=4.5, width=8)





#no dolphin

full_cor_df_30_1_EQ<-summarize_cor_output("protein/coevol/results/no_dolphin/GI_EQ_body_mass_30species/exon1_GI_EQ_body_mass_dsomggc1.cov") %>%
  mutate(model="GI_EQ_body_mass")

