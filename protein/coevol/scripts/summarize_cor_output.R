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

prettify_table<-function(table){
  table<-table %>%
    dplyr::arrange(param2, desc(marginal_cor)) %>%
    mutate(marginal_cor=paste0(marginal_cor," (", marginal_pp, ")"),
           partial_cor=paste0(partial_cor," (", partial_pp, ")")) %>%
    dplyr::select(-combi, -marginal_pp, -partial_pp)
  
  table_pretty<-table
  table_pretty$param1<-gsub("_"," ",table_pretty$param1)
  table_pretty$param2<-gsub("_"," ",table_pretty$param2)
  names(table_pretty)<-c('Parameter 1', "Parameter 2", "Marginal Correlation (Posterior Probability)", "Partial Correlation (Posterior Probability)")
  return(table_pretty)
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


# THE FINAL ONE: brain, body mass, GI together ####
full_cor_df_31_3<-summarize_cor_output("protein/coevol/results/all/body_mass_brain_mass_GI_31species/exon1_body_mass_brain_mass_GI_dsomggc1.cov")

full_cor_df_31_3$param2<-factor(full_cor_df_31_3$param2, levels=c("omega","dS","brain_mass","body_mass"))
#prettify
full_cor_df_pretty<-prettify_table(full_cor_df_31_3)
print(xtable(full_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_GI_brain_mass_body_mass_31sp.txt")




#the individual cors
#GI alone ####
full_cor_df_31_GI<-summarize_cor_output("protein/coevol/results/all/GI_31species/exon1_GI_dsomggc1.cov")
full_cor_df_31_GI$param2<-factor(full_cor_df_31_GI$param2, levels=c("omega","dS"))
#prettify
GI_cor_df_pretty<-prettify_table(full_cor_df_31_GI)
print(xtable(GI_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_GI_31sp.txt")




#brain mass alone ####
full_cor_df_31_brain<-summarize_cor_output("protein/coevol/results/all/brain_mass_31species/exon1_brain_mass_dsomggc1.cov")
full_cor_df_31_brain$param2<-factor(full_cor_df_31_brain$param2, levels=c("omega","dS"))

#prettify
brain_cor_df_pretty<-prettify_table(full_cor_df_31_brain)
print(xtable(brain_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_brain_mass_31sp.txt")




#body mass alone ####
full_cor_df_31_body<-summarize_cor_output("protein/coevol/results/all/body_mass_31species/exon1_body_mass_dsomggc1.cov")
full_cor_df_31_body$param2<-factor(full_cor_df_31_body$param2, levels=c("omega","dS"))

#prettify
body_cor_df_pretty<-prettify_table(full_cor_df_31_body)
print(xtable(body_cor_df_pretty), include.rownames=FALSE, file="protein/coevol/output_xtables/full_cor_body_mass_31sp.txt")

sep_cors_together<-bind_rows(GI_cor_df_pretty, brain_cor_df_pretty, body_cor_df_pretty)
print(xtable(sep_cors_together), hline.after = c(-1,0,3,6,9),include.rownames=FALSE, file="protein/coevol/output_xtables/sep_cors_combined_31sp.txt")

