
library(tidyverse)
library(xtable)


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
  
  return(full_cor_df)
}
