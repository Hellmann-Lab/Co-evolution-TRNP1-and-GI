# summarize coevol results (separate runs for each pheno)


#check convergence (rel_diff<0.1, effsize>300)####
setwd("/data/share/htp/TRNP1/paper_data/protein/coevol/results/all/tracecomp_summary/")
list.tracecomp<-list.files(recursive = TRUE)

for (file in list.tracecomp){
  f<-read.table(paste(file), header = T)
  print(table(f$effsize>300 | f$rel_diff<0.1))
}

#all converged




#comparing ln marginal likelihoods of the models####
#converged, so I'll randomly take 1 per model (f.e.,run3)

coevol_groups<-c("GI","GI_EQ_body_mass","brain_mass","GI_body_mass","GI_brain_mass","GI_EQ","body_mass","EQ","EQ_body_mass")

cor_out<-"/data/share/htp/TRNP1/paper_data/protein/coevol/results/all/correlation_output/"

ln_probs_coevol_31sp<-data.frame(ln_prob_mean=as.numeric(),
                                 ln_prob_sd=as.numeric(),
                                 combi=as.character())

for (i in coevol_groups){

  ln_prob_run1<-read.table(paste0(cor_out,i,"_31sp/ln_prob_",i,"_run1.txt"))$V8
  ln_prob_run2<-read.table(paste0(cor_out,i,"_31sp/ln_prob_",i,"_run2.txt"))$V8
  ln_prob_run3<-read.table(paste0(cor_out,i,"_31sp/ln_prob_",i,"_run3.txt"))$V8
  
  ln_prob_mean_runs<-data.frame(ln_prob_mean=mean(c(ln_prob_run1,ln_prob_run2,ln_prob_run3)),
                                ln_prob_sd=sd(c(ln_prob_run1,ln_prob_run2,ln_prob_run3)),
                                combi=i)
  ln_probs_coevol_31sp<-bind_rows(ln_probs_coevol_31sp,ln_prob_mean_runs)
}

#null model (run3)
null_ln_prob_run1<-(-4789.88)
null_ln_prob_run2<-(-4736.85)
#null_ln_prob_run3<-(-4765.26) --> exclude this one cause it was prematurely stopped and has not fully converged with the others

ln_probs_coevol_31sp<-ln_probs_coevol_31sp %>%
  add_row(ln_prob_mean=mean(c(null_ln_prob_run1,null_ln_prob_run2)),
          ln_prob_sd=sd(c(null_ln_prob_run1,null_ln_prob_run2)),
          combi="null")

# computing bayes factor:
# 1) based on marginal likelihoods of the model meant for estimating the joint posterior probability of the model parameters
# 2) does not matter which one is first (only for interpretation); these models do not need to be nested
# 3) BF=exp(ln(D|M0)-ln(D|M1)) --> above 10 already suggests strong support for the first model (in log space: above 2.3) 
# source https://revbayes.github.io/tutorials/model_selection_bayes_factors/bf_intro.html



# make a function
calculate_BF<-function(ln_probs_coevol_31sp,model1, model2){
  df<-data.frame(model1=paste0(model1),
                 model2=paste0(model2),
                 ln_prob_m1=ln_probs_coevol_31sp$ln_prob_mean[ln_probs_coevol_31sp$combi==paste(model1)],
                 ln_prob_m2=ln_probs_coevol_31sp$ln_prob_mean[ln_probs_coevol_31sp$combi==paste(model2)])
  df<-df %>%
    mutate(BF=exp(ln_prob_m1 - ln_prob_m2),
           model1_support=case_when(BF<1 ~ "negative",
                     BF>1 & BF<3.2 ~ "weak",
                     BF>3.2 & BF<10 ~ "substantial",
                     BF>10 & BF<100 ~ "strong",
                     BF>100 ~ "decisive"))
  return(df)
}




BF_summary<-bind_rows(calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ_body_mass",model2="GI_EQ"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ_body_mass",model2="GI_body_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ_body_mass",model2="EQ_body_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ",model2="GI"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ",model2="EQ"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_body_mass",model2="GI"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_body_mass",model2="body_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="EQ_body_mass",model2="EQ"),
                      calculate_BF(ln_probs_coevol_31sp,model1="EQ_body_mass",model2="body_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_brain_mass",model2="GI"),
                       calculate_BF(ln_probs_coevol_31sp,model1="GI_brain_mass",model2="brain_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_body_mass",model2="GI_EQ"),
                       calculate_BF(ln_probs_coevol_31sp,model1="GI_body_mass",model2="GI_brain_mass"),
                      calculate_BF(ln_probs_coevol_31sp,model1="GI_EQ",model2="GI_brain_mass"),
                       calculate_BF(ln_probs_coevol_31sp,model1="brain_mass",model2="null"),
                      #need to insert the next ones after running the null model
                      calculate_BF(ln_probs_coevol_31sp,model1="GI",model2="null"),
                      calculate_BF(ln_probs_coevol_31sp,model1="EQ",model2="null"),
                      calculate_BF(ln_probs_coevol_31sp,model1="body_mass",model2="null"))











#summarize correlation across runs ####
setwd("/data/share/htp/TRNP1/paper_data/protein/coevol/results/all/correlation_output")

omega_pheno_cor_func<-function(name_phenos,phenos){
  
  omega_pheno<-list()
  for (run in 1:3){
  
  cors<-read.table(paste0("./",name_phenos,"_30sp/cors_",name_phenos,"_run",run,".txt"), skip = 1)
  rownames(cors)<-read.table(paste0("./",name_phenos,"_30sp/names_",name_phenos,"_run",run,".txt"))$V1
  colnames(cors)<-read.table(paste0("./", name_phenos, "_30sp/names_",name_phenos,"_run",run,".txt"))$V1
  cor_df<-cors %>% rownames_to_column("param1") %>% 
    gather("param2","rho",2:ncol(.)) %>%
    #check only the omega vs pheno for now
    filter(param1=="omega" & param2 %in% phenos) 
  
  pps<-read.table(paste0("./",name_phenos,"_30sp/pp_",name_phenos,"_run",run,".txt"), skip = 1)
  rownames(pps)<-read.table(paste0("./",name_phenos,"_30sp/names_",name_phenos,"_run",run,".txt"))$V1
  colnames(pps)<-read.table(paste0("./",name_phenos,"_30sp/names_",name_phenos,"_run",run,".txt"))$V1
  pp_df<-pps %>% rownames_to_column("param1") %>% 
    gather("param2","pp",2:ncol(.)) %>%
    filter(param1=="omega" & param2 %in% phenos) 
  
  omega_pheno[[run]]<-inner_join(cor_df, pp_df) 
}
  omega_pheno_df<-bind_rows(omega_pheno, .id="run")
  return(omega_pheno_df)
}

omega_pheno_cor_func(name_phenos = "GI_body_mass",phenos = c("GI","body_mass") )

GI_runs<-omega_pheno_cor_func(name_phenos = "GI",phenos = c("GI") )
brain_mass_runs<-omega_pheno_cor_func(name_phenos = "brain_mass",phenos = c("brain_mass") )
EQ_runs<-omega_pheno_cor_func(name_phenos = "EQ",phenos = c("EQ") )
body_mass<-omega_pheno_cor_func(name_phenos = "body_mass",phenos = c("body_mass") )

sep_cors<-bind_rows(GI_runs,brain_mass_runs, EQ_runs, body_mass)
ggplot(sep_cors, aes(x=reorder(param2, rho), y=rho, color=as.numeric(pp)))+
  geom_jitter(width=0.2, size=2.5)+
  theme_bw()+
  theme(axis.title.x = element_blank())
ggsave("cor_summary.pdf", height=3, width=5)
