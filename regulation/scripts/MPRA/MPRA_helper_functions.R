#helper functions for MPRA analysis

#LOAD ALL THE PACKAGES

`%ni%`<-Negate('%in%')

#wrapping up PGLS output summary (for predicting EQ or GI using region activities)
wrap_summary_table<-function(model_output){
  repl_summary<-summary(model_output)
  repl_summary_table<-as.data.frame(repl_summary$tTable) %>%
    rownames_to_column(var="region")%>%
    dplyr::rename(p.value="p-value", t.value="t-value") %>%
    mutate(Value=round(Value,3),
           p.value=as.numeric(formatC(p.value, format = "e", digits = 1)),
           modelStruct_param=round(as.numeric(model_output$modelStruct), digits=2),
           logLik=model_output$logLik)
  return(repl_summary_table)
}




PGLS_pheno_vs_MPRAactivity<-function(activity_overlap_summary, tree, cell_lines, regions, phenos, pheno_data_table,corStructs){
  
  model_output<-list()
  model_sel_output<-list()
  
  for (i in regions){
    for (clines in cell_lines){
      for (pheno in phenos){
        
        print(i)
        print(clines)
        print(pheno)
        
        activity_overlap_summary_wide<-activity_overlap_summary %>%
          dplyr::filter(region==paste0(i) & cell_line==clines) %>%
          left_join(pheno_data_table[,c(paste0(pheno), "species")]) %>%
          drop_na()
        
        mpra_tree<-ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% activity_overlap_summary_wide$species])
        
        rownames(activity_overlap_summary_wide)<-activity_overlap_summary_wide$species
        name.check(mpra_tree,activity_overlap_summary_wide)
       
        form1<-paste0("log2(",pheno,")~log2_total_activity")
        form_null<-paste0("log2(",pheno,")~1")
      
        # BROWNIAN MOTION CORRELATION STRUCTURE #

          #full model
          mod.totalA.BM.pre<-paste("gls(",form1,", data=activity_overlap_summary_wide, correlation=corBrownian(1,mpra_tree,form=~species),method='ML')")
          mod.totalA.BM<-eval(parse(text = mod.totalA.BM.pre))  
          mod_BM<-wrap_summary_table(mod.totalA.BM)
          mod_BM$region<-gsub("log2_total_activity", paste0("log2(",i,")"), mod_BM$region)
          mod_BM$model<-paste0("log2(",pheno,")~log2(",i,")")
          mod_BM$cor_structure<-paste0("BM")
          mod_BM$df<-Dim(summary(mod.totalA.BM)$modelStruct$corStruct)$N
          
          
          #null model
          mod.totalA.BM.null.pre<-paste("gls(",form_null,", data=activity_overlap_summary_wide, correlation=corBrownian(1, mpra_tree,form=~species), method='ML')")
          mod.totalA.BM.null<-eval(parse(text = mod.totalA.BM.null.pre))

          #compare full and null model using anova & MLs of the models
          aov_res_BM<-anova(mod.totalA.BM.null,mod.totalA.BM)
          mod_BM$mod.sel.lrt<-aov_res_BM$`p-value`[2]
          mod_BM$cell_line<-paste0(clines)
          mod_BM$pheno<-paste0(pheno)
          #mod_BM$R2<-rr2::R2.resid(mod.totalA.BM,mod.totalA.BM.null)
          model_output[[paste0(i,clines,pheno)]]<-mod_BM
          #model_output<-bind_rows(model_output, mod_BM)
          
          model_sel_output[[paste0(i,clines,pheno)]]<-as.data.frame(aov_res_BM) %>%
            rownames_to_column("CorStr") %>%
            mutate(CorStr=gsub("mod.totalA.","",CorStr),
                   region=i,
                   cell_line=clines,
                   pheno=pheno)
      }
    }
  }
  return(list(model_output, model_sel_output))
}




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



#zoom into the intron tree - subclades ####
PGLSzoomin<-function(df=intron_activity_GI,tree=intron_tree_full,clade_type=NULL, clade_subset=NULL, cells="human1", subset_name="all"){
  
  df<-df %>% filter(.data[[clade_type]] %in% clade_subset,
                    cell_line %in% cells)
  subset_tree<-ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% df$species])
  rownames(df)<-df$species
  print(name.check(subset_tree,df))
  
  mod.GI.intron<-gls(log2(GI)~log2_total_activity, data=df, correlation=corBrownian(1,subset_tree,form=~species), method="ML")
  mod.GI.null<-gls(log2(GI)~1, data=df, correlation=corBrownian(1,subset_tree,form=~species), method="ML")
  R2_value<-round(rr2::R2.resid(mod.GI.intron,mod.GI.null),digits=3)
  
  lrt<-anova(mod.GI.null,mod.GI.intron)
  return(list(as.data.frame(lrt %>% mutate(call=lrt$call, Species=paste0(subset_name), `Cell line`=paste0(cells))),
              wrap_summary_table_BM(mod.GI.intron) %>% mutate(call=lrt$call, Species=paste0(subset_name), `Cell line`=paste0(cells)),mod.GI.null,R2_value))
  
}


PGLSzoomin_coeff<-function(df=intron_activity_GI,tree=intron_tree_full,clade_type=NULL, clade_subset=NULL, cells="human1", subset_name="all"){
  
  df<-df %>% filter(.data[[clade_type]] %in% clade_subset,
                    cell_line %in% cells)
  subset_tree<-ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% df$species])
  rownames(df)<-df$species
  print(name.check(subset_tree,df))
  
  mod.GI.intron<-gls(log2(GI)~log2_total_activity, data=df, correlation=corBrownian(1,subset_tree,form=~species), method="ML")
  mod.GI.null<-gls(log2(GI)~1, data=df, correlation=corBrownian(1,subset_tree,form=~species), method="ML")
  R2_value<-round(rr2::R2.resid(mod.GI.intron,mod.GI.null),digits=3)
  
  lrt<-anova(mod.GI.null,mod.GI.intron)
  
  return(list(as.data.frame(lrt %>% mutate(call=lrt$call, Species=paste0(subset_name), `Cell line`=paste0(cells))),
              mod.GI.intron,mod.GI.null,R2_value))
  
}





# topGO on the TFs ####

do_topGO_motifs_ENSEMBL<-function(all_genes_table, cutoff, nodesize=10){
  rownames(all_genes_table)<-all_genes_table$ENSEMBL
  sign<-all_genes_table$score>=cutoff
  List<-rep(1,dim(all_genes_table)[1])
  names(List)<-rownames(all_genes_table)
  List[sign]<-0.01
  topDiffGenes<-function(x){return(x==0.01)}
  newClass<-new("topGOdata",description="sign",
                ontology="BP", allGenes=List,
                geneSel=topDiffGenes, annot=annFUN.org,
                nodeSize=nodesize,ID="ENSEMBL",
                mapping="org.Hs.eg.db")
  numSigGenes(newClass)
  GO.res<-runTest(newClass,algorithm = "elim",statistic = "fisher")
  sum(score(GO.res)<=0.01)
  topGO_table<-GenTable(newClass,Fisher=GO.res,orderBy="Fisher",ranksOf="Fisher",topNodes=30)
  return(list(newClass,topGO_table,GO.res))
}






# predict activity or GI using motif scores - dropping them one by one if anova does not prefer the full model over the reduced one

wrap_summary_table_motifs<-function(model_output){
  repl_summary<-summary(model_output)
  repl_summary_table<-as.data.frame(repl_summary$tTable) %>%
    rownames_to_column(var="motif") %>% 
    dplyr::rename(p.value="p-value", t.value="t-value") %>%
    mutate(Value=round(Value, 3),
           p.value=as.numeric(formatC(p.value, format = "e", digits = 1)),
           logLik=model_output$logLik)
  return(repl_summary_table)
}


select_model_by_dropping<-function(data_wide,pheno="log2(GI)", top_list, tree, cor="corBrownian", lambda=1, cutoff=0.05){
  
  #initiate collecting-lists & the first - full model
  l<-list()
  l_aov<-list()
  aov_res<-as.data.frame(matrix(1,2,9))
  
  top_list<-top_list %>% filter(motif!="(Intercept)")
  
  fun_mod_full <- paste0(pheno,"~", paste(top_list[,"motif"],  collapse="+"))
  fun_mod_full2<-paste("gls(",fun_mod_full,", data=",data_wide,", correlation=",cor,"(",lambda,",", tree,",form=~species), method='ML')")
  res_full<-eval(parse(text = fun_mod_full2))
  mod_full_df<-wrap_summary_table_motifs(res_full)
  
  
  while (aov_res[2,9]>cutoff){
    
    #idenitfy the least significant motif
    least_sign<-mod_full_df %>% 
      filter(motif!="(Intercept)") %>%
      top_n(n=1, wt=p.value)
    
    top_list_red<-mod_full_df %>% 
      filter(motif!="(Intercept)" & !motif %in% least_sign$motif)
    
    if(dim(top_list_red)[1]==0){
      top_list_red[1,]<-1
    }
    
    fun_mod_red <- paste0(pheno,"~", paste(top_list_red[,"motif"],  collapse="+"))
    fun_mod_red2<-paste("gls(",fun_mod_red,", data=",data_wide,", correlation=",cor,"(",lambda,",", tree,",form=~species), method='ML')")
    res_red<-eval(parse(text = fun_mod_red2))
    mod_red_df<-wrap_summary_table_motifs(res_red)
    
    #compare the 2 models
    aov_res<-anova(res_full, res_red)
    
    #if the p-value is above the cutoff, drop the predictor in question by setting the reduced model to be the full one
    if (aov_res[2,9]>cutoff){
      res_full<-res_red
      mod_full_df<-mod_red_df
      
    }
  } 
  return(mod_full_df)
}




select_indiv_motifs<-function(data_wide,pheno,tf_list,tree,cutoff){
  newlist_aov<-list()
  newlist_res<-list()
  fulllist_aov<-list()
  
  for (i in tf_list){
    
    fun_mod_full <- paste0(pheno,"~", paste(i))
    fun_mod_full2<-paste("gls(",fun_mod_full,", data=",data_wide,", correlation=corBrownian(1,", tree,",form=~species), method='ML')")
    res_full<-eval(parse(text = fun_mod_full2))
    
    fun_mod_null <- paste0(pheno,"~1")
    fun_mod_null2<-paste("gls(",fun_mod_null,", data=",data_wide,", correlation=corBrownian(1,", tree,",form=~species), method='ML')")
    res_null<-eval(parse(text = fun_mod_null2))
    res_full_df<-wrap_summary_table_motifs(res_full)
    
    aov_res<-anova(res_null,res_full)
    fulllist_aov[[i]]<-aov_res
    
    if (aov_res[2,9]<cutoff){
      print(paste(i, "significant"))
      print(aov_res)
      print(res_full_df)
      newlist_aov[[i]]<-aov_res
      newlist_res[[i]]<-res_full_df
      
    } else 
    {print(paste0("drop ",i))}
    
  }
  return(list(newlist_aov,newlist_res,fulllist_aov))
  }









#run ctrain and cbust####
run_ctrain_cbust<-function(tf_motif_vector,JASPAR_folder_name, query_fasta,add_weights=TRUE, gap=NULL, cluster_cutoff=5, motif_cutoff=6){
  
  expr_tf_list_NPCs<-list()
  for (i in tf_motif_vector){
    path<-paste0("regulation/data/TFs/JASPAR_2020/JASPAR_collection/",i,".jaspar")
    expr_tf_list_NPCs[[i]]<-path
  }
  system(paste0('rm -r regulation/data/TFs/JASPAR_2020/JASPAR_collection_transposed/',JASPAR_folder_name))
  system(paste0('mkdir regulation/data/TFs/JASPAR_2020/JASPAR_collection_transposed/',JASPAR_folder_name))
  
  expr_tf_list_NPCs<-unlist(expr_tf_list_NPCs, use.names = F)
  tf_list_NPCs<-list()
  
  for (i in 1:length(expr_tf_list_NPCs)){
    jaspar_example<-read.table(paste0(expr_tf_list_NPCs[i]), skip=1)
    jaspar_transposed<-t(jaspar_example[,-c(1,2,length(jaspar_example))])
    jaspar_name<-read.table(paste0(expr_tf_list_NPCs[i]), nrows=1, sep=";")
    jaspar_name$V1<-gsub("\t"," ", jaspar_name$V1)
    jaspar_name$V1<-gsub(">","", jaspar_name$V1)
    tf_list_NPCs[[i]]<-jaspar_name$V1
    write.table(jaspar_transposed, paste0("regulation/data/TFs/JASPAR_2020/JASPAR_collection_transposed/",JASPAR_folder_name,"/",jaspar_name), 
                row.names = F, col.names = F, quote = F, na="")
  }
  
  system(paste0('rm -r regulation/data/TFs/ClusterBuster/results/',JASPAR_folder_name))
  system(paste0('mkdir regulation/data/TFs/ClusterBuster/results/',JASPAR_folder_name))
  
  system(paste('sbatch scripts/TFs/run_clusterbuster.sh',JASPAR_folder_name,query_fasta,add_weights,cluster_cutoff,motif_cutoff,gap), wait = T)
  
}


#gathering cbust scores####
gather_scores<-function(JASPAR_folder_name){
  #gather scores into one file --all species
  # load in scores
  score.files <- list.files(path=paste0("regulation/data/TFs/ClusterBuster/results/",JASPAR_folder_name,"/subset_cbust/scores"), recursive = F, full.names=T, pattern=".txt")
  score.list <- lapply(score.files, function(x) read.table(x, skip=1))
  names(score.list) <- paste0(sapply(strsplit(score.files, "/"), "[[", 9))
  names(score.list)<-gsub(".txt","",names(score.list))
  # get the column names
  names_scores<-read.table(paste0("regulation/data/TFs/ClusterBuster/results/",JASPAR_folder_name,"/subset_cbust/scores/Homo_sapiens.txt"), skip = 1, nrows = 1, comment.char = "", stringsAsFactors = F)[,-1]
  scores_df<-plyr::ldply(score.list, .id="species")
  colnames(scores_df)[2:ncol(scores_df)]<-names_scores
  saveRDS(scores_df, paste0("regulation/data/TFs/ClusterBuster/results/",JASPAR_folder_name,"/subset_cbust/scores/scores_NPC_expr_df.rds"))
  return(scores_df)
}


