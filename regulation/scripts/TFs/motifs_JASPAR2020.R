# PREPARE JASPAR 2020 DATA

library(tidyverse)

# 1 DOWNLOAD AND TRANSPOSE,COMBINE PWMs ####

setwd("/data/share/htp/TRNP1/paper_data/")
system("fold=regulation/data/TFs/JASPAR_2020; mkdir -p $fold; mkdir -p $fold/JASPAR_collection; 
       mkdir -p $fold/JASPAR_collection_transposed; mkdir -p $fold/JASPAR_collection_transposed/all;
       cd $fold/JASPAR_collection; wget http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.zip; 
       unzip JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.zip; 
       rm JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.zip")
#also download PWM clustering data
system('sbatch --wrap="wget http://jaspar.genereg.net/static/clustering/JASPAR_2020_clusters/vertebrates/radial_trees/JASPAR_2020_matrix_clustering_vertebrates_archive.zip -P /data/share/htp/TRNP1/paper_data/regulation/data/TFs/JASPAR_2020/matrix_clustering"')       



allfiles <- list.files(path="regulation/data/TFs/JASPAR_2020/JASPAR_collection", full.names=T) #746 motifs

#need to transpose each Position Frequency Matrix: currently bases are the rows, I need them as columns for ClusterBuster
tf_list<-list()
for (i in 1:length(allfiles)){
  jaspar_example<-read.table(paste0(allfiles[i]), skip=1)
  jaspar_transposed<-t(jaspar_example[,-c(1,2,length(jaspar_example))])
  jaspar_name<-read.table(paste0(allfiles[i]), nrows=1, sep=";")
  jaspar_name$V1<-gsub("\t"," ", jaspar_name$V1)
  jaspar_name$V1<-gsub(">","", jaspar_name$V1)
  tf_list[[i]]<-jaspar_name$V1
  write.table(jaspar_transposed, paste0("regulation/data/TFs/JASPAR_2020/JASPAR_collection_transposed/all/",jaspar_name$V1), row.names = F, col.names = F)
}

system('mkdir regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster; cd regulation/data/TFs/JASPAR_2020/JASPAR_collection_transposed/all/; for i in *; do echo ">$i";cat "$i"; done >  ../../ready_for_ClusterBuster/TF_psfm_jaspar2020.txt')

tf_df<-plyr::ldply(tf_list)
tf_df<-tidyr::separate(tf_df,V1, c("ID", "SYMBOL"), sep=" ")
saveRDS(tf_df,"regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/TF_IDs_jaspar2020.rds")




# 2 GET ENSEMBL IDs ####
#IDs are based on standardized Entrez symbols. Concatenated names (with ::) indicate hetero-dimers; versions are given in case different splice forms of the same protein bind to different motifs. http://jaspar.genereg.net/docs/
#The JASPAR CORE contains a curated, non-redundant set of profiles, derived from published and experimentally defined transcription factor binding sites for eukaryotes. It should be used, when seeking models for specific factors or structural classes, or if experimental evidence is paramount.

#JASPAR 2020 CORE TF IDs vs SYMBOLS 
tf_df<-readRDS("regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/TF_IDs_jaspar2020.rds") #746 TFs in total
tf_df<-mutate_all(tf_df, .funs=toupper)
tf_df<-tf_df %>%
  mutate(SYMBOL_noVar=gsub("[(]VAR.(2|3)[)]","", SYMBOL),
         SYMBOL1=word(SYMBOL_noVar, 1,1,"::"),
         SYMBOL2=word(SYMBOL_noVar, 2,2,"::"),
         SYMBOL3=word(SYMBOL_noVar, 3,3,"::"))

#translate from ENSEMBL to symbols (expressed_ensembl)
library(org.Hs.eg.db)
tf_df_ENSEMBL<-AnnotationDbi::select(org.Hs.eg.db,keys=unique_TFs,keytype="SYMBOL",columns=c("ENSEMBL","SYMBOL")) %>%
  left_join(tf_df, by=c("SYMBOL"="SYMBOL1")) %>%
  dplyr::rename(motif=ID, motif_SYMBOL=SYMBOL.y) %>%
  dplyr::select(SYMBOL, motif, ENSEMBL, motif_SYMBOL)
saveRDS(tf_df_ENSEMBL, "regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/TF_ENSEMBL_IDs_jaspar2020.rds")
