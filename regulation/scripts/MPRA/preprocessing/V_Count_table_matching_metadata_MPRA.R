### Script to read in many tables and combine them
library(tidyverse)
library(readr)
library(stringr)
library(reshape2)

#adjust the path to the output of script II_filtering.sh!!
files<-list.files(path="/data/share/htp/TRNP1/MPRA_full_300919/MPRA_data/deML/demult/count_summaries/counts/", pattern = ".txt")

meta_data <- read.delim2("regulation/data/MPRA/input/complete_library_meta_information.txt")


revcomp = function(x){
  comp = chartr('ACGT', 'TGCA', x)
  revcomp = sapply(comp, FUN = function(y){
    element = unlist(strsplit(y, split =''))
    element.rev = rev(element)
    return(paste(element.rev, collapse =''))
  })
  return(revcomp)
}

meta_data$barcode_revcomp<-revcomp(meta_data$barcode)


cnt_mtrx_MPRA<-data.frame(SeqID=as.character(meta_data$SeqID))

for (i in 1:length(files)){
  tmp<-read.table(paste0(files[i]))
  colnames(tmp)<- c("read_count", "barcode")
  tmp<-mutate(tmp,lib = str_split_fixed(files[i], "_", n=6)[,4])
  tmp<-mutate(tmp,rep = str_split_fixed(str_split_fixed(files[i], "_r1", n=2)[,1],"_",n=5)[,5])
  tmp<-mutate(tmp,lane = str_split_fixed(files[i], "_", n=6)[,3])
  tmp2<-left_join(meta_data, tmp, by=c("barcode_revcomp"="barcode")) %>%
    dplyr::select(SeqID,read_count)
  colnames(tmp2)[2]<-paste0(tmp$lib,"_",tmp$rep,"_",tmp$lane)[1]
  cnt_mtrx_MPRA<-left_join(cnt_mtrx_MPRA,tmp2)
}

colnames(cnt_mtrx_MPRA)

head(cnt_mtrx_MPRA)
cnt_mtrx_MPRA[is.na(cnt_mtrx_MPRA)]<-0


####################################
# combine counts from all lanes ####
####################################

cnt_mtrx_MPRA_long<-cnt_mtrx_MPRA %>% gather(sample,counts,2:55) %>%
  mutate(lane=word(sample,3,3,"_"),
         sample_biorep=word(sample,1,2,"_"),
         technrep=case_when(grepl("39B2_II_lane",sample) ~1,
                            grepl("39B2_II_2", sample) ~ 2,
                            grepl("39B2_II_3",sample) ~3))
cnt_mtrx_MPRA_long$technrep[is.na(cnt_mtrx_MPRA_long$technrep)]<-0

cnt_mtrx_MPRA_long<-cnt_mtrx_MPRA_long %>%
  group_by(sample_biorep,technrep, SeqID) %>%
  summarise(counts=sum(counts)) %>%
  ungroup() %>%
  filter(!grepl("unknown",sample_biorep))

table(cnt_mtrx_MPRA_long$sample_biorep,cnt_mtrx_MPRA_long$technrep)

write.table(cnt_mtrx_MPRA_long,"regulation/data/MPRA/input/all_counts_MPRA_long.txt",sep="\t",quote = F, row.names = T,col.names = T)

