
#selecting candidate genes
#we want to have ~ <1000 AA long proteins, preferentially from 1 exon, present in 30 species, overlap with human sequence at least 50%

setwd("protein/other_protein_alignments/blat_output_full/")
overlap_list<-list.files(recursive = T, pattern = ".rds")

rbb_list<-lapply(overlap_list, function(x){readRDS(x)})
names(rbb_list)<-overlap_list
rbb_df<-bind_rows(rbb_list, .id="species") %>%
  mutate(species=word(species,1,1,sep="/"))

# select IDs that have only one exon in total
# filter rbb for these
# select length <1000
# check that overlap across all species is high enough (>=0.5) for each candidate


setwd("protein/other_protein_alignments/")
#read in the ccds ( retrieved from https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/, command in README )
CCDS<-read_delim("CCDS/CCDS.current.txt", delim="\t")
#select CCDS with 1 exon
library(stringr)
CCDS_filt<- CCDS %>% mutate(lines=str_count(cds_locations,"-")) %>% filter(lines==1)


rbb_df_filt<-rbb_df %>%
  dplyr::filter(name.ref %in% CCDS_filt$ccds_id, width<1000) %>%
  dplyr::filter(percOverlap_ref>0.50) %>%
  dplyr::filter(!seqnames %in% c("chrX", "chrY")) %>%
  dplyr::group_by(name.ref) %>% 
  dplyr::mutate(n_species=length(unique(species))) %>% 
  dplyr::ungroup() %>%
  dplyr::filter(n_species==29)


#subset CCDS
CCDS<-readAAStringSet("CCDS/CCDS_protein.current.faa.gz")
CCDS_322<-CCDS[word(names(CCDS),1,1,"[|]") %in% rbb_df_filt$name.rblat]
writeXStringSet(CCDS_322, "CCDS_protein.current.322.fa")



