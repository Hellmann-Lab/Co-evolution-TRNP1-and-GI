#generate additional_genomes_full file

setwd("/data/share/htp/TRNP1/paper_data/")

dir<-"/data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/downl_genomes/"

genomes<-bind_rows(read.table("protein/data/NCBI_genomes.txt"),read.table("protein/data/ENSEMBL_genomes.txt"))
genomes$fa_name<-word(genomes$V2,-1,-1,"/")
genomes$path<-paste0(dir,genomes$V1)


ucsc_genomes<-read.table("protein/data/additional_genomes.txt") %>%
  filter(grepl("UCSC",V2)) %>%
  dplyr::rename(path=V2, fa_name=V3)

comb<-bind_rows(genomes[,c("V1","path","fa_name")],
                ucsc_genomes) %>%
  #mutate(fa_name=gsub("([.]gz)$","",fa_name)) %>%
  dplyr::rename(species=V1)

#for these species use the resequenced contigs
exclude_for_TRNP1<-c("Cercopithecus_mona", "Mandrillus_sphinx", "Papio_anubis", "Mustela_putorius", "Homo_sapiens")

write.table(comb %>% filter(!species %in% exclude_for_TRNP1), "protein/data/additional_genomes_full.txt", col.names = F, row.names = F, quote = F)




#which ones have GI info?
#also, still need to add human and ferret
pheno_data_new<-readRDS("/data/share/htp/TRNP1/paper_data/data_tables/pheno_data/pheno_data_new.rds")


#make a selection of the 31 species for searching other protein sequences in those. 

comb %>% filter(species %in% (pheno_data_new %>% filter(!is.na(GI) & GI!="NaN"))$species | 
                  species=="Cercopithecus_mona") %>%
  write.table("protein/data/coevol_species_genomes.txt", col.names = F, row.names = F, quote = F)

#once just NCBI for testing
comb %>% filter(species %in% (pheno_data_new %>% filter(!is.na(GI) & GI!="NaN"))$species | 
                  species=="Cercopithecus_mona") %>%
  filter(grepl("genomic.fna.gz",fa_name)) %>%
  write.table("protein/data/coevol_NCBI_genomes.txt", col.names = F, row.names = F, quote = F)




