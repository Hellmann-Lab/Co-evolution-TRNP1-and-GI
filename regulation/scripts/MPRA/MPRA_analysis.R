# MPRA data analysis by Zane Kliesmete

#load libraries 
libs<-c("rBLAST","tidyverse","cowplot","GenomicRanges","data.table","reshape2", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","grid","gtable","xtable","rr2","Biostrings")
sapply(libs, require, character.only=T)


################
# LOAD DATA ####
################

#load helper functions 
source("regulation/scripts/MPRA/MPRA_helper_functions.R")

#load the upgraded mammal tree from Bininda-Emonds 2007
mammaltree<-read.tree("protein/trees/mammaltree.txt") 


#load the combined phenotype data
pheno_data<-read_rds("pheno_data/pheno_data.rds") %>%
  dplyr::select("species","brain_mass","GI") 

#load MPRA mpra_metadata, select only the relevant information
mpra_metadata<-read.delim("regulation/data/MPRA/input/complete_library_meta_information.txt") %>%
  #also exclude ancestral sequences and mouse-specific elements
  filter(!grepl("node|unique", SeqID)) %>%
  dplyr::select(SeqID, species, region, enhancer_number, enhancer_sequence, enhancer_length) %>%  
  dplyr::rename(enhancer=enhancer_number, sequence=enhancer_sequence,length=enhancer_length) %>%
  #adjust mismatching species names
  mutate(species=gsub("Mustela putorius furo","Mustela putorius", species),
         species=gsub("Equus ferus caballus","Equus ferus", species),
         species=gsub("Canis lupus familiaris","Canis lupus", species),
         species=gsub(" ","_", species))


#annotate mpra_metadata with clade information
mpra_metadata_tree<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% mpra_metadata$species])
length(mpra_metadata_tree$tip.label) #75, all there.
ggtree(mpra_metadata_tree)+geom_tiplab(size=3)+xlim(0,200) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
write.table(mpra_metadata_tree$tip.label, "regulation/trees/all_mpra_species.txt", col.names = F, row.names = F, quote = F)


primates_mpra_metadata<-extract.clade(mpra_metadata_tree, node=91)
rodents_mpra_metadata<-extract.clade(mpra_metadata_tree, node=83)
carnivores_mpra_metadata<-extract.clade(mpra_metadata_tree, node=138)

#primate clades
OWMs_mpra_metadata<-extract.clade(mpra_metadata_tree, node=95)
apes_mpra_metadata<-extract.clade(mpra_metadata_tree, node=111)
NWMs_mpra_metadata<-extract.clade(mpra_metadata_tree, node=118)



mpra_metadata<-mpra_metadata %>%
  mutate(clade=case_when(species %in% primates_mpra_metadata$tip.label ~ "Primate",
                         species %in% rodents_mpra_metadata$tip.label ~ "Rodent",
                         species %in% carnivores_mpra_metadata$tip.label ~ "Carnivore",
                         T ~ "Other"),
         primate_clade= case_when(species %in% OWMs_mpra_metadata$tip.label ~ "Old World monkey",
                                  species %in% apes_mpra_metadata$tip.label ~ "Great ape",
                                  species %in% NWMs_mpra_metadata$tip.label ~ "New World monkey",
                                  T ~ "Other"))
saveRDS(mpra_metadata, "regulation/data/MPRA/input/mpra_metadata.rds")


#########################
# 0 PREPROCESS DATA #####
#########################

all_counts_MPRA <-read.delim("regulation/data/MPRA/input/all_counts_MPRA_long.txt") %>%
  #exclude unassigned barcodes, virus libraries and technical replicates 
  dplyr::filter(!grepl("unknown|virus",sample_biorep) & technrep!=3 & technrep!=1) %>%
  #also exclude ancestral sequences and mouse-specific elements
  dplyr::filter(!grepl("node|unique", SeqID)) %>%
  dplyr::select(-technrep) %>%
  dplyr::rename(sample_names=sample_biorep) %>%
  #annotate
  dplyr::mutate(cell_line=case_when(grepl("29",sample_names) ~ "human1",
                                  grepl("30", sample_names) ~ "human2",
                                  grepl("39", sample_names) ~ "macaque",
                                  grepl("plasmid", sample_names) ~ "plasmid"),
         biorep=case_when(grepl("_I$", sample_names) ~ 1,
                          grepl("_II$", sample_names) ~ 2,
                          grepl("_III$", sample_names) ~ 3)) 


# library size
lib_size<-all_counts_MPRA %>% group_by(sample_names) %>%
  summarise(sum=sum(counts)) %>%
  ungroup()
summary(lib_size$sum)

ggplot(lib_size,aes(x=sample_names, y=sum)) +
  geom_bar(stat="identity") +
  coord_flip()

# library complexity and chunk distribution across regions ####
#how many are there at all?
total_SeqID<-all_counts_MPRA %>%
  filter(grepl("plasmid", cell_line)) %>%
  filter(counts>0) %>%
  group_by(SeqID) %>%
  dplyr::summarise(n=length(SeqID)) %>%
  filter(n>0)

dim(total_SeqID)[1]/length(unique(all_counts_MPRA$SeqID)) # reliably captured 86% of chunks (4251)


#how are they distributed across the regions and species? (prepare the data for plotting)
chunk_distr<-data.frame(SeqID=unique(all_counts_MPRA$SeqID)) %>%
  left_join(mpra_metadata) %>%
  mutate(present=.$SeqID %in% total_SeqID$SeqID) %>%
  group_by(species, region, present) %>%
  summarise(n=length(species)) %>%
  spread(present, n, fill=0) %>%
  dplyr::rename(absent="FALSE", present="TRUE") %>%
  mutate(total=sum(c(absent, present)),
         fract_present=present/total) %>%
  left_join(mpra_metadata[,c("species","clade")])

chunk_distr$region<-droplevels(as.factor(chunk_distr$region))
chunk_distr$species<-droplevels(as.factor(chunk_distr$species))
chunk_distr$region<-factor(chunk_distr$region,
                                          levels=c("upstream1","upstream2","upstream3","exon1","intron","exon2","downstream"))

saveRDS(chunk_distr, "regulation/data/MPRA/output/tile_coverage_distr.rds")


 

#filtering ####
#exclude chunks that have been detected in 1/3 plasmid samples or less
keep_SeqID<-all_counts_MPRA %>%
  filter(grepl("plasmid", cell_line)) %>%
  filter(counts>0) %>% #normally:0 
  group_by(SeqID) %>%
  dplyr::summarise(n=length(SeqID)) %>%
  filter(n>1)
dim(keep_SeqID)[1]/length(unique(all_counts_MPRA$SeqID)) # keep 85% of chunks (4202)

# filter SeqID and calculate counts per milion for each library
filtered_counts_MPRA<-all_counts_MPRA %>%
  filter(SeqID %in% keep_SeqID$SeqID) %>%
  group_by(sample_names) %>%
  mutate(cpm=1e6*counts/sum(counts, na.rm=T), 
         log2_cpm=log2(cpm+0.01)) %>%
  ungroup()

saveRDS(filtered_counts_MPRA, "regulation/data/MPRA/output/filtered_counts_MPRA.rds")



# correlation across samples
cpms_MPRA_wide_log2 <- dcast(filtered_counts_MPRA, SeqID ~ cell_line+biorep, value.var = "log2_cpm", fill=0) %>%
  column_to_rownames("SeqID")

sampleCors <- cor(cpms_MPRA_wide_log2, method="pearson")

sampleCors_long<-as.data.frame(sampleCors) %>%
  rownames_to_column(var="sample2")%>%
  gather(sample,cor,2:ncol(.)) %>%
  dplyr::filter(sample!=sample2) %>%
  rowwise() %>%
  mutate(combi=paste(sort(c(sample,sample2)), collapse=",")) %>%
  dplyr::distinct(combi, .keep_all=T) %>%
  ungroup() %>%
  mutate(mean_cor=mean(cor), sd_cor=sd(cor)) %>%
  rowwise() %>%
  mutate(outlier=cor<(mean_cor-1.5*sd_cor))

#exclude macaque replicate 3 as it is a clear outlier
filtered_counts_MPRA<-filtered_counts_MPRA %>%
  filter(!(cell_line=="macaque" & biorep==3))
table(filtered_counts_MPRA$cell_line, filtered_counts_MPRA$biorep)




#####################################################################
# 1 CALCULATE NORMALIZED MEDIAN ACTIVITY PER CHUNK PER CELL LINE ####
#####################################################################

# first, calculate median input CPM from plasmid
input <-filtered_counts_MPRA %>%
  filter(grepl("plasmid", cell_line)) %>%
  dplyr::select(SeqID, cpm, biorep) %>%
  group_by(SeqID) %>%
  dplyr::summarise(median_input_cpm=median(cpm))

median_activity_MPRA<-filtered_counts_MPRA %>%
  filter(cell_line!="plasmid") %>%
  #calculate median CPM for each SeqID per cell line
  group_by(cell_line, SeqID) %>% 
  summarize(median_cpm=median(cpm)) %>%
  ungroup() %>%
  left_join(mpra_metadata) %>%
  left_join(input) %>%
  mutate(activity = median_cpm/median_input_cpm,
         log2_activity=log2(activity+0.01),
         activity_per_bp=activity/length,
         enhancer=ifelse(is.na(enhancer),1, enhancer))

dim(median_activity_MPRA)[1]/3


median_activity_MPRA<-median_activity_MPRA %>%
  dplyr::mutate(species=gsub(" ","_",species),
                enhancer=as.numeric(enhancer)) %>%
  #exclude regions that were accidentally included twice during library construction
  filter(!grepl("Saguinus_oedipus_trinity",SeqID))

saveRDS(median_activity_MPRA, "regulation/data/MPRA/output/median_activity_MPRA.rds")



##########################################
# 2 CALCULATE THE PER-REGION ACTIVITY ####
##########################################

system("bash regulation/scripts/MPRA/region_activity_calc/run_MPRA_summarization.sh", wait = T) 

#after all regions have run through, summarize
mpra_metadata<-readRDS("regulation/data/MPRA/input/mpra_metadata.rds")
pheno_data<-read_rds("pheno_data/pheno_data.rds") %>%
  dplyr::select("species","body_mass","brain_mass","EQ","GI") 

mpra_metadata_short<-mpra_metadata %>% 
  dplyr::select(species,clade,primate_clade) %>%
  dplyr::distinct(species, .keep_all=T) 


activity_overlap_summary<-data.frame(cell_line=as.character(),
                                species=as.character(),
                                sum_activity=as.numeric(),
                                length=as.numeric(), 
                                region=as.character(), 
                                log2_total_activity=as.numeric())

regions<-c("intron", "upstream3", "upstream2", "upstream1", "exon2", "downstream", "exon1")

for (i in regions){
  summarized_output<-readRDS(paste0("regulation/data/MPRA/output/summarized_activity/summarized_activity_",i,".rds"))
  activity_overlap_summary<-bind_rows(activity_overlap_summary, summarized_output %>% 
                                      mutate(region=paste0(i),log2_total_activity=log2(sum_activity+0.01)))
}

activity_overlap_summary<-activity_overlap_summary %>%
  left_join(mpra_metadata_short) 
  
saveRDS(activity_overlap_summary,"regulation/data/MPRA/output/activity_overlap_summary.rds")



#number of species per region with pheno data
activity_overlap_summary_brain<-activity_overlap_summary %>% 
  left_join(pheno_data) %>% 
  filter(!(is.na(brain_mass)))
table(activity_overlap_summary_brain$region, activity_overlap_summary_brain$cell_line)

activity_overlap_summary_GI<-activity_overlap_summary %>% left_join(pheno_data) %>% 
  filter(!(is.na(GI)))
table(activity_overlap_summary_GI$region, activity_overlap_summary_GI$cell_line)



#3 DO PGLS ####
cell_lines<-c("human1")
regions<-c("upstream1","upstream2","upstream3","exon1","intron","exon2","downstream")
phenos<-c("GI","brain_mass")
corStructs<-c("BM")
pheno_data<-readRDS("pheno_data/pheno_data.rds")

pgls_model_output_list<-PGLS_pheno_vs_MPRAactivity(activity_overlap_summary, mammaltree, cell_lines, regions, phenos,pheno_data,corStructs)



# LRT model selection on regions as predictors ####
pgls_LRT_output<-pgls_model_output_list[[2]] %>% bind_rows() %>%
  mutate(call=gsub("(gls[(]model = )","",call),
         call=gsub("\\,.*","",call),
         call=case_when(grepl("log2_total_activity",call) ~paste0(call,region,")"),
                        T ~ call),
         call=gsub("log2_total_activity","log2(",call),
         call=gsub("brain_mass","brain size",call),
         df=as.factor(df)) %>%
  dplyr::select(-pheno,-region,-cell_line,-CorStr)




#get the estimates for either
pgls_model_output<-pgls_model_output_list[[1]] %>% bind_rows() %>%
  filter(region!="(Intercept)") %>%
  mutate(model=gsub("brain_mass","brain size",model),
         model=gsub("~"," ~ ", model),
         df=as.factor(1)) %>%
  left_join(pgls_LRT_output[,c("call","L.Ratio")], by=c("model"="call")) %>%
  dplyr::select(model, Value, Std.Error,df,L.Ratio, mod.sel.lrt) %>%
  dplyr::rename(Model=model,`LRT p.value`=mod.sel.lrt) %>%
  dplyr::arrange(Model)
  
print(xtable(pgls_model_output, digits=c(1,1,2,3,1,2,3)),hline.after = c(-1,0,7,14),include.rownames=FALSE,file="regulation/data/MPRA/xtables/LRT_GI_brainW_allSp_coeff.tex")





# ZOOM INTO THE INTRON TREE ####
#intron activity in OWMs & apes across cell lines ####
pgls_hum1<-PGLSzoomin(clade_type = "primate_clade",clade_subset = c("Old World monkey","Great ape"), cells = "human1")
pgls_hum2<-PGLSzoomin(clade_type = "primate_clade",clade_subset = c("Old World monkey","Great ape"), cells = "human2")
pgls_mac<-PGLSzoomin(clade_type = "primate_clade",clade_subset = c("Old World monkey","Great ape"), cells = "macaque")


#LRTs
pgls_LRT_OWMs_apes<-bind_rows(pgls_hum1[[1]],pgls_hum2[[1]],pgls_mac[[1]]) %>%
  dplyr::select(-Species) %>%
  mutate(call=gsub("(gls[(]model = )","",call),
         call=gsub("\\,.*","",call),
         call=gsub("log2_total_activity","log2(intron)",call),
         df=as.factor(df))

#coefficients
pgls_intron_cell_lines<-bind_rows(pgls_hum1[[2]] %>% mutate(R2=pgls_hum1[[4]], `Cell line`="human1"),
                                       pgls_hum2[[2]] %>% mutate(R2=pgls_hum2[[4]], `Cell line`="human2"),
                                       pgls_mac[[2]] %>% mutate(R2=pgls_mac[[4]], `Cell line`="macaque")) %>%
  mutate(Predictor=gsub("log2_total_activity","log2(intron)",Predictor))


#combine both 
comb_lines<-pgls_intron_cell_lines %>% filter(Predictor!="(Intercept)") %>%
  dplyr::select(Value, Std.Error, `Cell line`) %>%
  right_join(pgls_LRT_OWMs_apes %>% filter(call!="log2(GI) ~ 1")) %>%
  mutate(df=as.factor(1)) %>%
  dplyr::select(call, Value, Std.Error, df, L.Ratio, `p-value`, `Cell line`) %>%
  dplyr::rename(Model=call, `LRT p-value`=`p-value`) 

print(xtable(comb_lines, digits=c(1,1,2,3,1,2,3,1)),include.rownames=FALSE,file="regulation/data/MPRA/xtables/LRT_GI_intron_lines_coeffs.tex")
