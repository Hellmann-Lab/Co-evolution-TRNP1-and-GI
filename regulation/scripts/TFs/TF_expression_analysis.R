# TF expression and binding to the DHS, particularly intron

#load libraries 
libs<-c("tidyverse","cowplot","GenomicRanges","data.table","reshape2", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","org.Hs.eg.db","DESeq2","flashClust","topGO", "biomaRt","vsn","GGally","pheatmap","RColorBrewer")
sapply(libs, require, character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/")


#load MPRA data 
activity_overlap_summary<-readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds")

#load the upgraded mammal tree from Bininda-Emonds 2007
mammaltree<-read.tree("protein/trees/mammaltree.txt") 


#load helper functions 
source("regulation/scripts/MPRA/MPRA_helper_functions.R")


#################################################################
# 1 Infer expressed TFs in the transduced NPCs used for MPRA ####
#################################################################

#prepare the data ####
npc_counts<-readRDS("regulation/data/TFs/expression/input/MPRNAseq_NPC.dgecounts.rds")

#select intron & exon counts
umi_counts<-as.data.frame(as.matrix(npc_counts$umicount$inex$all))
annot<-read.table("regulation/data/TFs/expression/input/MPRA_info.txt", header = T) %>%
  mutate(individual=case_when(clone == "29B5" ~ "human1",
                              clone == "30B1" ~ "human2",
                              clone == "39B2" ~ "macaque"))

#filtering####
# filter, visualize
annot<-annot[annot$XC %in% colnames(umi_counts),]
rownames(annot)<-annot$XC
umi_counts <- umi_counts[,match(rownames(annot),colnames(umi_counts))]
colnames(umi_counts) == rownames(annot)
umi_counts <- umi_counts[rowSums(umi_counts)>0,]
summary(colSums(umi_counts))


#plot the number of detected genes in each sample
annot<- annot %>%
  left_join(data.frame(XC=colnames(umi_counts), 
                       expressed_genes=apply(umi_counts, 2, function(x){sum(x>0)}),
                       total_umi_counts=apply(umi_counts,2,sum)))

ggplot(annot, aes(x=sample, y=expressed_genes, fill=total_umi_counts))+
  geom_bar(stat="identity")+
  coord_flip()+
  facet_grid(clone~., scales = "free_y")+
  scale_fill_continuous(labels = function(x) format(x, scientific = TRUE))
ggsave("regulation/data/TFs/expression/figures/1_expressed_genes.pdf", height=4, width=6)

#very similar number of expressed genes and total counts across samples (only 3x variation in the latter) 


#dropouts and average expression
#Dropouts
quick_zero <- rowMeans(umi_counts==0)
s<-summary(quick_zero)
boxplot(quick_zero)
dropout_umi = sum(umi_counts==0)/(nrow(umi_counts)*ncol(umi_counts))
s
#median: 0.0, mean: 0.26

#average expression
rwm <- rowMeans(umi_counts)
hist(rwm[rwm<100], freq=FALSE, xlab="Average UMI count per gene", col="grey80", cex.lab=1, cex.axis=0.8, main="", breaks = 50, xaxt="n")
axis(side=1, at=seq(0,100, 10), labels=seq(0,100, 10))
fr2<-rwm[rwm<=2]
fr2_2<-length(fr2)/length(rwm)
#38.4% genes show <= 2 UMI counts in average
fr7<-rwm[rwm<=7]
fr7_2<-length(fr7)/length(rwm)
#53% genes show <= 7 UMI counts in average


#filtering: keep only genes that are expressed in at least 70% of samples (so in 7/9) 
freq_umicounts <- apply(umi_counts > 0, 1, mean, na.rm = T)
freq_expressed <- 0.7
expressed_genes <- freq_umicounts > freq_expressed
table(expressed_genes)
umi_counts_filt <- umi_counts[expressed_genes,]

#also exclude genes with lower than 7 average expression
include_genes<-rownames(umi_counts_filt) %in% names(fr7)
table(include_genes)
umi_counts_filt<-umi_counts_filt[!include_genes,]
dim(umi_counts_filt) 
#17306  genes 

#is TRNP1 there?
umi_counts_filt[rownames(umi_counts_filt)=="ENSG00000253368",] %>%
  gather(XC, TRNP1_counts) %>%
  left_join(annot)
#yes, a bit higher in human samples


#find expressed TFs####
tf_df_ENSEMBL<-readRDS("regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/TF_ENSEMBL_IDs_jaspar2020.rds")
umi_counts_filt_ENS<-rownames_to_column(umi_counts_filt,"ENSEMBL")

NPC_all_TF_counts<-umi_counts_filt_ENS %>% 
  filter(umi_counts_filt_ENS$ENSEMBL %in% unique(tf_df_ENSEMBL$ENSEMBL)) %>%
  column_to_rownames("ENSEMBL") #392 TF ENSEMBL gene IDs


#select the expressed TF motifs
NPC_expr_tf_df<-tf_df_ENSEMBL %>%
  filter(ENSEMBL %in% umi_counts_filt_ENS$ENSEMBL) #462
length(unique(NPC_expr_tf_df$motif)) #462
length(unique(NPC_expr_tf_df$ENSEMBL)) #392

saveRDS(NPC_expr_tf_df,"regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/NPC_expressed_TF_IDs_jaspar2020.rds")


#normalization####
# normalize counts of NPCs 
rownames(annot)<-annot$XC
colnames(umi_counts_filt)==rownames(annot)
dds<-DESeq2::DESeqDataSetFromMatrix(countData=as.matrix(umi_counts_filt), 
                            colData=data.frame(clone=annot$clone),
                            design=~0+clone)
dds<-DESeq2::DESeq(dds)

#dispersion estimates
pdf("regulation/data/TFs/expression/figures/2_dispersion.pdf", height=4.5, width=6)
plotDispEsts(dds)
dev.off()

# do variance stabilizing transformation, export normalized counts
vsd<-DESeq2::varianceStabilizingTransformation(dds)
vsdMat<-assay(vsd)
saveRDS(vsdMat,"regulation/data/TFs/expression/vsdMat.rds")

#check normalization
normBoxplot<-function(norm_table){
  k<-reshape2::melt(as.matrix(norm_table))
  colnames(k)<-c("Ensembl","XC","var_stabilized_counts")
  plotdat_k <- left_join(k, annot, by="XC")
  p2 <- ggplot(data = plotdat_k, aes(x=sample, y=var_stabilized_counts)) + geom_boxplot() + coord_flip()+ theme(axis.text.y = element_text(size=6))
  plot(p2)
}

pdf("regulation/data/TFs/expression/figures/3_meanSdPlot.pdf", height=3, width=5)
meanSdPlot(vsdMat)
dev.off()

normBoxplot(vsdMat)
ggsave("regulation/data/TFs/expression/figures/4_normalization_boxplot.pdf", height=4, width=4)



# PCA on mv 500 genes ####
mat <- vsdMat
rowVar <- apply(mat, 1, var)
mv500 <- order(rowVar, decreasing = T)[1:500] # transpose matrix and do analysis
mat <- t(mat[mv500, ])
pc <- prcomp(mat, scale = T)
pc.sum <- summary(pc)$importance 
pc.sum[, 1:5]


#PCA_500mv genes
varExp <- round(pc.sum[2, ] * 100, 2)
pcs <- data.frame(pc$x, clone=annot$clone, species=annot$species)

p12_point_bb <- ggplot(pcs, aes(x = PC1, y = PC2, col =clone, shape=species)) +geom_point(size=4, alpha=0.7)+ xlab(paste("PC1 (", varExp[1], "%)")) + ylab(paste("PC2 (",varExp[2], "%)"))
p12_point_bb
ggsave("regulation/data/TFs/expression/figures5_PCA_500mv.pdf", height=3.5, width=5)





#TRNP1 expression ####

norm_counts_trnp1<-vsdMat[rownames(vsdMat)=="ENSG00000253368",] %>%
  as.data.frame(.) %>%
  rownames_to_column("XC") %>%
  left_join(annot) %>%
  dplyr::rename(TRNP1_norm_expr=".")

ggplot(norm_counts_trnp1, aes(x=clone, y=TRNP1_norm_expr, color=species))+
  geom_point(size=3.5)+
  ylab("variance-stabilized expression of TRNP1")

ggsave("regulation/data/TFs/expression/figures6_TRNP1_expression.pdf", height=4, width=6)



# DE results
# TRNP1 ENSG00000253368
hum1_vs_hum2<-results(dds,contrast=c("clone","29B5","30B1")) 
hum1_vs_hum2_table<-as.data.frame(hum1_vs_hum2) %>%
  rownames_to_column(var="ENSEMBL") #TRNP1 nonsignificant

hum1_vs_mac<-results(dds,contrast=c("clone","29B5","39B2"))
hum1_vs_mac_table<-as.data.frame(hum1_vs_mac) %>%
  rownames_to_column(var="ENSEMBL")


hum2_vs_mac<-results(dds,contrast=c("clone","30B1","39B2"))
hum2_vs_mac_table<-as.data.frame(hum2_vs_mac) %>%
  rownames_to_column(var="ENSEMBL")

#gather results on TRNP1

DE_comb<-bind_rows(hum1_vs_hum2_table %>% arrange(padj) %>% 
                     mutate(contrast="human1 vs human2", rank = 1:nrow(.)),
                   hum1_vs_mac_table %>% arrange(padj) %>% 
                     mutate(contrast="human1 vs macaque", rank = 1:nrow(.)),
                   hum2_vs_mac_table %>% arrange(padj) %>% 
                     mutate(contrast="human2 vs macaque", rank = 1:nrow(.))) %>%
  filter(ENSEMBL=="ENSG00000253368")




vsdMat_long<-vsdMat %>%
  as.data.frame(.) %>%
  rownames_to_column("ENSEMBL") %>%
  pivot_longer(-ENSEMBL, names_to = "XC", values_to = "vsd_expr") %>%
  left_join(annot) %>%
  group_by(XC) %>%
  arrange(-vsd_expr) %>%
  mutate(rank = 1:length(XC)) %>%
  ungroup()





# 2 Intron-enriched TF expression ####
# Used ClusterBuster for motif finding on the intron sequence. Considering only the expressed TFs in our NPCs and OWM+ape sequences 

NPC_expr_tf_df<-readRDS("regulation/data/TFs/JASPAR_2020/ready_for_ClusterBuster/NPC_expressed_TF_IDs_jaspar2020.rds")
unique_motifs_vec<-unique(NPC_expr_tf_df$motif)

system("mkdir -p regulation/data/TFs/ClusterBuster/results")
system("mkdir -p regulation/data/TFs/ClusterBuster/results/slurms")
run_ctrain_cbust(tf_motif_vector = unique_motifs_vec,JASPAR_folder_name =  "OWMs_apes_intron_weighted",query_fasta = "OWMs_apes_intron_correct_withGI.fa", add_weights = TRUE,cluster_cutoff = 5, motif_cutoff = 3) #takes around 7 min to run ctrain and cbust.


#ctrain weights
ctrain_weights<-read.table("regulation/data/TFs/ClusterBuster/results/OWMs_apes_intron_weighted/ctrain_output.txt", skip = 19, nrows = length(unique_motifs_vec), sep=",")
ctrain_weights<-ctrain_weights %>%
  tidyr::separate(V1, into=c("score","motif"), sep="MA",extra="merge") %>%
  tidyr::separate(motif, into=c("motif","motif_SYMBOL"), sep="\\s") %>%
  mutate(motif=paste0("MA",motif),
         score=as.numeric(score),
         motif_SYMBOL=toupper(motif_SYMBOL),
         SYMBOL_clean=gsub("::|-","_",motif_SYMBOL),
         SYMBOL_clean=gsub("[(]VAR.2[)]","_VAR2",SYMBOL_clean)) %>%
  left_join(NPC_expr_tf_df, by=c("motif", "motif_SYMBOL"))

intron_motifs_1_OWMs_expr<-ctrain_weights %>% 
  filter(score>=1) %>%
  left_join(vsdMat %>% as.data.frame() %>%rownames_to_column("ENSEMBL"))


#visualize their expression (and add TRNP1)
norm_counts_trnp1_2<-norm_counts_trnp1 %>%
  mutate(SYMBOL="TRNP1",
         ENSEMBL="ENSG00000253368") %>%
  dplyr::rename(norm_expr=TRNP1_norm_expr)

intron_tf_expr_NPCs_long<-intron_motifs_1_OWMs_expr %>%
  gather(XC, norm_expr, 7:ncol(intron_motifs_1_OWMs_expr)) %>%
  left_join(annot) %>%
  bind_rows(norm_counts_trnp1_2)
saveRDS(intron_tf_expr_NPCs_long,"regulation/data/TFs/expression/vsdMat_22TFs.rds")



# topGO analysis on the TFs ####
cbust_weights_for_topGO<-ctrain_weights %>%
  drop_na() %>%
  group_by(ENSEMBL) %>%  
  top_n(n=1,wt=score) %>%
  ungroup() #392 genes

length(unique(cbust_weights_for_topGO$ENSEMBL)) #392 

#this function does topGO enrichment analysis for the genes that have a ctrain score larger than the specified cutoff
#and returns a list containing 1) initialized topGO object, 2) table containing the top 30 enriched categories
topGO_ctrain_weights1<-do_topGO_motifs_ENSEMBL(cbust_weights_for_topGO, cutoff=1,nodesize = 20)
topGO_ctrain_weights1_table<-topGO_ctrain_weights1[[2]] 

topGO_ctrain_weights1_table_sign<-topGO_ctrain_weights1_table %>%
  filter(Fisher<0.05) %>%
  mutate(`Signif/Expected` = round(Significant/Expected, digits=2)) %>%
  dplyr::rename(`Fisher's P`=Fisher)

topGO_ctrain_weights1_table_sign$Term[topGO_ctrain_weights1_table_sign$GO.ID=="GO:0042127" & topGO_ctrain_weights1_table_sign$Term=="regulation of cell population proliferat..."]<-"regulation of cell population proliferation"

topGO_ctrain_weights1_table_sign$Term[topGO_ctrain_weights1_table_sign$GO.ID=="GO:0008285" & topGO_ctrain_weights1_table_sign$Term=="negative regulation of cell population p..."]<-"negative regulation of cell population proliferation"

topGO_ctrain_weights1_table_sign$Term[topGO_ctrain_weights1_table_sign$GO.ID=="GO:1901615" & topGO_ctrain_weights1_table_sign$Term=="organic hydroxy compound metabolic proce..."]<-"organic hydroxy compound metabolic process"


saveRDS(topGO_ctrain_weights1_table_sign, "regulation/data/TFs/expression/topGO_table_ctrain_weights1.rds")
print(xtable(topGO_ctrain_weights1_table_sign[,1:6]),include.rownames=FALSE,file="regulation/data/TFs/xtables/topGO_ctrain_weights1_05.tex")



#get gene names from the top 3 enriched categories
top10_genes_OWMs<-genesInTerm(topGO_ctrain_weights1[[1]], topGO_ctrain_weights1_table_sign$GO.ID[1:4]) #%>%
  
all_g<-unique(unlist(top10_genes_OWMs))
df_x<-plyr::ldply(top10_genes_OWMs,function(s){data.frame(unlist(s))}) %>%
  dplyr::rename(ENSEMBL="unlist.s.", GO.ID=".id") %>%
  left_join(cbust_weights_for_topGO %>% filter(score>=1)) %>%
  left_join(topGO_ctrain_weights1_table) %>%
  filter(!is.na(SYMBOL))

saveRDS(df_x, "regulation/data/TFs/expression/topGO_genes_ctrain_weights1.rds")



#summarise the full info
intron_22TFs_summarized<-ctrain_weights %>% 
  filter(score>=1) %>%
  mutate(GO_term=case_when(SYMBOL %in% df_x$SYMBOL[df_x$GO.ID=="GO:0042127"] & SYMBOL %in% df_x$SYMBOL[df_x$GO.ID=="GO:0010817"] ~ "hormone levels and\ncell proliferation",
                           SYMBOL %in% df_x$SYMBOL[df_x$GO.ID=="GO:0042127"] ~ "cell proliferation",
                           SYMBOL %in% SYMBOL %in% df_x$SYMBOL[df_x$GO.ID=="GO:0010817"] ~ "regulation of hormone levels",
                           T ~ "none"),
         activity_predictor=case_when(SYMBOL %in% c("CTCF","ZBTB26","SOX8") ~ "yes",
                                      T ~ "no"),
         GI_predictor= case_when(SYMBOL %in% "CTCF" ~ "yes",
                                 T ~ "no"))

saveRDS(intron_22TFs_summarized, "regulation/data/TFs/expression/summarized_info_top22TFs_ctrain.rds")






#get the motif distances####
#system("cd regulation/data/TFs/JASPAR_2020/matrix_clustering/; unzip JASPAR_2020_matrix_clustering_vertebrates_archive.zip 'JASPAR_2020_matrix_clustering_vertebrates_tables/distance_table.tab' -d  .")

motif_distances<-read.table("regulation/data/TFs/JASPAR_2020/matrix_clustering/JASPAR_2020_matrix_clustering_vertebrates_tables/distance_table.tab") 
colnames(motif_distances)<-word(colnames(motif_distances), start=5, end=6, sep="_")
rownames(motif_distances)<-word(rownames(motif_distances), start=5, end=6, sep="_")
dim(motif_distances) #all motifs present, good.
#adjust names
colnames(motif_distances)<-gsub("_",".",colnames(motif_distances))
rownames(motif_distances)<-gsub("_",".",rownames(motif_distances))


#select the relevant 22
motif_distances_22<-motif_distances
motif_distances_22<-motif_distances_22[,match(intron_motifs_1_OWMs$motif, colnames(motif_distances_22))]
motif_distances_22<-motif_distances_22[match(intron_motifs_1_OWMs$motif, rownames(motif_distances_22)),]

#change to symbols
table(intron_motifs_1_OWMs$motif==rownames(motif_distances_22))
table(intron_motifs_1_OWMs$motif==colnames(motif_distances_22))

motif_distances_22_symb<-motif_distances_22
rownames(motif_distances_22_symb)<-intron_motifs_1_OWMs$SYMBOL_clean
colnames(motif_distances_22_symb)<-intron_motifs_1_OWMs$SYMBOL_clean
saveRDS(motif_distances_22_symb, "regulation/data/TFs/ClusterBuster/results/motif_distances_22_symb.rds")
