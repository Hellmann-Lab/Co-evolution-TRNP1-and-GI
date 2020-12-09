# This script contains TRNP1 protein-coding orthologue search by Zane Kliesmete

setwd("/data/share/htp/TRNP1/paper_data/")

libs <- c("IRanges", "Biostrings", "tidyverse", "reshape2", "reshape", "GenomicRanges", "data.table", "phangorn", "ape", "ggrepel", "ggtree", "caper", "phytools", "rBLAST", "stringr","xtable", "geiger","nlme","phytools")

sapply(libs, require, character.only=T)


# load resequenced promoter-exon1 region from 18 species ####
# data generation in /data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/scripts/Rscripts/processingL_final.R by Lucas Wange
reseq.prom_exon1.fas<-readDNAStringSet('protein/fastas/resequenced_promexon1_primates.fa')


#load human TRNP1 protein sequence####
protein.fa <- readAAStringSet("protein/fastas/trnp1_protein_seq_hum.fasta")


#create a blast function
do_rblast<-function(subject,query=protein.fa,type="tblastn", min_perc_identity=70, min_alignment_length=190){
  dir<-tempdir()
  writeXStringSet(subject, filepath = file.path(dir,"seqs.fasta"))
  db<-makeblastdb(file.path(dir, "seqs.fasta"), dbtype = "nucl")
  res<-blast(file.path(dir, "seqs.fasta"),type=paste0(type))
  #switch off masking (normally masks overrepresented or low-complexity sequences), slightly lower the gap penalties
  blast_output<-stats::predict(res,query, BLAST_args="-soft_masking false -seg no")
  blast_output_filt<-blast_output %>% 
    filter(Perc.Ident>min_perc_identity & Alignment.Length>min_alignment_length) %>%
    group_by(SubjectID) %>% top_n(1, wt=Alignment.Length) %>% ungroup() %>%
    dplyr::distinct(SubjectID, .keep_all=T)
  
  #select the relevant sequences
  subject_seq<-subject[names(subject) %in% blast_output_filt$SubjectID]
  subject_seq<-subject_seq[match(blast_output_filt$SubjectID, names(subject_seq))]
  
  #cut the region and do revCompl where necessary
  subject_cut<-subject_seq
  for (i in 1:nrow(blast_output_filt)){
    if (blast_output_filt$S.end[i]<blast_output_filt$S.start[i]){
     subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.end[i], end=blast_output_filt$S.start[i])
     subject_cut[i]<-Biostrings::reverseComplement(subject_cut[i])
   } else {
      subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.start[i], end=blast_output_filt$S.end[i])
   }}
  
  return(list(subject_cut,blast_output_filt)) 
}

#resequenced ####
#blast against the human protein
blast.reseq.promexon1.out.all<-do_rblast(reseq.prom_exon1.fas) 
reseq.coding.fas<-blast.reseq.promexon1.out.all[[1]]
names(reseq.coding.fas)<-word(names(reseq.coding.fas),1,2,"_")
blast.reseq.promexon1.out<-blast.reseq.promexon1.out.all[[2]]
(reseq.coding.fas, 'protein/fastas/resequenced_coding_TRNP1.fa')


#NCBI####
# additional species from NCBI --> already pre-selected regions: blast separately
NCBI_Delphinapterus_leucas_TRNP1.fas<-readDNAStringSet("protein/fastas/TRNP1_coding/NCBI_Delphinapterus_leucas_TRNP1.fa")
names(NCBI_Delphinapterus_leucas_TRNP1.fas)<-"Delphinapterus_leucas"

NCBI_Odocoileus_virginianus_TRNP1.fas<-readDNAStringSet("protein/fastas/TRNP1_coding/NCBI_Odocoileus_virginianus_TRNP1.fa")
names(NCBI_Odocoileus_virginianus_TRNP1.fas)<-"Odocoileus_virginianus"

NCBI_species_TRNP1.fas<-c(NCBI_Delphinapterus_leucas_TRNP1.fas, NCBI_Odocoileus_virginianus_TRNP1.fas)
#blast
blast.NCBI.out.all<-do_rblast(NCBI_species_TRNP1.fas)
NCBI.coding.fas<-blast.NCBI.out.all[[1]]
blast.NCBI.out<-blast.NCBI.out.all[[2]]  #looking great




# ferret from Dr.  Miriam Esgleas (Mustela_putorius)####
mputorius.coding.fa<-DNAStringSet("atgccgggctgccgcatcagcgcctgcggcccgggggcccaggaagggaccgcggaaccggggtccccgccgccgccgccccgggagctcgtgtcgtcccctcagcccccgcccccatctccgaccttgactccgaccccggcttcggtctcggcgcccgccgactcagccccggcgtgggcgggctcggcagaggggcaggagctgcagcgctggcgccagggcgctaacgggggcgcgggggctaccgcgccggcagggggcgcggcggcggcgggggcagccgggggccgagcgctagagctggcggaagcgcgccggcgactgctggaggtggagggccgcaggcgcctggtgtcggagctggagagccgtgtgctgcagctgcaccgcgtcttcttggcggccgagctgcgcctggcgcaccgcgccgagagcctgggccgcctcagcggcggcgtggcttttgccaaaaatggctttgccagggatggagccaggttacaacctttaactggccccaagagggcagagtggcaaggtcacatgccagcactggacagaactaggtcccaggcctgtgtgtgtcttccccatggtgtctccttcctgcccaggggcaggggcta")
names(mputorius.coding.fa)<-"Mustela_putorius"

blast.ferret.out<-do_rblast(mputorius.coding.fa, min_perc_identity = 50, min_alignment_length = 150)
ferret.coding.fas<- subseq(mputorius.coding.fa, start=blast.ferret.out[[2]]$S.start, end=blast.ferret.out[[2]]$S.end )




#more seqs from various sources (these have been blasted already) ####
more.coding.files<-list.files("protein/fastas/TRNP1_coding/coding_seqs_sep2_full", full.names =T, pattern = ".fa")
more.coding.fas.list <- lapply(more.coding.files, function(x) readDNAStringSet(x))
more.coding.fas<-DNAStringSet(do.call(c, unname(unlist(more.coding.fas.list))))
more.coding.fas<-more.coding.fas[!names(more.coding.fas) %in% names(NCBI.coding.fas)]

#combine all ####
trnp1.coding.all.fa<-c(reseq.coding.fas,
                       NCBI.coding.fas,
                       ferret.coding.fas,
                       more.coding.fas)

writeXStringSet(trnp1.coding.all.fa, filepath='protein/fastas/TRNP1_coding_seqs_45sp_full.fa', format = 'fasta', append=F)







#adjust the phylogenetic tree####
#load the upgraded mammal tree from Bininda-Emonds 2007
mammaltree<-read.tree("protein/trees/mammaltree.txt")
mammaltree2<-mammaltree
mammaltree2$tip.label<-gsub("Propithecus_verreauxi", "Propithecus_coquereli", mammaltree2$tip.label )
mammaltree2$tip.label<-gsub("Cercocebus_galeritus", "Cercocebus_atys", mammaltree2$tip.label )
mammaltree2$tip.label<-gsub("Odocoileus_virginianus_texanus", "Odocoileus_virginianus", mammaltree2$tip.label )
mammaltree2$tip.label<-gsub("Procolobus_badius", "Piliocolobus_tephrosceles", mammaltree2$tip.label )


tree.exon1.coding.all<-drop.tip(mammaltree2, mammaltree2$tip.label[!mammaltree2$tip.label %in% names(trnp1.coding.all.fa)])

ggtree(tree.exon1.coding.all)+geom_tiplab()+xlim(0,200)
tree.exon1.coding.all<-multi2di(tree.exon1.coding.all)
is.binary.tree(tree.exon1.coding.all)
length(tree.exon1.coding.all$tip.label)
write.tree(tree.exon1.coding.all,"protein/trees/tree_TRNP1_coding_45sp.txt") 





#prepare data for Coevol####
#shorten names for Coevol & select species with phenotype data

trnp1.coding.all.fa<-readDNAStringSet('protein/fastas/TRNP1_coding_seqs_45sp_full.fa')
pheno_data_with_GI<-readRDS("data_tables/pheno_data/pheno_data.rds") %>%
  filter(!is.na(GI))


trnp1.coding.31sp.fa<-trnp1.coding.all.fa[names(trnp1.coding.all.fa) %in% pheno_data_with_GI$species]
trnp1.coding.31sp.fa.forCoevol<-trnp1.coding.31sp.fa
names(trnp1.coding.31sp.fa.forCoevol)<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",names(trnp1.coding.31sp.fa.forCoevol))
names(trnp1.coding.31sp.fa.forCoevol)<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",names(trnp1.coding.31sp.fa.forCoevol))
names(trnp1.coding.31sp.fa.forCoevol)<-gsub("_", "", names(trnp1.coding.31sp.fa.forCoevol))
names(trnp1.coding.31sp.fa.forCoevol)<-substr(names(trnp1.coding.31sp.fa.forCoevol), 1, 10)

writeXStringSet(trnp1.coding.31sp.fa.forCoevol, filepath='protein/fastas/TRNP1_coding_seqs_31species_forCoevol.fa', format = 'fasta', append=F)



#adjust the tree
tree.exon1.31sp.forCoevol<-drop.tip(tree.exon1.coding.all, tree.exon1.coding.all$tip.label[!tree.exon1.coding.all$tip.label %in% names(trnp1.coding.31sp.fa)])
write.tree(tree.exon1.31sp.forCoevol,"protein/trees/tree_TRNP1_coding_31sp.txt") 

tree.exon1.31sp.forCoevol$tip.label<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops", tree.exon1.31sp.forCoevol$tip.label)
tree.exon1.31sp.forCoevol$tip.label<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",tree.exon1.31sp.forCoevol$tip.label)
tree.exon1.31sp.forCoevol$tip.label<-gsub("_", "", tree.exon1.31sp.forCoevol$tip.label)
tree.exon1.31sp.forCoevol$tip.label<-substr(tree.exon1.31sp.forCoevol$tip.label, 1, 10)
ggtree(tree.exon1.31sp.forCoevol)+geom_tiplab()+xlim(0,200)

write.tree(tree.exon1.31sp.forCoevol,"protein/trees/tree_TRNP1_coding_31sp_forCoevol.txt") 




#adjust the species names in the new pheno data
pheno_31sp_forCoevol<-pheno_data_with_GI %>% 
  filter(species %in% names(trnp1.coding.31sp.fa)) %>%
  dplyr::select(species, body_mass, brain_mass, EQ, GI)

saveRDS(pheno_31sp_forCoevol,"protein/coevol/pheno_data/pheno_data_31species.rds")
 
pheno_31sp_forCoevol$species_short<-pheno_31sp_forCoevol$species
pheno_31sp_forCoevol$species_short<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",pheno_31sp_forCoevol$species_short)
pheno_31sp_forCoevol$species_short<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",pheno_31sp_forCoevol$species_short)
pheno_31sp_forCoevol$species_short<-gsub("_","", pheno_31sp_forCoevol$species_short)
pheno_31sp_forCoevol$species_short<-substr(pheno_31sp_forCoevol$species_short, 1, 10)


insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}

wrap_pheno_coevol<-function(pheno_table, trait_vector){
  pheno_table_selected<-pheno_table[,c("species_short",paste(trait_vector))]
  pheno_table_selected<-insertRow(pheno_table_selected, c("#TRAITS",rep(NA, times=length(trait_vector))),1)
  pheno_table_selected<-insertRow(pheno_table_selected, c(paste0(dim(pheno_table)[1], " ", length(trait_vector)), trait_vector),2)
write.table(pheno_table_selected, paste0("protein/coevol/pheno_data/",paste(trait_vector, collapse = "_"),"_exon1_",dim(pheno_table)[1],"species.lht"), row.names = F, col.names = F, quote = F, na="")
}


wrap_pheno_coevol(pheno_31sp_forCoevol, c("GI"))
wrap_pheno_coevol(pheno_31sp_forCoevol, c("brain_mass"))
wrap_pheno_coevol(pheno_31sp_forCoevol, c("body_mass"))
wrap_pheno_coevol(pheno_31sp_forCoevol, c("body_mass","brain_mass","GI"))





# prepare data for PAML ####
trnp1.coding.all.forPAML<-trnp1.coding.all.fa
names(trnp1.coding.all.forPAML)<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",names(trnp1.coding.all.forPAML) )
names(trnp1.coding.all.forPAML)<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",names(trnp1.coding.all.forPAML) )
names(trnp1.coding.all.forPAML)<-substr(names(trnp1.coding.all.forPAML), 1, 10)

writeXStringSet(trnp1.coding.all.forPAML, filepath='protein/fastas/TRNP1_coding_seqs_45species_forPAML.fa', format = 'fasta', append=F)



#adjust tree tip.labels for PAML
tree.exon1.45sp.forPAML<-tree.exon1.coding.all
tree.exon1.45sp.forPAML$tip.label<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",tree.exon1.45sp.forPAML$tip.label)
tree.exon1.45sp.forPAML$tip.label<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",tree.exon1.45sp.forPAML$tip.label)
tree.exon1.45sp.forPAML$tip.label<-substr(tree.exon1.45sp.forPAML$tip.label, 1, 10)
ggtree(tree.exon1.45sp.forPAML)+geom_tiplab()+xlim(0,200)

write.tree(tree.exon1.45sp.forPAML,"protein/trees/tree_TRNP1_coding_45sp_forPAML.txt") 


