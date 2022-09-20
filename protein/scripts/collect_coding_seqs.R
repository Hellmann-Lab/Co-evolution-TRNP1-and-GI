# This script contains TRNP1 protein-coding orthologue search by Zane Kliesmete

setwd("/data/share/htp/TRNP1/paper_data/")

libs <- c("IRanges", "Biostrings", "tidyverse", "reshape2", "reshape", "GenomicRanges", "data.table", "phangorn", "ape", "ggrepel", "ggtree", "caper", "phytools", "rBLAST", "stringr","xtable", "geiger","nlme","phytools")

sapply(libs, require, character.only=T)


# load resequenced promoter-exon1 region from 18 species ####
# data generation in protein/scripts/Rscripts/processingL_final.R by Lucas Wange
reseq.prom_exon1.fas<-readDNAStringSet('protein/fastas/resequenced_promexon1_primates.fa')


#load human TRNP1 protein sequence####
protein.fa <- readAAStringSet("protein/fastas/trnp1_protein_seq_hum.fasta")


#create a blast function
do_rblast<-function(subject,query=protein.fa,type="tblastn", min_perc_identity=70, min_alignment_length=190, extend=TRUE){
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
      if (extend==TRUE){
     subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.end[i]-3, end=blast_output_filt$S.start[i])
      } else {
        subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.end[i], end=blast_output_filt$S.start[i])
      }
     subject_cut[i]<-Biostrings::reverseComplement(subject_cut[i])
     
   } else {
     if (extend==TRUE){
       subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.start[i], end=blast_output_filt$S.end[i]+3)
     } else {
       subject_cut[i]<-subseq(subject_seq[i], start=blast_output_filt$S.start[i], end=blast_output_filt$S.end[i])
     }}}
  
  return(list(subject_cut,blast_output_filt)) 
}


##################
# resequenced ####
##################

names(reseq.prom_exon1.fas)<-word(names(reseq.prom_exon1.fas),1,2,"_")
blast.reseq.promexon1.out.all<-do_rblast(reseq.prom_exon1.fas, extend=FALSE) 
reseq.coding.fas<-blast.reseq.promexon1.out.all[[1]]
# writeXStringSet(reseq.coding.fas,"protein/fastas/resequenced_promexon1_6full.fa")
blast.reseq.promexon1.out<-blast.reseq.promexon1.out.all[[2]]



##############################################
# ferret sequence from later resequencing ####
##############################################

mputorius.coding.fa<-readDNAStringSet("protein/fastas/TRNP1_coding/Mustela_putorius_reseq.fa")
names(mputorius.coding.fa)<-"Mustela_putorius"

blast.ferret.out<-do_rblast(mputorius.coding.fa, min_perc_identity = 50, min_alignment_length = 150, extend=FALSE)
ferret.coding.fas<- subseq(mputorius.coding.fa, start=blast.ferret.out[[2]]$S.start, end=blast.ferret.out[[2]]$S.end )



###########
# NCBI ####
###########

# additional species from NCBI --> already pre-selected regions using the online blast tool: for refinement, blast again
NCBI_Delphinapterus_leucas_TRNP1.fas<-readDNAStringSet("protein/fastas/TRNP1_coding/NCBI_Delphinapterus_leucas_TRNP1.fa")
names(NCBI_Delphinapterus_leucas_TRNP1.fas)<-"Delphinapterus_leucas"
 
NCBI_Odocoileus_virginianus_TRNP1.fas<-readDNAStringSet("protein/fastas/TRNP1_coding/NCBI_Odocoileus_virginianus_TRNP1.fa")
names(NCBI_Odocoileus_virginianus_TRNP1.fas)<-"Odocoileus_virginianus"
 
NCBI_species_TRNP1.fas<-c(NCBI_Delphinapterus_leucas_TRNP1.fas, NCBI_Odocoileus_virginianus_TRNP1.fas)
#blast
blast.NCBI.out.all<-do_rblast(NCBI_species_TRNP1.fas, extend=FALSE)
NCBI.coding.fas<-blast.NCBI.out.all[[1]]
blast.NCBI.out<-blast.NCBI.out.all[[2]]  #looking great



###################################################################################
# more seqs from UCSC, NCBI, ENSEMBL genomes (these have been blasted already) ####
###################################################################################

more.coding.files<-list.files("protein/fastas/TRNP1_coding/coding_seqs_sep2_extended/", full.names =T, pattern = ".fa")
more.coding.fas.list <- lapply(more.coding.files, function(x) readDNAStringSet(x))
more.coding.fas<-DNAStringSet(do.call(c, unname(unlist(more.coding.fas.list))))
#remove the 24 bp extension for Dasypus_novemcinctus and Papio_hamadryas which have Ns and nonsense there
to_cut<-more.coding.fas[names(more.coding.fas) %in% c("Dasypus_novemcinctus", "Papio_hamadryas")]
for (i in 1:length(to_cut)){
  to_cut[i]<-subseq(to_cut[i], start=1, end=width(to_cut[i])-24)
  more.coding.fas[names(more.coding.fas)==names(to_cut[i])]<-to_cut[i]
}




#throw in the actual human coding sequence to know where the prot starts and ends
blast.reseq.promexon1.hum<-do_rblast(reseq.prom_exon1.fas[names(reseq.prom_exon1.fas)=="Homo_sapiens"], extend=TRUE) 
reseq.coding.hum.fas<-blast.reseq.promexon1.hum[[1]]
# writeXStringSet(reseq.coding.fas,"protein/fastas/resequenced_promexon1_6full.fa")
blast.reseq.promexon1.hum.out<-blast.reseq.promexon1.hum[[2]]

more.coding.fas.withHum<-c(reseq.coding.hum.fas, more.coding.fas)
writeXStringSet(more.coding.fas.withHum, filepath='protein/fastas/TRNP1_coding_seqs_39sp_extended.fa', format = 'fasta', append=F)


#ALIGN WITH MAFFT
system("/home/zane/TRNP1/mafft7/mafft-linux64/mafft.bat /data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_39sp_extended.fa > /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/mafft_TRNP1_coding_39species_extended.fa")


#cut out the sequences from the mafft alignment, align with prank afterwards
aln39<-readDNAMultipleAlignment("protein/fastas/prank_output/mafft_TRNP1_coding_39species_extended.fa")
aln39_trimmed <- DNAMultipleAlignment(aln39,start=27,end=746) #without stop codon
aln39_trimmed

aln39_trimmed2<-gsub("-","", aln39_trimmed)
aln39_trimmed2_noGaps<-as.data.frame(aln39_trimmed2) 

aln39_trimmed2_noGaps.fa<-DNAStringSet(aln39_trimmed2_noGaps$aln39_trimmed2)
names(aln39_trimmed2_noGaps.fa)<-names(DNAStringSet(aln39))


#combine with other sequences ####
aln39_trimmed2_noGaps.fa<-aln39_trimmed2_noGaps.fa[!names(aln39_trimmed2_noGaps.fa) %in% c(names(NCBI.coding.fas), names(reseq.coding.fas))]


#################
#combine all ####
#################

trnp1.coding.all.fa<-c(reseq.coding.fas,
                        NCBI.coding.fas,
                        ferret.coding.fas,
                        aln39_trimmed2_noGaps.fa)


writeXStringSet(trnp1.coding.all.fa, filepath='protein/fastas/TRNP1_coding_seqs_45sp_longer.fa', format = 'fasta', append=F)






###################################
# adjust the phylogenetic tree ####
###################################

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





##############################
# prepare data for Coevol ####
##############################
# shorten names for Coevol (max10 characters) & select species with phenotype data


trnp1.coding.all.fa2<-readDNAStringSet('protein/fastas/TRNP1_coding_seqs_45sp_longer.fa')
pheno_data_with_GI<-readRDS("Co-evolution-TRNP1-and-GI/pheno_data/pheno_data.rds") %>%
  filter(!is.na(GI), species!="Procavia_capensis")

trnp1.coding.30sp.fa<-trnp1.coding.all.fa2[names(trnp1.coding.all.fa2) %in% pheno_data_with_GI$species]
trnp1.coding.30sp.fa.forCoevol<-trnp1.coding.30sp.fa
names(trnp1.coding.30sp.fa.forCoevol)<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",names(trnp1.coding.30sp.fa.forCoevol))
names(trnp1.coding.30sp.fa.forCoevol)<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",names(trnp1.coding.30sp.fa.forCoevol))
names(trnp1.coding.30sp.fa.forCoevol)<-gsub("_", "", names(trnp1.coding.30sp.fa.forCoevol))
names(trnp1.coding.30sp.fa.forCoevol)<-substr(names(trnp1.coding.30sp.fa.forCoevol), 1, 10)

writeXStringSet(trnp1.coding.30sp.fa.forCoevol, filepath='protein/fastas/TRNP1_coding_seqs_30species_forCoevol_longer.fa', format = 'fasta', append=F)


#adjust the tree
tree.exon1.30sp.forCoevol<-drop.tip(mammaltree2, mammaltree2$tip.label[!mammaltree2$tip.label %in% names(trnp1.coding.30sp.fa)])
tree.exon1.30sp.forCoevol<-multi2di(tree.exon1.30sp.forCoevol)
write.tree(tree.exon1.30sp.forCoevol,"protein/trees/tree_TRNP1_coding_30sp.txt") 

tree.exon1.30sp.forCoevol$tip.label<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops", tree.exon1.30sp.forCoevol$tip.label)
tree.exon1.30sp.forCoevol$tip.label<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",tree.exon1.30sp.forCoevol$tip.label)
tree.exon1.30sp.forCoevol$tip.label<-gsub("_", "", tree.exon1.30sp.forCoevol$tip.label)
tree.exon1.30sp.forCoevol$tip.label<-substr(tree.exon1.30sp.forCoevol$tip.label, 1, 10)
ggtree(tree.exon1.30sp.forCoevol)+geom_tiplab()+xlim(0,200)

write.tree(tree.exon1.30sp.forCoevol,"protein/trees/tree_TRNP1_coding_30sp_forCoevol.txt") 



#adjust the species names in the new pheno data
pheno_30sp_forCoevol<-pheno_data_with_GI %>% 
  filter(species %in% names(trnp1.coding.30sp.fa)) %>%
  dplyr::select(species, body_mass, brain_mass, EQ, GI)
saveRDS(pheno_30sp_forCoevol,"protein/coevol/pheno_data/pheno_data_30species.rds")

pheno_30sp_forCoevol$species_short<-pheno_30sp_forCoevol$species
pheno_30sp_forCoevol$species_short<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",pheno_30sp_forCoevol$species_short)
pheno_30sp_forCoevol$species_short<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",pheno_30sp_forCoevol$species_short)
pheno_30sp_forCoevol$species_short<-gsub("_","", pheno_30sp_forCoevol$species_short)
pheno_30sp_forCoevol$species_short<-substr(pheno_30sp_forCoevol$species_short, 1, 10)


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

#30 species
wrap_pheno_coevol(pheno_30sp_forCoevol, c("GI"))
wrap_pheno_coevol(pheno_30sp_forCoevol, c("brain_mass"))
wrap_pheno_coevol(pheno_30sp_forCoevol, c("body_mass"))
wrap_pheno_coevol(pheno_30sp_forCoevol, c("body_mass","brain_mass","GI"))



############################
# prepare data for PAML ####
############################

trnp1.coding.all.forPAML<-trnp1.coding.all.fa
names(trnp1.coding.all.forPAML)<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",names(trnp1.coding.all.forPAML) )
names(trnp1.coding.all.forPAML)<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",names(trnp1.coding.all.forPAML) )
names(trnp1.coding.all.forPAML)<-substr(names(trnp1.coding.all.forPAML), 1, 10)
 
writeXStringSet(trnp1.coding.all.forPAML, filepath='protein/fastas/TRNP1_coding_seqs_45species_forPAML_longer.fa', format = 'fasta', append=F)


#adjust tree tip.labels for PAML
tree.exon1.45sp.forPAML<-tree.exon1.coding.all
tree.exon1.45sp.forPAML$tip.label<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",tree.exon1.45sp.forPAML$tip.label)
tree.exon1.45sp.forPAML$tip.label<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",tree.exon1.45sp.forPAML$tip.label)
tree.exon1.45sp.forPAML$tip.label<-substr(tree.exon1.45sp.forPAML$tip.label, 1, 10)
ggtree(tree.exon1.45sp.forPAML)+geom_tiplab()+xlim(0,200)
 
write.tree(tree.exon1.45sp.forPAML,"protein/trees/tree_TRNP1_coding_45sp_forPAML.txt") 
 
 

#make a species-name decoder for the short versions
trnp1.coding.all.fa2<-readDNAStringSet('protein/fastas/TRNP1_coding_seqs_45sp_longer.fa')
species_names<-data.frame(species=names(trnp1.coding.all.fa2)) %>%
  mutate(speciesCoevol=gsub("Chlorocebus_aethiops", "Chloroc_aethiops",species),
         speciesCoevol=gsub("Chlorocebus_sabeus", "Chloroc_sabeus",speciesCoevol),
         speciesPAML=substr(speciesCoevol,1,10),
         speciesCoevol=gsub("_", "", speciesCoevol),
         speciesCoevol=substr(speciesCoevol, 1, 10))

saveRDS(species_names, "data_tables/species_names.rds")

