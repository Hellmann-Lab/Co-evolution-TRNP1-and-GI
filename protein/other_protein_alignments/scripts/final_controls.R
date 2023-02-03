library(Biostrings)
library(tidyverse)
library(plyranges)

discrepancies_list <- list.files("protein/other_protein_alignments/self_matches_diff", recursive = T, full.names = T)[-25] #exclude procavia

discrepancies_rds<-lapply(discrepancies_list, function(x){readRDS(x)})
discrepancies_df<-bind_rows(discrepancies_rds)

unique(discrepancies_df$protein) #ok lets just exclude these


CCDS_322<-readAAStringSet("protein/other_protein_alignments/CCDS/CCDS_protein.current.322.fa")

toKeep<-data.frame(protein=names(CCDS_322)) %>%
  filter(!protein %in% unique(discrepancies_df$protein)) #274




#bind the lists from self matches fa
prot_seqs_sp <- list.files("protein/other_protein_alignments/self_matches_fa", recursive = T, full.names = T)[-25]

prot_seqs_sp_rds<-lapply(prot_seqs_sp, function(x){readRDS(x)})
names(prot_seqs_sp_rds)<-word(prot_seqs_sp,9,9,"/") %>% gsub(".rds","",.)

#select names in toKeep
prot_seqs_sp_rds_filt<-lapply(prot_seqs_sp_rds, function(x){x<-x[names(x) %in% toKeep$protein]})

#keep the first sequence since they are supposed to be identical
prot_seqs_sp_rds_filt2<-prot_seqs_sp_rds_filt
for (i in names(prot_seqs_sp_rds_filt2)){
  prot_seqs_sp_rds_filt2[[i]]<-lapply(prot_seqs_sp_rds_filt2[[i]], function(x){x<-x[1]})
}



#adjust list names (concatenate species and protein name)
for (i in names(prot_seqs_sp_rds_filt2)) {
  names(prot_seqs_sp_rds_filt2[[i]]) <- paste(i, names(prot_seqs_sp_rds_filt2[[i]]), sep="_")
}
dna_list <- prot_seqs_sp_rds_filt2
names(dna_list) <- NULL # otherwise do.call() will return the input list
all.seqs.fas <- do.call(c, dna_list)

#set list names to be the sequence names
for (i in 1:length(all.seqs.fas)) {
  names(all.seqs.fas[[i]]) <- names(all.seqs.fas)[[i]]
}

#generate a final DNAString set
all.seqs.fas.final <- DNAStringSet(do.call(c, unname(unlist(all.seqs.fas))))


#add the human CCDS nucleotide sequences (download from UCSC)
humCCDS.fa<-readDNAStringSet("protein/other_protein_alignments/CCDS/CCDS_nucleotide.current.fna.gz")
humCCDS.fa.ctrls<-humCCDS.fa[names(humCCDS.fa) %in% toKeep$protein]
#remove the stop codons because they are missing in the orthologous seqs because of prot-nuc blat
humCCDS.fa.ctrls<-subseq(humCCDS.fa.ctrls, start=1, end=width(humCCDS.fa.ctrls)-3)
names(humCCDS.fa.ctrls)<-paste("Homo_sapiens",names(humCCDS.fa.ctrls), sep="_")

all.seqs.fas.final<-c(all.seqs.fas.final,humCCDS.fa.ctrls)



#make a more stringent version regarding sequence lengths (must be longer than 0.5 that of human and shorter than 2x that of human) -- only lost 9 additional proteins, so sticking with this one


for (prot in toKeep$protein){
  
  #filter for the protein
  prot.fas.final<-all.seqs.fas.final[word(names(all.seqs.fas.final),3,3,"_")==prot]
  
  #change names to species names and shorten for Coevol
  names(prot.fas.final)<-word(names(prot.fas.final),1,2,"_")
  names(prot.fas.final)<-gsub("Chlorocebus_aethiops", "Chloroc_aethiops",names(prot.fas.final))
  names(prot.fas.final)<-gsub("Chlorocebus_sabeus", "Chloroc_sabeus",names(prot.fas.final))
  names(prot.fas.final)<- gsub("Cercopithecus_mona","Cercopithecus_mitis", names(prot.fas.final))
  names(prot.fas.final)<-gsub("_", "", names(prot.fas.final))
  names(prot.fas.final)<-substr(names(prot.fas.final), 1, 10)
  
  #save human seq length
  hum.width<-width(prot.fas.final[names(prot.fas.final)=="Homosapien"])
  
  if ((min(width(prot.fas.final))>0.5*hum.width | max(width(prot.fas.final))<2*hum.width) & length(prot.fas.final)==30) {

    writeXStringSet(prot.fas.final, filepath =paste("protein/other_protein_alignments/blat_322cand/prank_input/","",word(prot,1,1,"[|]"), ".fa", sep=""), format = "fasta", append = F)
  } else {
    print(prot)
  }
}
# final number: 257

