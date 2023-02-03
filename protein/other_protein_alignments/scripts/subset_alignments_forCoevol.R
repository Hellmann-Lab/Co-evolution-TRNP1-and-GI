
#check out and subset the good alignments

aln_QC<-readr::read_tsv("protein/other_protein_alignments/prank_out/Coevol_proteins_aln1 - protein_list_all.tsv")


good<-aln_QC %>% filter(decision=="use") %>% filter(!notes %in% c("short","sus_gap"))
write_tsv(good[,"protein"],"protein/other_protein_alignments/prank_out/good_noCuration.tsv", col_names = F)


#ones that might be worth curating
curate<-aln_QC %>% filter(decision=="curate" | notes=="sus_gap") %>% filter(!grepl("short",notes))
write_tsv(curate,"protein/other_protein_alignments/better_genomes/needCuration.tsv")






#now I have looked at the curated sequences ####

#(workflow: folder better_genomes) and made a decision regarding them
#output in better_genomes/curated folder

aln_QC_modified<-readr::read_tsv("protein/other_protein_alignments/better_genomes/prank_out/Coevol_proteins_aln1_modified - protein_list_all.tsv")

#these are good as they are
aln_QC_modified_use<-aln_QC_modified %>% dplyr::filter(decision3=="use")

#copy these into prank_input (this seems the easiest way currently, also consistent with the following)
for (i in aln_QC_modified_use$protein){
    seq<-readDNAStringSet(paste0("protein/other_protein_alignments/better_genomes/prank_input/",i))
    writeXStringSet(seq,paste0("protein/other_protein_alignments/better_genomes/curated/prank_input/",i))
}



#these ones need curation 
aln_QC_modified_curate<-aln_QC_modified %>%
  dplyr::filter(decision3=="curate") %>%
  separate(to_cut, into=c("range1","range2"),sep="; ") %>%
  separate(range1, into=c("start1","end1"),sep="-") %>%
  separate(range2, into=c("start2","end2"),sep="-") %>%
  #mutate(start1=as.numeric(start1), start2=as.numeric(start2), 
         #end1=as.numeric(end1), end2=as.numeric(end2)) %>%
  #need to remember - coordinates are given in AA space, but here we're working with nucleotides
  rowwise() %>%
  mutate(start1=ifelse(start1==1,as.numeric(start1),as.numeric(start1)*3),
         end1=as.numeric(end1)*3,
         start2=as.numeric(start2)*3,
         end2=as.numeric(end2)*3) %>%
  #adjust names
  mutate(species2=stringr::str_to_title(species2),
         species2=gsub("_","",species2),
         species2=substr(species2,1,10))


#there are different types of action to take:

# 1) for these I shall cut the species in species2 (better_genomes/prank_input for decision2 curate, blat_322cand/prank_input for decision2 take_1) --> realign
c1<-aln_QC_modified_curate %>% filter(decision2 %in% c("curate","take_1"), type_seq=="species") 

for (i in c1$protein){
  
  if(c1$decision2[c1$protein==i]=="curate"){
  seq<-readDNAStringSet(paste0("protein/other_protein_alignments/better_genomes/prank_input/",i))
  
  } else {
    
    if (c1$decision2[c1$protein==i]=="take_1"){
      seq<-readDNAStringSet(paste0("protein/other_protein_alignments/prank_input/",i))
    }
  }
  
  sp_to_curate<-c1[c1$protein==i,]
  seq_to_curate<-seq[grepl(sp_to_curate$species2,names(seq))]
  
  seq_curated<-subseq(seq_to_curate, start=sp_to_curate$start1[1], end=sp_to_curate$end1[1])
  #if there's a second region
  if (is.na(sp_to_curate$start2[1])==FALSE){
    seq_curated2<-subseq(seq_to_curate, start=sp_to_curate$start2[1], end=sp_to_curate$end2[1])
    #concatenate
    seq_curated<-DNAStringSet(paste0(seq_curated, seq_curated2))
  }
  
  seq[grepl(sp_to_curate$species2,names(seq))]<-seq_curated
  writeXStringSet(seq,paste0("protein/other_protein_alignments/better_genomes/curated/prank_input/",i))
}




#2) these I shall cut on the level of alignment (better_genomes/prank_out) or blat_322cand/prank_out
c2<-aln_QC_modified_curate %>% filter(decision2 %in% c("curate","take_1") & type_seq=="aln")

for (i in c2$protein){
  
  if(c2$decision2[c2$protein==i]=="curate"){
    
    seq<-readDNAStringSet(paste0("protein/other_protein_alignments/better_genomes/prank_out/aln_forChecking/",i, ".best.nuc.fas"))
    
  } else {
    
    if (c2$decision2[c2$protein==i]=="take_1"){
  
  seq<-readDNAStringSet(paste0("protein/other_protein_alignments/prank_out/aln_forChecking/",i, ".best.nuc.fas"))
    }
  }
  
  sp_to_curate<-c2[c2$protein==i,]

  seq_curated<-subseq(seq, start=sp_to_curate$start1[1], end=sp_to_curate$end1[1])
  #if there's a second region
  if (is.na(sp_to_curate$start2[1])==FALSE){
    seq_curated2<-subseq(seq, start=sp_to_curate$start2[1], end=sp_to_curate$end2[1])
    #concatenate
    seq_curated<-DNAStringSet(paste0(seq_curated, seq_curated2))
  }

    seq_curated_noGaps<-DNAStringSet(gsub("-|N|n","",seq_curated))
    names(seq_curated_noGaps)<-names(seq)
    #save unaligned and realign
    writeXStringSet(seq_curated_noGaps,paste0("protein/other_protein_alignments/better_genomes/curated/prank_input/",i))
}





#Additional very complicated case where I have to use the old gorilla and mandrill seqs and cut the current Cercocebus seq --> realign
c3<-aln_QC_modified_curate %>% filter(decision2=="take_oldGorilla,Mandrill")

seq_add<-readDNAStringSet(paste0("protein/other_protein_alignments/better_genomes/prank_input/",c3 %>% pull(protein)))

#first, cut Cercocebus
seq_add[names(seq_add)=="Cercocebus"]<-subseq(seq_add[names(seq_add)=="Cercocebus"], start=c3$start1[1], end=c3$end1[1])

#now pull out the previous Gorilla and Mandrill sequences and exchange the current with the old ones
seq_add_old<-readDNAStringSet(paste0("protein/other_protein_alignments/blat_322cand/prank_input/",c3 %>% pull(protein)))

seq_add[names(seq_add)=="Gorillagor"]<-seq_add_old[names(seq_add)=="Gorillagor"]
seq_add[names(seq_add)=="Mandrillus"]<-seq_add_old[names(seq_add)=="Mandrillus"]

writeXStringSet(seq_add,paste0("protein/other_protein_alignments/better_genomes/curated/prank_input/",c3 %>% pull(protein)))



