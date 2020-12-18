
libs<-c("rBLAST","tidyverse","GenomicRanges","data.table","broom", "ape", "ggtree", "readr","geiger","nlme","phytools","grid","gtable","xtable","phangorn","Biostrings")
sapply(libs, require, character.only=T)



##########################################################
# separate sequences from library according to region ####
##########################################################

minlength.seqs.set<-readDNAStringSet("regulation/data/MPRA/input/minlength.library.seqs.set.fa", format = 'fasta')
mammaltree<-read.tree("protein/trees/mammaltree.txt") 

#adjust names to the tree
names(minlength.seqs.set)<-gsub("Mustela_putorius_furo","Mustela_putorius",names(minlength.seqs.set))
names(minlength.seqs.set)<-gsub("Equus_ferus_caballus","Equus_ferus",names(minlength.seqs.set))
names(minlength.seqs.set)<-gsub("Canis_lupus_familiaris","Canis_lupus",names(minlength.seqs.set))

#remove the double Saguinus oedipus sequence that was accidentally included (keep the filtered ones)
names(minlength.seqs.set[grepl("Saguinus_oedipus", names(minlength.seqs.set))])
minlength.seqs.set<-minlength.seqs.set[!grepl("Saguinus_oedipus_trinity",names(minlength.seqs.set))]
names(minlength.seqs.set)<-gsub("_prom_","",names(minlength.seqs.set))


#generate output directories
system('mkdir -p regulation/fastas/mafft_regions/ regulation/fastas/separated_regions_correct/ regulation/fastas/separated_regions/  regulation/trees/')


#subset sequences from each region, save; generate trees with the respective species.
reg_regions<-c("intron", "upstream3", "upstream2", "upstream1", "exon2", "downstream", "exon1")

for (i in reg_regions){
  seqs.fa<-minlength.seqs.set[grepl(pattern=paste0("*",i,"$"), names(minlength.seqs.set))]
  seqs.fa<-seqs.fa[!grepl(pattern="*node*", names(seqs.fa))]
  names(seqs.fa)<-gsub("(exon1|upstream1|intron|upstream2|upstream3|unique1|unique2|unique3|exon2|downstream|trinity|dhs|peaks|rbh|filt)", "", names(seqs.fa))
  names(seqs.fa)<-gsub("(^_*|*__$|*_$|*___$)", "", names(seqs.fa))
  writeXStringSet(seqs.fa, filepath=paste0("regulation/fastas/separated_regions/TRNP1_",i,"_seqs.fa"), format = 'fasta', append=F)
  
  tree<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% names(seqs.fa)])
  write.tree(tree,paste0("regulation/trees/tree_TRNP1_",i,".txt"))
  
}


###########################
# alignment, reversing ####
###########################

#there is reverse input in some cases- need to turn it around, which mafft is very good at
system('sbatch regulation/scripts/MPRA/mafft_adjust_direction.sh', wait = T) 


###################################
#collect and restore sequences ####
###################################

#first, get the correct sequences in the right direction 
#one should stick to the direction that human has.
regregs.fa<-list()
comb.regregs.fa<-list()

for (i in reg_regions){
  
  mafft_adj<-readDNAStringSet(paste0("regulation/fastas/mafft_regions/mafft_",i,"_adj.fa"))
  
  #so to have the same strand for all regions, we orientate based on the human seq (which should be the + strand in all cases) --> basically, if human was turned around, we do reverse compl on everything to turn it back
  mafft_adj_hum<-mafft_adj[grepl("Homo_sapiens",names(mafft_adj))]
  if (grepl("_R_",names(mafft_adj_hum))==TRUE){
    
    mafft_adj2<-Biostrings::reverseComplement(mafft_adj)
  } else {
    mafft_adj2<-mafft_adj
  }
  
  mafft_adj_noGaps<-gsub("-|N|n","", mafft_adj2)
  mafft_adj_noGaps<-as.data.frame(mafft_adj_noGaps) 
  mafft_adj_noGaps$species<-gsub("_R_","", names(mafft_adj2))
  
  mafft_adj_noGaps.fa<-DNAStringSet(mafft_adj_noGaps$mafft_adj_noGaps)
  names(mafft_adj_noGaps.fa)<-mafft_adj_noGaps$species
  
  writeXStringSet(mafft_adj_noGaps.fa, paste0("regulation/fastas/separated_regions_correct/",i,"_correct_full_seqs.fa"))
  
  comb.regregs.fa[[i]]<-mafft_adj_noGaps.fa
  #add region to the species name
  names(comb.regregs.fa[[i]])<-paste0(names(comb.regregs.fa[[i]]),"_",  names(comb.regregs.fa[i]))
}

all.regregs.correct.fa<-DNAStringSet(do.call(c, unname(unlist(comb.regregs.fa))))
writeXStringSet(all.regregs.correct.fa, paste0("regulation/fastas/separated_regions_correct/all_correct_full_seqs.fa"))
