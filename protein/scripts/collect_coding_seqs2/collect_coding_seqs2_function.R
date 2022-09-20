#!/opt/R/4.0.0/lib/R/bin/exec/R

#load libraries
libs <- c( "Biostrings","tidyverse", "rBLAST", "stringr")
sapply(libs,require,character.only=T)

setwd("/data/share/htp/TRNP1/paper_data/")

#read arguments
args <- commandArgs()
print(args)

input_species <-args[6]
print(input_species)

input_seq_dir <-args[7]
print(input_seq_dir)

print("first readin")
print(paste0(input_seq_dir, "/", args[8]))

input_seq <-Biostrings::readDNAStringSet(paste0(input_seq_dir, "/", args[8]))
head(input_seq)

output_dir <- args[9]
print(output_dir)

print("before function")


#load human TRNP1 protein sequence
protein.fa <- readAAStringSet("protein/data/trnp1_protein_seq_hum.fasta")


#create a blast function
do_rblast_full<-function(subject=input_seq,query=protein.fa,type="tblastn", min_perc_identity=70, min_alignment_length=150, alt_min_perc_identity=90, species=input_species){
  print("within function")
  
  # the genomes from NCBI have a very weird long gene names, that are being cut off after the space (which is fine because the first word is the seq id)
  if (grepl("*.fna*",args[8])==TRUE){
  #generate blast db of the genome -- omg, need to unzip for this..
    system(paste0("gunzip --keep ", input_seq_dir, "/", args[8]))
    db<-makeblastdb(gsub(".gz","",paste0(input_seq_dir, "/", args[8])), dbtype = "nucl")
    res<-blast(gsub(".gz","",paste0(input_seq_dir, "/", args[8])),type="tblastn")
    #switch off masking (normally masks overrepresented or low-complexity sequences)
    blast_output<-stats::predict(res,query, BLAST_args="-soft_masking false -seg no")
    system(paste("rm", paste0(input_seq_dir, "/", "{*.nsq*,*.nin*,*.nhr*,*.fna}")))
  
  } 
  else {
    #whereas in other cases, the full name (which I concatenate together) is required to uniquely identify the sequences
    dir<-paste(input_seq_dir)
    names(subject)<-gsub(" ","",names(subject))
    writeXStringSet(subject, filepath = file.path(paste0(dir,"/",species,".Newseqs.fasta")))
    db<-makeblastdb(file.path(paste0(dir,"/",species,".Newseqs.fasta")), dbtype = "nucl")
    res<-blast(file.path(paste0(dir,"/",species,".Newseqs.fasta")),type=paste0(type))
    #switch off masking (normally masks overrepresented or low-complexity sequences)
    blast_output<-stats::predict(res,query, BLAST_args="-soft_masking false -seg no")
    system(paste0("rm ", dir,"/",species,".Newseqs.fasta*"))
  }
  
  #filter for high identity and alignment length
  blast_output_filt<-blast_output %>% filter((Perc.Ident>min_perc_identity & Alignment.Length>min_alignment_length) | Perc.Ident>alt_min_perc_identity)
  
  #select the relevant sequence
  #subject_seq<-subject[names(subject)==blast_output_filt$SubjectID[1]]
  subject_seq<-subject[grepl(paste0(blast_output_filt$SubjectID[1]," |",blast_output_filt$SubjectID[1],"$"),names(subject))]
  
  #cut the region and do revCompl where necessary
  if (blast_output_filt$S.end[1]<blast_output_filt$S.start[1]){
    subject_cut<-subseq(subject_seq, start=blast_output_filt$S.end[1]-24, end=blast_output_filt$S.start[1]+24)
    subject_cut<-Biostrings::reverseComplement(subject_cut)
  } else {
    subject_cut<-subseq(subject_seq, start=blast_output_filt$S.start[1]-24, end=blast_output_filt$S.end[1]+24)
  }
  
  #make an exception for cases where there is a assembly-weirdness that leads to frameshift within the trnp1 coding sequence, such as in the case of M.fascicularis --> add the 2nd sequence part
  if (dim(blast_output_filt)[1]==2 & blast_output_filt$SubjectID[1]==blast_output_filt$SubjectID[2]){
    if (blast_output_filt$S.end[2]<blast_output_filt$S.start[2]){
      subject_cut_2<-subseq(subject_seq, start=blast_output_filt$S.end[2], end=blast_output_filt$S.start[2])
      subject_cut_2<-Biostrings::reverseComplement(subject_cut_2)
    } else {
      subject_cut_2<-subseq(subject_seq, start=blast_output_filt$S.start[2], end=blast_output_filt$S.end[2])
    }
    #concatenate the 2 seqs
    if (blast_output_filt$Q.end[1]<blast_output_filt$Q.start[2]){
      subject_cut<-DNAStringSet(paste0(subject_cut, subject_cut_2))
    }
  }
  
  names(subject_cut)<-paste0(species)
  writeXStringSet(subject_cut, filepath=file.path(output_dir,paste0(species,"_TRNP1_coding.fa")), format = 'fasta', append=F)
}



do_rblast_full()

