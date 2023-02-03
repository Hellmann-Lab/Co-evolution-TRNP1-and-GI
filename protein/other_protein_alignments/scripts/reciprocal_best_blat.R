#!/usr/bin/env Rscript
 
require(optparse)
 
option_list <- list(
   make_option(c("--ref2bit"), type="character", default=NULL, 
               help="genome 2bit file of", metavar="character"),
   make_option(c("--other2bit"), type="character", default=NULL, 
               help="genome 2bit file of", metavar="character"),
   make_option(c("--refOOC"), type="character", default=NULL, 
               help="genome ooc file of the reference species", metavar="character"),
   make_option(c("--otherOOC"), type="character", default=NULL, 
               help="genome ooc file of the other species", metavar="character"),
   make_option(c("--refprotein"), type="character", default=NULL, 
               help="fasta file with the protein sequences of the reference genome", metavar="character"),
   make_option(c("--refCDSbed"), type="character", default=NULL, 
               help="bed file with coordinates for the refprotein in refgenome", metavar="character"),
   make_option(c("--outbase"), type="character", default=NULL, 
               help="basename for all files", metavar="character"),
   make_option(c("--outdir"), type="character", default=NULL, 
               help="directory to which the output should be written", metavar="character"),
   make_option(c("--prop_prot_length"), type="numeric", default=0.3, 
               help="filter for the min proportion of the human protein length that should be covered", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

library(methods)
library(plyranges)
library(dplyr)
library(Biostrings)
 
 pslxColumns<-c("match","mismatch","rep_match","Ns",
                "Q_gap_count","Q_gap_bases", "T_gap_count", "T_gap_bases", "strand", 
                "Q_name", "Q_size", "Q_start", "Q_end",
                "T_name", "T_size", "T_start", "T_end",
                "block_count","block_sizes","Q_starts","T_starts",
                "Q_seq","T_seq")
 
 
 ###### costfun adapt to protein this is Still for DNA check out if you find a blosummatrix then the sequences
costfun<-function(match,mismatch, gaps, gapN){
   cost<- list( m=1, mm = 1, go=4, ge=1)
   return( match* cost$m - mismatch*cost$mm - gaps*cost$go - log(gapN+1)*cost$ge )
 }
 
 
costfunAA<-function(aa1,aa2, n, bp){
  aa1<-Biostrings::AAString( gsub(",","",aa1) )
  aa2<-Biostrings::AAString( gsub(",","",aa2) )
  sc<- Biostrings::pairwiseAlignment(aa1, aa2, substitutionMatrix = "BLOSUM62", gapOpening = 3, gapExtension = 1) %>% 
       Biostrings::score()
  
  return(sc- n*3 - bp)
}
########## makeFasta From pslx
makeFasta<-function( tab , filename){
   tab<- tab %>% dplyr::mutate( Q_name = paste0(">",Q_name),
                                T_seq= gsub(",","", T_seq) ) 
   write(unlist(t(tab)),file=filename,ncolumns = 1,sep="\n")
}
 
 
# blat sequences
blat_c1<- paste("/opt/bin/blat", opt$other2bit,
                    opt$refprotein,
                   "-out=pslx -q=prot -t=dnax -mask=lower", # -minIdentity=80",
                    #paste0("-ooc=",opt$otherOOC),
                    paste0(opt$outdir,"/",opt$outbase,"_blat_to_other.outpslx") )
   print("blat no 1")
   print(blat_c1)
   system(blat_c1,wait = T)

#I guess this is where it broke the first time!

####################
## filter highest scoring alignment
pslxH<- data.table::fread(paste0(opt$outdir,"/",opt$outbase,
                                  "_blat_to_other.outpslx"),
                             skip = 5, col.names = pslxColumns) %>%
   rowwise() %>%
   #the cost function (BLOSUM matrix) does not like ambiguous bases--> just change these to "any" at this point (X)
   dplyr::mutate(T_seq2=gsub("B|J|Z|O|U","X",T_seq),
                 Q_seq2=gsub("B|J|Z|O|U","X",Q_seq),
                 score = costfunAA( Q_seq2, T_seq2, Q_gap_count, Q_gap_bases)) %>%
   dplyr::group_by( Q_name) %>%
   dplyr::filter( score == max(score), match/Q_size>opt$prop_prot_length) %>% ungroup() %>%
    dplyr::mutate(index_other=1:nrow(.))

write.table(pslxH,file = paste0(opt$outdir,"/",opt$outbase,"_blat_to_other_filt.outpslx"))

#pslxH<-read.table(file = paste0(opt$outdir,"/",opt$outbase,"_blat_to_other_filt.outpslx"))

## write the other species sequences that matched to the reference out as fasta
makeFasta(  pslxH %>% ungroup() %>% mutate(Q_name=paste0(index_other,"_",Q_name)) %>%
              dplyr::select(Q_name,T_seq),
             filename = paste0(opt$outdir,"/",opt$outbase,"_blat_to_other.fa") )


rm(pslxH)

## reciprocal blat
blat_c2<- paste("/opt/bin/blat", opt$ref2bit,
                paste0(opt$outdir,"/",opt$outbase,"_blat_to_other.fa"),
                "-out=pslx  -q=prot -t=dnax -mask=lower",
                #paste0("-ooc=",opt$refOOC),
                paste0(opt$outdir,"/",opt$outbase,"_blat_to_REF.outpsl") )

print("blat no 2")
system(blat_c2, wait = T)
 
########################
# read in reciprocal blat, i.e. the ref. coordinates for best hits
pslR<- data.table::fread(paste0(opt$outdir,"/",opt$outbase,"_blat_to_REF.outpsl"),
                         skip=5, col.names = pslxColumns) %>%
  rowwise() %>%
  dplyr::mutate(T_seq2=gsub("B|J|Z|O|U","X",T_seq),
                Q_seq2=gsub("B|J|Z|O|U","X",Q_seq),
                score = costfunAA( T_seq2, Q_seq2, Q_gap_count, Q_gap_bases)) %>%
  tidyr::separate(col = "Q_name", into = c("index_other", "Q_name"),sep = "_" ) %>%
  dplyr::group_by( Q_name, Q_size) %>%
  dplyr::filter( score == max(score)) %>% 
  ungroup()

# example: "/data/share/htp/TRNP1/other_protein_alignments/blat_output_full/Pongo_abelii/Pongo_abelii_blat_to_REF.outpsl"


#read original cds coordinates of REF
library(readr)
library(tidyverse)

CCDS<-read_tsv("protein/other_protein_alignments/CCDS/CCDS.current.bed", col_names = F) 
colnames(CCDS)<-c("seqnames","start","end","strand","name","gene_name","accession")
CCDS <-CCDS %>%
  filter(!start=="-")
CCDS$start<-as.numeric(CCDS$start)
CCDS$end<-as.numeric(CCDS$end)
refBED<-CCDS %>%
  mutate(index_ref=1:nrow(.), start_ref=start, end_ref=end) %>%
  drop_na() %>%
  as_granges() %>%
  mutate(exon.width.ref = width(.))

# # make a granges object out of the reciprocal blat results
rblat <- pslR %>%
  distinct()%>%
  mutate(seqnames=T_name,
         start=T_start,
         end=T_end,
         strand=substr(strand,2,2),
         name = word(Q_name,1,1,"[|]"),
         index_rblat=1:nrow(.),
         start_rblat=T_start,
         end_rblat=T_end) %>%
  dplyr::select(seqnames, start, end, strand, name,index_rblat, start_rblat, end_rblat, Q_start, Q_end, Q_seq2, T_seq2, score, index_other) %>%
  as_granges() %>%
  mutate(exon.width.rblat=width(.))


#find and quantify overlaps 
hits <- findOverlaps(refBED, rblat)
overlaps <- pintersect(refBED[queryHits(hits)], rblat[subjectHits(hits)])
hits_overlaps<-data.frame(hits, width(overlaps))

refGTFreci<-join_overlap_intersect(refBED, rblat, suffix = c(".ref", ".rblat")) %>%
  as_tibble() %>%
  inner_join(hits_overlaps, by=c("index_ref"="queryHits", "index_rblat"="subjectHits")) %>%
  dplyr::rename(overlap_width=`width.overlaps.`) %>%
  mutate(percOverlap_ref=overlap_width/exon.width.ref) %>%
  filter(name.ref==name.rblat, percOverlap_ref>0.1) %>% 
  distinct() %>%
  #directly filter for the max overlap per ref sequence
  group_by(name.ref) %>%
  slice_max(order_by=percOverlap_ref, n=1) %>%
  ungroup()

saveRDS(refGTFreci,paste0(opt$outdir,"/",opt$outbase,"_blat_overlap.rds"))




