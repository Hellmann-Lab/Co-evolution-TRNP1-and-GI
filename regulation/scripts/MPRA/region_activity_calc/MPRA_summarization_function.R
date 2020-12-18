#write a new summarization function

#1) blast back to the respective region in the respective species to know the exact positions covered by each chunk
#2) make a decoder that spits out per base coverage, aka for each base which chunks that are present in the library are covering and what was their per base activity
#3) sum the per bp mean activity


#mark how long the function runs
old <- Sys.time() # get start time


#load libraries
libs <- c("IRanges", "Biostrings", "tidyverse", "reshape2", "reshape", "GenomicRanges", "data.table", "phangorn", "caper", "phytools", "rBLAST", "stringr")
sapply(libs, require, character.only=T)

#read arguments
args <- commandArgs()
print(args)

region_to_get <-args[6]
print(region_to_get)



#1) blast chunks to the regions####

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


#get the included chunks
median_activity_MPRA<-readRDS("regulation/data/MPRA/output/median_activity_MPRA.rds") %>%
  left_join(mpra_metadata)


#exclude the regions consisting of only one tile since many of them seem hard to blast and they are not overlapping with anything (will add them back at the end)
low_length_regs<-mpra_metadata %>%
  group_by(species, region) %>%
  summarise(n=length(enhancer), SeqID=paste0(SeqID, collapse=",")) %>%
  dplyr::filter(n==1) %>%
  dplyr::filter(region==paste0(region_to_get)) %>%
  inner_join(median_activity_MPRA) 


region_chunk_seqs<-median_activity_MPRA %>%
  #select just one cell line (need the chunk-sequences just once)
  dplyr::filter(region==paste0(region_to_get) & cell_line=="human1") %>%
  dplyr::select(SeqID, species, sequence) %>%
  dplyr::filter(!(SeqID %in% low_length_regs$SeqID))

region_chunks_fa<-DNAStringSet(region_chunk_seqs$sequence)
names(region_chunks_fa)<-region_chunk_seqs$SeqID


#get the full region sequences
region_seq_fa<-readDNAStringSet(paste0("regulation/fastas/separated_regions/TRNP1_",region_to_get,"_seqs.fa"))


#blast (only chunks from same species to their original sequences)
blast_region_chunks<-list()
for (i in names(region_seq_fa)){
  sequence<-region_seq_fa[grep(paste0(i), names(region_seq_fa))]
  chunks<-region_chunks_fa[grep(paste0(i), names(region_chunks_fa))]
  dir<-tempdir()
  writeXStringSet(sequence, filepath = file.path(dir,"seqs.fasta"))
  db<-makeblastdb(file.path(dir, "seqs.fasta"), dbtype = "nucl")
  res<-blast(file.path(dir, "seqs.fasta"),type="blastn")
  blastn_region_seqs<-stats::predict(res,chunks, BLAST_args="-task blastn-short -word_size 10 -soft_masking false -dust no")
  blast_region_chunks[[i]]<-blastn_region_seqs
}


blast_region_chunks_df<-do.call(rbind.data.frame, blast_region_chunks) %>%
  dplyr::filter(Perc.Ident>99) %>%
  group_by(QueryID) %>%
  dplyr::top_n(1,wt=Alignment.Length) %>%
  ungroup() %>%
  inner_join(median_activity_MPRA[,c("SeqID", "cell_line", "activity_per_bp")], by=c("QueryID"="SeqID")) %>%
  dplyr::select(-Mismatches, -Gap.Openings,-E, -Bits) %>%
  dplyr::mutate(cell_line_species=paste0(cell_line,"_",SubjectID))

print(paste("Lost chunks during blasting:", dim(blast_region_chunks_df)[1]/3-length(region_chunks_fa))) 


# 2) per base coverage ####

calculate_per_bp_activity<-function(blast_output){
  
  activity_df<-data.frame(cell_line=as.character(),
                                species=as.character(),
                                position=as.integer(),
                                SeqID=as.character(),
                                activity_per_bp=as.numeric())
  
  for (sp in unique(blast_output$cell_line_species)){  
    
    activity_table<-blast_output[blast_output$cell_line_species==paste0(sp),]
    
    for (i in min(activity_table$S.start): max(activity_table$S.end)){
      
      chunk_activity<-as.vector(activity_table$activity_per_bp[(activity_table$S.start<=i & activity_table$S.end>=i)])
      names(chunk_activity)<-as.vector(activity_table$QueryID[(activity_table$S.start<=i & activity_table$S.end>=i)])
      
      if (length(chunk_activity)>0){
        
        activity_df<-activity_df %>%
          dplyr::add_row(cell_line=unique(word(activity_table$cell_line_species, 1,1, sep = "_")),
                  species=unique(word(activity_table$cell_line_species,2,3, sep="_")), 
                  position=i, 
                  SeqID=names(chunk_activity), 
                  activity_per_bp=chunk_activity)
      }
      
    }}
  return(activity_df)
}

positional_activity_df<-calculate_per_bp_activity(blast_region_chunks_df)
saveRDS(positional_activity_df, paste0("regulation/data/MPRA/output/positional_activity/positional_activity_",region_to_get,".rds"))

print(paste("Lost sequences:", length(unique(blast_region_chunks_df$SubjectID))-length(unique(positional_activity_df$species)))) 


# 3) summarized activity ####

#first, select the relevant ones
low_length_reg_add<-low_length_regs %>%
  dplyr::select(cell_line,species, length,activity) %>%
  dplyr::rename(sum_activity=activity)

summarized_output<-positional_activity_df %>%
  group_by(cell_line,species,position) %>%
  dplyr::summarise(mean_activity_per_bp=mean(activity_per_bp)) %>%
  group_by(cell_line, species) %>%
  dplyr::summarise(sum_activity=sum(mean_activity_per_bp),
            length=length(mean_activity_per_bp)) %>%
  ungroup() %>%
  dplyr::bind_rows(low_length_reg_add)

saveRDS(summarized_output, paste0("regulation/data/MPRA/output/summarized_activity/summarized_activity_",region_to_get,".rds"))


# print elapsed time
new <- Sys.time() - old # calculate difference
print(new)

