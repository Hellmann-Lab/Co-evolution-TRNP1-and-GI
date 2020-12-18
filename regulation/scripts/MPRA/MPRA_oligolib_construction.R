
# SETTINGS ----------------------------------------------------------------

libs<-c("tidyverse","GenomicRanges","data.table","Biostrings","microRNA", 'reshape2')
sapply(libs, require, character.only=T)


system('mkdir -p output/mpra_chunk')

# READ IN SEQUENCES -------------------------------------------------------

all.seqs.set <- readDNAStringSet("regulation/data/MPRA/input/minlength.library.seqs.set.fa", format = 'fasta')

# SLIDING WINDOW ----------------------------------------------------------

# the final set contains targeted resequenced primate species, reciprocal best hits results, inferred ancestral sequences, real mouse and dhs sequences.
all.seqs.set
# make a sliding window per sequence entry, moving by 50 bases, for the once that are longer than 90 bases
minlength.seqs.set <- all.seqs.set[width(all.seqs.set)>94]
start.table <- list()
for (i in 1:length(minlength.seqs.set)) {
  seqwidth <- width(minlength.seqs.set[i])
  starts <- seq(from=1, to=seqwidth, by=40)
  ends <- starts+94
  real.ends <- ends[!ends>seqwidth]
  real.starts <- starts[!ends>seqwidth]
  laststart <- seqwidth-94
  real.starts <- c(real.starts, laststart)
  start.table[[i]] <- real.starts
}

all.seqs.chunks <- list()
for (i in 1:length(minlength.seqs.set)) {
  print(i)
  my_pos <- start.table[[i]]
  my_seq <- minlength.seqs.set[i]
  tmpchunk <- lapply(my_pos, function(pos) {
    subseq(x = my_seq, start = pos, width=94)
  })
  chunkset <- DNAStringSet(do.call(c, unname(unlist(tmpchunk))))
  names(chunkset) <- paste0(names(chunkset), "_chunk", seq(from = 1, to = length(chunkset), by=1))
  print(chunkset)
  # writeXStringSet(x = chunkset, filepath = paste0("mpra_chunks/", i, ".fa"), format = "fasta", append = F)
  all.seqs.chunks[[i]] <- chunkset
}

all.95chunks.seq <- DNAStringSet(do.call(c, unlist(unname(all.seqs.chunks))))

# all.95chunks.seq <- readDNAStringSet(filepath = "mpra_chunks/", format = "fasta", use.names = T)
alphabetFrequency(all.95chunks.seq)
all.seqs.95bases.set <-  c(DNAStringSet(all.95chunks.seq), DNAStringSet(all.seqs.set[!width(all.seqs.set)>94]))

# filter out the sequences that contain more than 10% N and thta are perfect duplicates!
acgt.count <- data.frame(alphabetFrequency(all.seqs.95bases.set), stringsAsFactors = F)
row.names(acgt.count) <- names(all.seqs.95bases.set) 
acgt.freq <- ( acgt.count / rowSums(acgt.count) ) *100
row.names(acgt.freq) <- names(all.seqs.95bases.set) 
acgt.freq.minN <- acgt.freq %>% tibble::rownames_to_column(var = "SeqID") %>% filter(N <5)

all.seqs.95bases.filt.set <- all.seqs.95bases.set[names(all.seqs.95bases.set) %in% acgt.freq.minN$SeqID]
all.seqs.95bases.filt.set <- unique(all.seqs.95bases.filt.set)

# save(all.seqs.95bases.filt.set, file = "95basechunks_oligolib.RData")

new.freq <- ( alphabetFrequency(all.seqs.95bases.filt.set) / rowSums(alphabetFrequency(all.seqs.95bases.filt.set)) ) * 100

head(new.freq)
summary(new.freq[,"N"])

# MPRA OLIGO LIBRARY ------------------------------------------------------

# the construct should be universal primer 1, enhancer variant seq, Kpn1 cut site, barcode seq, universal primer 2
uniprimer1 <- DNAString("ACTGGCCGCTTCACTG")
Kpn1ressite <- DNAString("GGTACCTCTAGA")
uniprimer2 <- DNAString("AGATCGGAAGAGCGTCG")

# barcode pool
# set how many sequences you want to produce
numOfSeqs <- 40000
# initialize empty object
seqs <- rep (NA, numOfSeqs)
# populate the object by shuffling and joining your bases
for (i in 1:numOfSeqs){
  # make a random sample of 10 * N bases
  tmpstring <- sample(DNA_ALPHABET[1:4], size=10, replace=TRUE)
  print(tmpstring)
  # filtering:
  # at least once an A , C, G and T
  if(length(match(c("A", "C", "G", "T"), tmpstring))<3) {
    next
  } else {
    # no kmer stretches of bases (AAAA, CCCC, TTTT, GGGG)
    tmpseq <- paste(tmpstring, collapse = '')
    kmerpattern <- "AAAA|TTTT|CCCC|GGGG"
    if(grepl(pattern=kmerpattern, tmpseq)) {
      next
    } else {
      # no KpnI recognition sequence : (5'-GGTACG-3', 3'-CCATGG-5')
      kpnI.rec <- c("GGTACG|CCATGG")
      if(grepl(pattern=kpnI.rec, tmpseq)) {
        next
      } else {
        print(tmpseq)
        seqs[i] <- tmpseq
      }
    } 
  }
}

table(is.na(seqs))

# kick out NAs
seqs <- seqs[!is.na(seqs)]
head(seqs)

# convert into DNAStringSet 
barcode.seqs <- DNAStringSet(seqs)
barcode.seqs
barcode.seqs.char <- as.character(barcode.seqs)
names(barcode.seqs.char) <- 1:length(barcode.seqs)

# no miRNA seed
# human
data(hsSeqs)
hSeedReg = seedRegions(hsSeqs)
comphSeed = as.character(reverseComplement(RNAStringSet(hSeedReg)))
comph = RNA2DNA(comphSeed)
mx.human = matchSeeds(comph, barcode.seqs.char)
mx.human.dat <- melt(mx.human)
miRNA.seq.hit.no.human <- unique(mx.human.dat$L2)
# mouse
data(mmSeqs)
mmSeedReg = seedRegions(mmSeqs)
compmmSeed = as.character(reverseComplement(RNAStringSet(mmSeedReg)))
compmm = RNA2DNA(compmmSeed)
mx.mouse = matchSeeds(compmm, barcode.seqs.char)
mx.mouse.dat <- melt(mx.mouse)
miRNA.seq.hit.no.mouse <- unique(mx.mouse.dat$L2)

# combine mouse and human miRNA hits on barcode seq
miRNA.seq.hit.no <- as.numeric(unique(c(miRNA.seq.hit.no.human, miRNA.seq.hit.no.mouse)))
length(miRNA.seq.hit.no)
head(miRNA.seq.hit.no)

# remove the hits from the barcode.seqs
barcode.seqs.final <- barcode.seqs[-c(miRNA.seq.hit.no)]
length(barcode.seqs.final)

# at least once an A , C, G and T 
tmp <- barcode.seqs.final
tmp.df <- alphabetFrequency(tmp)[,c(1:4)]
acgt.test <- apply(tmp.df, MARGIN = 1, function(x) all(x >= 1))
table(acgt.test)
barcode.seqs.final.acgt <- barcode.seqs.final[acgt.test]
tmp.test <- apply(alphabetFrequency(barcode.seqs.final.acgt)[,c(1:4)], MARGIN = 1, function(x) all(x >= 1))
table(tmp.test)
# remove duplicated entries
barcode.seqs.final.acgt.unique <- unique(barcode.seqs.final.acgt)
table(duplicated(barcode.seqs.final.acgt.unique))

# OLIGO LIBRARY -----------------------------------------------------------

# uniprimer1
uniprimer1.char <- as.character(uniprimer1)
uniprimer1.char <- rep(uniprimer1.char, length(all.seqs.95bases.filt.set))
# variable enhancer seq
enhancer.seq.char <- as.character(all.seqs.95bases.filt.set)
# kpnI cut site
Kpn1ressite.char <- as.character(Kpn1ressite)
Kpn1ressite.char <- rep(Kpn1ressite.char, length(all.seqs.95bases.filt.set))
# barcode seq
barcode.char <- barcode.seqs.final.acgt.unique[c(1:length(enhancer.seq.char))]
barcode.char <- as.character(barcode.char)

# uniprimer2
uniprimer2.char <- as.character(uniprimer2)
uniprimer2.char <- rep(uniprimer2.char, length(all.seqs.95bases.filt.set))

# make a data frame to save
oligo.lib.df <- data.frame(uniprimer1=uniprimer1.char, 
                           enhancer.seq.95=enhancer.seq.char, 
                           kpnIsite=Kpn1ressite.char, 
                           barcode=barcode.char,
                           uniprimer2=uniprimer2.char,
                           SeqID=names(enhancer.seq.char),
                           stringsAsFactors=F)
head(oligo.lib.df)
# save(oligo.lib.df, file = "oligolib_MPRA.RData")
# write.table(oligo.lib.df, file = "oligolib_df_MPRA.txt", quote = F, sep = "\t", col.names = T, row.names = T)
# OLIGO LIBRARY
oligo.lib <- paste(uniprimer1.char, enhancer.seq.char, Kpn1ressite.char, barcode.char, uniprimer2.char, sep="")
head(oligo.lib)
tail(oligo.lib)
length(oligo.lib)
# save to file 
# write(oligo.lib, file = "oligolib_MPRA.txt")

