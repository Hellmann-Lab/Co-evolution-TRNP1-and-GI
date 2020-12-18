# This script contains the commands used to identify the longest contig per species, which covers the coding sequences. The input of this script has not been uploaded due to limited git space. But the output of this script can be found under protein/fastas/resequenced_promexon1_primates.fa


# SETTINGS ----------------------------------------------------------------

# load in needed packages
libs <- c("IRanges", "Biostrings","reshape2", "reshape", "GenomicRanges", "data.table", "tidyverse", "phangorn", "ape", "ggrepel", "ggtree", "caper", "phytools", "ggbio")

sapply(libs, require, character.only=T)
options(stringsAsFactors = F)


# PEAK BED FILES -------------------------------------------------------------------

# get file names -- these have NOT been uploaded due to size restrictions!
allfiles <- list.files(path="/data/share/htp/TRNP1/TRNP1_Beate/output/peak_bed_fa", recursive = T, full.names=T, include.dirs = T, pattern="*.bed$")
peakfiles <- allfiles[grepl(pattern = "*50*", allfiles)==F]
peakfiles
# get species names
sc_species <- (t(sapply(strsplit(sapply(strsplit(peakfiles, "/"), "[[", 6), "_"), function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species
peaks.bed <- lapply(peakfiles, function(x) read.table(x, stringsAsFactors = F, header=F, sep="\t"))
names(peaks.bed) <- sc_species
peaks.bed[1]
# CONTIG PEAK IDENTIFICATION ----------------------------------------------

# blat peak contig overlap bed files
# once with species specific peak blat
blatpeak.coord.files <- list.files(path="/data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/trinity_refguided/trinity_lucas", recursive = T, full.names=T, include.dirs = T, pattern="Trinity-GG.blat.peak.intersect.bed")
# trinity fasta files
trinity.fa.files <- list.files(path="/data/share/htp/TRNP1/TRNP1_exon_promoter_resequencing/trinity_refguided/trinity_lucas", recursive = T, full.names=T, include.dirs = T, pattern="Trinity-GG.fasta$")
# load in contig fasta files by trinity in a list
contig.fas <- lapply(trinity.fa.files, function(x) readDNAStringSet(x))
# name fasta by secies
names(contig.fas) <- paste0(sapply(strsplit(trinity.fa.files, "/"), "[[", 9))

# replace long trinity style fasta headers by shortened ones
for (i in 1:length(contig.fas)) {
  names(contig.fas[[i]]) <- sub(" .*", "", names(contig.fas[[i]]))
}

contig.fas.lengths <- contig.fas
for (i in 1:length(contig.fas.lengths)) {
  contig.fas.lengths[[i]] <- data.frame(contig.length=width(contig.fas.lengths[[i]]), contig.name=names(contig.fas[[i]]), stringsAsFactors = F)
}

# make recoder df based on blat contig
blatpeak.coord <- lapply(blatpeak.coord.files , function(x) read.table(x, header=F, sep="\t", stringsAsFactors = F))
names(blatpeak.coord) <- paste0(sapply(strsplit(blatpeak.coord.files, "/"), "[[", 9))
blatpeak.coord <- lapply(blatpeak.coord, function(x) {x <- dplyr::mutate(x, peak.length=V3-V2, contig.blat.length=V7-V6)})
# some of the contigs have unrealistic lengths (far longer than possible pcr amplicons) and horrible long blat lengths, discard them!
blatpeak.coord <- lapply(blatpeak.coord, function(x) {x <- dplyr::filter(x, contig.blat.length<4000)})

View(all.coord[["Pan_troglodytes_lucas_trinity"]])

# combine blat peak intersect results with contig lengths
all.coord <- blatpeak.coord
for (i in names(blatpeak.coord)) {
  all.coord[[i]] <- dplyr::inner_join(blatpeak.coord[[i]], contig.fas.lengths[[i]], by=c("V8"="contig.name"))
}

# retain per species only one contig per peak location (longest contig)
all.coord <- lapply(all.coord, function(x) {x <- x %>%  dplyr::group_by(V4) %>% dplyr::top_n(n = 1, contig.length) %>% dplyr::distinct(V4, .keep_all=T)})

# make a recoder input for each species based on peak contig match
recode.contig <- lapply(all.coord, function(x) {x <- x %>% dplyr::select(V4, V8) %>% dplyr::rename(peak.name=V4, contig.name=V8) } )

# rename fasta headers of contig fa files to peak names
recode.contig.peak <- recode.contig
for (i in names(recode.contig.peak)) {
  recode.contig.peak[[i]] <- dplyr::left_join(data.frame(contig.name=names(contig.fas[[i]]), stringsAsFactors = F), recode.contig.peak[[i]], by="contig.name")
  recode.contig.peak[[i]] <- na.omit(recode.contig.peak[[i]])
}

contig.peak.fas <- contig.fas
for (i in names(contig.peak.fas)) {
  contig.peak.fas[[i]] <- contig.peak.fas[[i]][recode.contig.peak[[i]]$contig.name]
  names(contig.peak.fas[[i]]) <- recode.contig.peak[[i]]$peak.name
  tmp <- factor(names(contig.peak.fas[[i]]), levels=c("upstream1", "upstream2","upstream3", "prom_exon1", "intron", "exon2", "downstream"))
  contig.peak.fas[[i]] <- contig.peak.fas[[i]][sort(tmp)]
}


# append species name to seq name to combine it to long dnastringset
contig.peakspecies.fas <- contig.peak.fas[unlist(lapply(contig.peak.fas, length))>0] # kick out species with no seq
for (i in names(contig.peakspecies.fas)) {
  names(contig.peakspecies.fas[[i]]) <- paste(i, names(contig.peakspecies.fas[[i]]), sep="_")
}
dna_list <- contig.peakspecies.fas
names(dna_list) <- NULL # otherwise do.call() will return the input list
all.contigs.fas <- do.call(c, dna_list)

#### NOTE I ran it only with species specific peak blat. since we have a lot of Ns in reference sequences, it might be good to do it once with human blat search!

# peak contigs per species, fasta file split by peak
contig.species.fas <- contig.peakspecies.fas

# seq of contig per peak, make human first sequence
all.prom_exon1.fas <- all.contigs.fas[grepl("prom_exon1",names(all.contigs.fas))]
names(all.prom_exon1.fas)
all.prom_exon1.fas <- all.prom_exon1.fas[c(6,1:5,7:18)]

# # write fasta seq to files 
writeXStringSet(all.prom_exon1.fas, filepath='protein/fastas/resequenced_promexon1_primates.fa', format = 'fasta', append=F)
