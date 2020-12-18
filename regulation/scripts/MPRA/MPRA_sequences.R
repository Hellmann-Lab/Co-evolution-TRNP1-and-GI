

# STEPS -------------------------------------------------------------------


# Convert fastq.gz to fq

# Cut Adapter sequences

# Map reads to the closest species genome with ngm

# Convert ngm unsorted sam output files to mapping coordinates sorted bam file

# Trinity reference guided assembly with bam ngm files

# Blat Trinity output fasta files (contigs) against reference genome and intersect with rbh peak bedfiles of closest species

# further processing in R


# NGM ---------------------------------------------------------------------

# settings
dir="TRNP1/dna_resequencing/ngm_mapped_ucsc";
slurms="TRNP1/dna_resequencing/slurms";
outfile="_ngm.sam";

cat /TRNP1/dna_resequencing/ngm_mapped_ucsc/fqfiles_genomes.txt | while read left right name genome_ncbi genome_ucsc;

do
# MAKING THE HEADER
echo '#!/bin/bash' >$name.sh  
echo '#SBATCH -n 1' >>$name.sh 
echo '#SBATCH --error='$name'.%J.err' >>$name.sh 
echo '#SBATCH --output='$name'.%J.out' >>$name.sh
#echo '#SBATCH --workdir='$slurms'' >> $name.sh
echo '#SBATCH --cpus-per-task=8' >> $name.sh
#echo '#SBATCH  --mem-per-cpu=9000' >> $name.sh

# THE ACTUAL COMMANDS
#echo "srun mkdir $dir/$name$outdir" >> $name.sh
echo "srun ngm.4.12 -r $genome_ucsc -1 $left -2 $right -i 0.65 -t 8 -o $dir/$name$outfile --skip-mate-check" >> $name.sh
# SUBMIT THE TASK 
sbatch $name.sh
done

# TRINITY -----------------------------------------------------------------

# settings
dir="TRNP1/dna_resequencing/trinity_refguided";
slurms="TRNP1/dna_resequencing/slurms";
outdir="_trinity";
mapdir="TRNP1/dna_resequencing/ngm_mapped_ucsc";
mapfile="_ngm.sort.bam";

cat /TRNP1/dna_resequencing/trinity_refguided/input_trinity.txt | while read left right name ncbi_genome uscsc_genome min_length;

do
# MAKING THE HEADER
echo '#!/bin/bash' >$name.sh  
echo '#SBATCH -n 1' >>$name.sh 
echo '#SBATCH --error='$name'.%J.err' >>$name.sh 
echo '#SBATCH --output='$name'.%J.out' >>$name.sh
#echo '#SBATCH --workdir='$slurms'' >> $name.sh
echo '#SBATCH --cpus-per-task=5' >> $name.sh
echo '#SBATCH  --mem-per-cpu=10000' >> $name.sh

# THE ACTUAL COMMANDS
echo "srun mkdir $dir/$name$outdir" >> $name.sh
echo "srun /home/vieth/progpackages/trinityrnaseq-2.0.6/Trinity --seqType fq --left $left --right $right --output $dir/$name$outdir  --min_contig_length $min_length --normalize_reads --CPU 5 --max_memory 10G --genome_guided_bam $mapdir/$name$mapfile --genome_guided_max_intron 1" >> $name.sh
# SUBMIT THE TASK 
sbatch $name.sh
done

# BLAT --------------------------------------------------------------------

# settings
oldfqdir="TRNP1/dna_resequencing/original_files";
fqdir="TRNP1/dna_resequencing/fqfiles";
mapdir="TRNP1/dna_resequencing/ngm_mapped_ucsc";
trinitydir="TRNP1/dna_resequencing/trinity_refguided";
slurms="TRNP1/dna_resequencing/slurms";
outdir="_trinity";
mapoutfile="_ngm.sam";
sortoutfile="_ngm.sort.bam";

cat TRNP1/dna_resequencing/scripts/trinity_refguided_contig_blat_input.txt | while read genome peakfile name;

do
# MAKING THE HEADER
echo '#!/bin/bash' >$name.sh  
echo '#SBATCH -n 1' >>$name.sh 
echo '#SBATCH --error='$name'.%J.err' >>$name.sh 
echo '#SBATCH --output='$name'.%J.out' >>$name.sh
echo '#SBATCH --workdir='$slurms'' >> $name.sh
echo '#SBATCH --cpus-per-task=2' >> $name.sh
#echo '#SBATCH  --mem-per-cpu=9000' >> $name.sh

# THE ACTUAL COMMANDS
echo "srun --chdir=$trinitydir/$name$outdir/ blat $genome Trinity-GG.fasta Trinity-GG.psl -t=DNA -q=DNA -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0 -extendThroughN" >> $name.sh
echo "srun pslToBed $trinitydir/$name$outdir/Trinity-GG.psl $trinitydir/$name$outdir/Trinity-GG.blat.bed" >> $name.sh
echo "srun  bedtools intersect -a $peakfile -b $trinitydir/$name$outdir/Trinity-GG.blat.bed -wao > $trinitydir/$name$outdir/Trinity-GG.blat.peak.intersect.bed" >> $name.sh
# SUBMIT THE TASK 
sbatch $name.sh
done


# SEQUENCE COLLECTION -----------------------------------------------------


# Peak bed files ----------------------------------------------------------

# get file names
allfiles <- list.files(path="TRNP1/output/peak_bed_fa", recursive = T, full.names=T, include.dirs = T, pattern="*.bed$")
peakfiles <- allfiles[grepl(pattern = "*50*", allfiles)==F]
peakfiles
# get species names
sc_species <- (t(sapply(strsplit(sapply(strsplit(peakfiles, "/"), "[[", 4), "_"), function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species
peaks.bed <- lapply(peakfiles, function(x) read.table(x, stringsAsFactors = F, header=F, sep="\t"))
names(peaks.bed) <- sc_species


# Contig peak identification ----------------------------------------------

# blat peak contig overlap bed files
# once with species specific peak blat
blatpeak.coord.files <- list.files(path="TRNP1/dna_resequencing/trinity_refguided", recursive = T, full.names=T, include.dirs = T, pattern="Trinity-GG.blat.peak.intersect.bed")
# trinity fasta files
trinity.fa.files <- list.files(path="TRNP1/dna_resequencing/trinity_refguided", recursive = T, full.names=T, include.dirs = T, pattern="Trinity-GG.fasta")
# load in contig fasta files by trinity in a list
contig.fas <- lapply(trinity.fa.files, function(x) readDNAStringSet(x))
# name fasta by secies
names(contig.fas) <- paste0(sapply(strsplit(trinity.fa.files, "/"), "[[", 8))

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
names(blatpeak.coord) <- paste0(sapply(strsplit(blatpeak.coord.files, "/"), "[[", 8))
blatpeak.coord <- lapply(blatpeak.coord, function(x) {x <- dplyr::mutate(x, peak.length=V3-V2, contig.blat.length=V7-V6)})
# some of the contigs have unrealistic lengths (far longer than possible pcr amplicons) and horrible long blat lengths, discard them!
blatpeak.coord <- lapply(blatpeak.coord, function(x) {x <- dplyr::filter(x, contig.blat.length<4000)})

View(blatpeak.coord[["Pan_troglodytes_trinity"]])

# combine blat peak intersect results with contig lengths
all.coord <- blatpeak.coord
for (i in names(blatpeak.coord)) {
  all.coord[[i]] <- dplyr::inner_join(blatpeak.coord[[i]], contig.fas.lengths[[i]], by=c("V8"="contig.name"))
}

# retain per species only one contig per peak location (longest contig)
all.coord <- lapply(all.coord, function(x) {x <- x %>%  dplyr::group_by(V4) %>% dplyr::top_n(n = 1, contig.blat.length) %>% dplyr::distinct(V4)})

# make a recoder input for each species based on peak contig match
recode.contig <- lapply(all.coord, function(x) {x <- x %>% dplyr::select(V4, V8) %>% dplyr::rename(peak.name=V4, contig.name=V8) } )

# rename fasta headers of contig fa files to peak names
recode.contig.peak <- recode.contig
for (i in names(recode.contig.peak)) {
  recode.contig.peak[[i]] <- dplyr::left_join(data.frame(contig.name=names(contig.fas[[i]]), stringsAsFactors = F), recode.contig.peak[[i]], by="contig.name")
  recode.contig.peak[[i]] <- na.omit(recode.contig.peak[[i]])
  recode.contig.peak[[i]] <- dplyr::filter(recode.contig.peak[[i]], peak.name!="prom_exon1") # prom_exon1 amplicons not yet in the mix!
}

contig.peak.fas <- contig.fas
for (i in names(contig.peak.fas)) {
  contig.peak.fas[[i]] <- contig.peak.fas[[i]][recode.contig.peak[[i]]$contig.name]
  names(contig.peak.fas[[i]]) <- recode.contig.peak[[i]]$peak.name
  tmp <- factor(names(contig.peak.fas[[i]]), levels=c("upstream1", "upstream2","upstream3", "prom_exon1", "intron", "exon2", "downstream"))
  contig.peak.fas[[i]] <- contig.peak.fas[[i]][sort(tmp)]
}

# change sequence in human to all plus strand, otherwise multiple alignment will fail! (can give directionality in mafft and he will take care by making reverse seq automatically)
all.coord[["Homo_sapiens_trinity"]][,"V10"] #(sequence 3, and 6 are on minus -> turn them on plus)
contig.peak.fas[["Homo_sapiens_trinity"]][3] <- reverseComplement(contig.peak.fas[["Homo_sapiens_trinity"]][3])
contig.peak.fas[["Homo_sapiens_trinity"]][6] <- reverseComplement(contig.peak.fas[["Homo_sapiens_trinity"]][6])

# append species name to seq name to combine it to long dnastringset
contig.peakspecies.fas <- contig.peak.fas[unlist(lapply(contig.peak.fas, length))>0] # kick out species with no seq
for (i in names(contig.peakspecies.fas)) {
  names(contig.peakspecies.fas[[i]]) <- paste(i, names(contig.peakspecies.fas[[i]]), sep="_")
}
dna_list <- contig.peakspecies.fas
names(dna_list) <- NULL # otherwise do.call() will return the input list
all.contigs.fas <- do.call(c, dna_list)

# continue with species-specific blat

# peak contigs per species, fasta file split by peak
contig.species.fas <- contig.peakspecies.fas

# peak contigs per species (only full sets!), fasta file split by peak
contigall.species.fas <- contig.peak.fas[unlist(lapply(contig.peak.fas, length))==6]
# peak contigs per species (only full sets!), concatenated fasta without split by peak
contigallconcat.species.fas <- lapply(contigall.species.fas, function(x) unlist(x))

# seq of contig per peak, make human first sequence
all.upstream1.fas <- all.contigs.fas[grepl("upstream1",names(all.contigs.fas))]
names(all.upstream1.fas)
all.upstream1.fas <- all.upstream1.fas[c(10,1:9,11:27)]

all.upstream2.fas <- all.contigs.fas[grepl("upstream2",names(all.contigs.fas))]
names(all.upstream2.fas)
all.upstream2.fas <- all.upstream2.fas[c(10,1:9,11:28)]

all.upstream3.fas <- all.contigs.fas[grepl("upstream3",names(all.contigs.fas))]
names(all.upstream3.fas)
all.upstream3.fas <- all.upstream3.fas[c(10,1:9,11:22)]

all.intron.fas <- all.contigs.fas[grepl("intron",names(all.contigs.fas))]
names(all.intron.fas)
all.intron.fas <- all.intron.fas[c(12,1:11,13:31)]

all.exon2.fas <- all.contigs.fas[grepl("exon2",names(all.contigs.fas))]
names(all.exon2.fas)
all.exon2.fas <- all.exon2.fas[c(13,1:12,14:32)]

all.downstream.fas <- all.contigs.fas[grepl("downstream",names(all.contigs.fas))]
names(all.downstream.fas)
all.downstream.fas <- all.downstream.fas[c(13,1:12,14:33)]

# # write fasta seq to files (NULL output = writing to disk is done)
sapply(names(all.contigs.fas), function(x) {
  y <- DNAStringSet(all.contigs.fas[[x]])
  names(y) <- x
  writeXStringSet(y, filepath =paste("TRNP1/output/trinity_contigs/per_species_per_peak/","",x, ".fa", sep=""), format = "fasta", append = F)
}  )
sapply(names(contig.species.fas), function(x) {
  writeXStringSet(contig.species.fas[[x]], filepath =paste("TRNP1/output/trinity_contigs/per_species/","",x, ".fa", sep=""), format = "fasta", append = F)
}  )
sapply(names(contigall.species.fas), function(x) {
  writeXStringSet(contigall.species.fas[[x]], filepath =paste("TRNP1/output/trinity_contigs/per_species/","", x, "_fullsetofpeaks.fa", sep=""), format = "fasta", append = F)
}  )


# Primate RBH sequences ---------------------------------------------------


# load in the results from previous rbh runs: extend by 500 to have adjacent sequences, extract sequences from fasta genomes
# dhs bed files
primate.bed.files <- list.files(path="TRNP1/output/peak_bed_fa/", pattern="*.bed$", recursive = T, full.names=T, include.dirs = T)
primate.bed.files <- primate.bed.files[grepl(pattern = "*500*",primate.bed.files)==F]
primate.bed.files
Species <- (t(sapply(strsplit(sapply(strsplit(primate.bed.files, "//"), "[[", 2), "_"),  function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species
Genome_short <- sapply(strsplit( c(sapply(strsplit(sapply(strsplit(primate.bed.files, "//"), "[[", 2), "_"),  function(x) tail(x, n=1))), "[.]"), "[[", 1)
primate.genomes <- paste0("/data/ngs/genomes/UCSC/", Species, "/", Genome_short, "*.fa")
primate.genomes 
primate.filename <- c(t(sapply(strsplit(sapply(strsplit(primate.bed.files, "//"), "[[", 2), "_"),  function(x) head(x, n=3))) %>% data.frame() %>% tidyr::unite(col=filename, X1:X3, sep="_"))$filename
primate.filename
primate.fa.files <- paste0("TRNP1/output/peak_bed_fa/",primate.filename,"_50.fa")
primate.fa.files
primate.bed.outfiles <- paste0("TRNP1/output/peak_bed_fa/",primate.filename,"_50.bed")
primate.bed.outfiles

for (i in 1:length(primate.filename)) {
  message(primate.filename[i])
  # make genome length file
  system(paste0("/usr/bin/fastalength ", primate.genomes[i]," > /TRNP1/output/peak_bed_fa/tmp"))
  system(paste0(' sed \'1i size chrom\' /TRNP1/output/peak_bed_fa/tmp > /TRNP1/output/peak_bed_fa/tmp2 '))
  system(paste0('awk -v OFS=\'\t\' \'{ print $2 "\t " $1}\' /TRNP1/output/peak_bed_fa/tmp2 > /TRNP1/output/peak_bed_fa/mygenome'))
  #bedtools slop extension
  system(paste0("/usr/bin/bedtools slop -b 50 -i ",primate.bed.files[i], " -g /TRNP1/output/peak_bed_fa/mygenome > ",primate.bed.outfiles[i]))
  # get the sequences
  system(paste0("/usr/bin/bedtools getfasta -name -fi ",primate.genomes[i]," -bed ",primate.bed.outfiles[i], " -fo ", primate.fa.files[i]))
}

rbh.primates.seq.files <- primate.fa.files
names(rbh.primates.seq.files) <- primate.filename
# load in fasta files
rbh.seqs.primate <- lapply(names(rbh.primates.seq.files), function(x) {
  y <- readDNAStringSet(rbh.primates.seq.files[x], format = "fasta",use.names = T)
  names(y) <- paste(x, names(y), sep="_")
  y
})

rbh.seqs.primate.set <- DNAStringSet(do.call(c, unname(unlist(rbh.seqs.primate))))
names(rbh.seqs.primate.set)


# Annotating trinity sequences --------------------------------------------

# blat
queryfiles <- list.files(path="TRNP1/output/trinity_contigs/per_species/", recursive = T, full.names=T, include.dirs = T)
queryfiles <- queryfiles[grepl(pattern = "*fullsetofpeaks*",queryfiles)==F]
queryfiles
peaks.fa <- list.files(path="TRNP1/output/peak_bed_fa/", recursive = T, full.names=T, include.dirs = T, pattern="*_50.fa$")

blatinfo <- read.table(file = "TRNP1/dna_resequencing//trinity_refguided/input_blat_contig", sep="\t")
fatmp <- paste(substr(blatinfo$V2, start = 1, stop=nchar(blatinfo$V2)-12), "50.fa", sep="_")
fatmp[13] <- "TRNP1/output/peak_bed_fa/Homo_sapiens_peaks_50.fa" # because the human genome is the only one not following the naming pattern
subjectfiles <- fatmp
names(subjectfiles) <- paste(blatinfo$V3, "trinity", sep="_")
filename = substr(sapply(strsplit(queryfiles, "//"), "[[", 2), start=1, stop=nchar(sapply(strsplit(queryfiles, "//"), "[[", 2))-3)
names(queryfiles) <- filename
output.blat = paste0("blat/blatresult_",filename,".psl")
names(output.blat) = filename

for (i in filename) {
  print(i)
  #blat
  system(paste0("/opt/bin/blat ",subjectfiles[i], " ", queryfiles[i]," ", output.blat[i]," -t=dna", " -q=dna", " minIdentity=10", " -out=psl", " -minScore=1"))
}

blatresult <- vector("list", length(queryfiles))
names(blatresult) <- filename
for (i in filename) {
  blat_table_names <- c("matches", "mismatches", "repmatches", "Ncount", "qNumInsert", "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize", "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount", "blockSizes", "qStarts", "tStarts")
  col_Classes.blat <- c(rep("numeric",8), rep("character",2), rep("numeric",3), "character", rep("numeric",4), rep("character",3))
  hit_table.blat <- tryCatch({
    if (file.size(output.blat[i]) > 427){
      data.table::fread( input  = output.blat[i], sep = "\t", header = FALSE, colClasses =col_Classes.blat, skip=5, verbose = F)
    } 
    else {
      data.table(t(c(rep(0,8), rep("empty",2), rep(0,3), "empty", rep(0,4), rep("empty",3))))
    }
  })
  data.table::setnames( x = hit_table.blat, old = paste0("V",1:length(blat_table_names)), new = blat_table_names)
  blatresult[[i]] <- hit_table.blat
}

checkup.blat <- do.call(rbind, lapply(blatresult, dim))

# make a list of sequence lengths
contig.lengths <- sapply(queryfiles, function(x) {
  fasta.seqlengths(x, use.names = T)
})
checkup.lengths <- do.call(rbind, lapply(contig.lengths, length))

subject.lengths <- sapply(subjectfiles, function(x) {
  fasta.seqlengths(x, use.names = T)
})

subject.alphabetfrq <- lapply(subjectfiles, function(x) {
  y <- readDNAStringSet(x)
  alphabetFrequency(y)
})

# kick out empty objects
kickout.blat <- sapply(blatresult, function(x) { any(grepl("empty", x)) })
blatresult <- blatresult[!kickout.blat]

# calculate the new query start end by adding 50 bases whenever possible
recoord.blat <- blatresult
for (i in names(blatresult)) {
  recoord.blat[[i]] <- blatresult[[i]] %>% mutate(qStart_new=qStart-50, qEnd_new=qEnd+50, qAlign=qEnd-qStart, tAlign=tEnd-tStart) %>% mutate(qFrac=qAlign/qSize*100, tFrac=tAlign/tSize*100) %>% dplyr::mutate(qStart_final=ifelse(qStart_new<0, 1, qStart_new), qEnd_final=ifelse(qEnd_new>qSize, qSize, qEnd_new)) %>% dplyr::mutate(qStart_FINAL=ifelse(qSize<tSize, 1, qStart_final), qEnd_FINAL=ifelse(qSize<tSize, qSize, qEnd_final-1)) %>% dplyr::mutate(qStart_FINAL=ifelse(qStart_FINAL==0, 1, qStart_FINAL)) %>% dplyr::group_by(qName) %>% dplyr::top_n(n = 1, tAlign) %>% dplyr::distinct(qName) 
}

recoord.blat.df <- do.call("rbind", recoord.blat)
openxlsx::write.xlsx(recoord.blat.df, file="recoord_blat_trinity_new.xlsx") # looks fine
# change the fasta files accordingly 
# make a list of sequences
contig.seqs <- lapply(names(recoord.blat), function(x) {
  y <- readDNAStringSet(queryfiles[x], format = "fasta",use.names = T)
  z <- y[names(y) %in% recoord.blat[[x]]$qName]
  tmp <- DNAStringSet(z, start=recoord.blat[[x]][,qStart_FINAL], end=recoord.blat[[x]][,qEnd_FINAL], use.names = T)
  tmp
})
names(contig.seqs) <- names(recoord.blat)

contig.seqs.set <- DNAStringSet(do.call(c, unname(unlist(contig.seqs))))
names(contig.seqs.set)


# Mammalian RBH sequences -------------------------------------------------

# rbh overlap bed files (best hit from human to other mammals, all possible blat hits for those in the human and intersected with real dhs sites)
rbh.coord.files <- list.files(path="/TRNP1/output/rbh", recursive = T, full.names=T, include.dirs = T, pattern="*.intersect.bed")
# read in rbh bed files
rbh.coord <- lapply(rbh.coord.files , function(x) read.table(x, header=F, sep="\t", stringsAsFactors = F))
names(rbh.coord) <- paste0(sapply(strsplit(rbh.coord.files, "/"), "[[", 8))
rbh.coord <- lapply(rbh.coord, function(x) {x <- dplyr::mutate(x, peak.length=V3-V2, rbh.length=V7-V6)})

# some of the rbhs have unrealistic lengths and horrible long blat lengths, discard them!
rbh.coord <- lapply(rbh.coord, function(x) {x <- dplyr::filter(x, rbh.length<10000)})

# retain per species only one rbh per peak location (the one that is the closest in length to the real peak length)
all.rbh.coord <- lapply(names(rbh.coord), function(x) {
  y <- rbh.coord[[x]]
  z <- y %>%  dplyr::mutate(peak.rbh.diff=abs(rbh.length-peak.length), rbhname=paste(x, "rbh", V4, sep="_")) %>% dplyr::group_by(V4) %>% dplyr::arrange(peak.rbh.diff) %>% dplyr::slice(1) %>% ungroup() %>% filter(V8!=".")
  z
})
names(all.rbh.coord) <- names(rbh.coord)

# get the sequences from genomes
# make a bed file of rbh to extract fasta files from genomes (remove mRNA location)
all.rbh.coord.bed <- lapply(all.rbh.coord, function(x) {
  y <- x
  y[,"Chromosome"] <- sapply(strsplit(y$V8, split = "[:]"), "[[", 1)
  tmp <- sapply(strsplit(y$V8, split = "[:]"), "[[", 2)
  y[,"Start"] <- sapply(strsplit(tmp, split = "[-]"), "[[", 1)
  y[,"End"] <- sapply(strsplit(tmp, split = "[-]"), "[[", 2)
  z <- y[,c("Chromosome", "Start", "End", "rbhname")]
  z <- z[grepl(z$rbhname, pattern="*mRNA")==F,]
  z
})
# save it as a physical bed file
sapply(names(all.rbh.coord.bed), function(x) {
  write.table(all.rbh.coord.bed[[x]], file =paste("TRNP1/output/rbh_bedfiles/","",x, ".bed", sep=""), col.names = F, row.names = F, sep = "\t", quote = F)
}  )


# dhs bed files
mammal.bed.files <- list.files(path="TRNP1/output/rbh_bedfiles/", pattern="*.bed$", recursive = T, full.names=T, include.dirs = T)
names(mammal.bed.files) <- substr(sapply(strsplit(mammal.bed.files, split = "//"), "[[", 2), start=1, stop=nchar(sapply(strsplit(mammal.bed.files, split = "//"), "[[", 2))-4)
mammal.genomes <- read.table(file = "/data/share/ngs/genomes/UCSC/MammalGenomes", sep="\t", stringsAsFactors = F)
mammal.genomes.fafiles <- paste0("/data/ngs/genomes/UCSC/", mammal.genomes$V2)
names(mammal.genomes.fafiles) <- mammal.genomes$V1
mammal.fa.outfiles <- paste0("TRNP1/output/rbh_fafiles/",mammal.genomes$V1,".fa")
names(mammal.fa.outfiles) <- mammal.genomes$V1

for (i in names(mammal.bed.files)) {
  message(i)
  # get the sequences
  system(paste0("/usr/bin/bedtools getfasta -name -fi ",mammal.genomes.fafiles[i]," -bed ",mammal.bed.files[i], " -fo ", mammal.fa.outfiles[i]))
}

rbh.seq.mammal.files <- mammal.fa.outfiles
# load in fasta files
rbh.seqs.mammal <- lapply(names(rbh.seq.mammal.files), function(x) {
  y <- readDNAStringSet(rbh.seq.mammal.files[x], format = "fasta",use.names = T)
  y
})

rbh.seqs.mammal.set <-  DNAStringSet(do.call(c, unname(unlist(rbh.seqs.mammal))))


# Union of RBH sequences --------------------------------------------------

# combine the rbh runs (primates form previous have priority as they are more congruent and contain more sequences)
rbh.seqs.mammal.set
rbh.seqs.primate.set

rbh.mammal <- data.frame(SeqNames=names(rbh.seqs.mammal.set), Species=(t(sapply(strsplit(names(rbh.seqs.mammal.set), "_"),function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species, PeakLocation=sapply(strsplit(names(rbh.seqs.mammal.set), "_"),function(x) tail(x, n=1)), stringsAsFactors = F)
rbh.mammal[,"CompareID"] <- paste(rbh.mammal[,"Species"], rbh.mammal[,"PeakLocation"],sep="_")
rbh.primate <- data.frame(SeqNames=names(rbh.seqs.primate.set), Species=(t(sapply(strsplit(names(rbh.seqs.primate.set), "_"),function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species, PeakLocation=sapply(strsplit(names(rbh.seqs.primate.set), "_"),function(x) tail(x, n=1)), stringsAsFactors = F)
rbh.primate[,"CompareID"] <- paste(rbh.primate[,"Species"], rbh.primate[,"PeakLocation"],sep="_")
rbh.seqs.combine.set <- dplyr::anti_join(rbh.mammal, rbh.primate, by="CompareID")
rbh.seqs.mammal.set.red <- rbh.seqs.mammal.set[names(rbh.seqs.mammal.set) %in% rbh.seqs.combine.set$SeqNames]

rbh.seqs.final.set <- c(rbh.seqs.primate.set, rbh.seqs.mammal.set.red)


# Union of RBH and Trinity sequences --------------------------------------

rbh.seqs.final.set
contig.seqs.set

rbh.set <- data.frame(SeqNames=names(rbh.seqs.final.set), Species=(t(sapply(strsplit(names(rbh.seqs.final.set), "_"),function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species, PeakLocation=sapply(strsplit(names(rbh.seqs.final.set), "_"),function(x) tail(x, n=1)), stringsAsFactors = F)
rbh.set[,"CompareID"] <- paste(rbh.set[,"Species"], rbh.set[,"PeakLocation"],sep="_")
head(rbh.set)

contig.set <- data.frame(SeqNames=names(contig.seqs.set), Species=(t(sapply(strsplit(names(contig.seqs.set), "_"),function(x) head(x, n=2))) %>% data.frame() %>% tidyr::unite(col=Species, X1:X2, sep="_"))$Species, PeakLocation=sapply(strsplit(names(contig.seqs.set), "_"),function(x) tail(x, n=1)), stringsAsFactors = F)
contig.set[,"CompareID"] <- paste(contig.set[,"Species"], contig.set[,"PeakLocation"],sep="_")
head(contig.set)

rbh.contig.combine.set <- dplyr::semi_join(rbh.set, contig.set, by="CompareID")

rbh.seqs.final.set.red <- rbh.seqs.final.set[!(names(rbh.seqs.final.set) %in% rbh.contig.combine.set$SeqNames)]

rbh.contig.final.set <- c(contig.seqs.set, rbh.seqs.final.set.red)
sort(names(rbh.contig.final.set))
table(grepl(names(rbh.contig.final.set), pattern = "*trinity*"))

# split by region
peakname <- c("upstream1", "upstream2", "upstream3", "intron", "exon2", "downstream", "prom_exon1")
peakpattern <- c("*upstream1$", "*upstream2$", "*upstream3$", "*intron$", "*exon2$", "*downstream$", "*prom_exon1$")
peaks <- data.frame(peakpattern=peakpattern, peakname=peakname, stringsAsFactors = F)
rbh.contig.all.seqs <- vector("list", length(peakname))
names(rbh.contig.all.seqs) <- peakname

for ( i in peaks$peakname ) {
  rbh.contig.all.seqs[[i]] <- rbh.contig.final.set[grep(pattern=peaks[peaks$peakname==i, "peakpattern"],names(rbh.contig.final.set))]
}
rbh.contig.all.seqs

# i realized after running msa once, that a few species had far too long sequences, remove them if there is no alternative
tmp.stat <- lapply(rbh.contig.all.seqs, function(x) {
  y <- data.frame(seqLength=width(x), seqNames=names(x))
  y
})
sapply(tmp.stat, summary)
tmp.stat[["intron"]]
needToRemove <- c("Microcebus_murinus_peaks_upstream1", "Otolemur_garnettii_peaks_upstream1", "Bos_taurus_rbh_intron", "Cavia_porcellus_rbh_intron", "Dipodomys_ordii_rbh_intron", "Ochotona_princeps_rbh_intron", "Ovis_aries_rbh_intron", "Rattus_norvegicus_rbh_intron", "Trichechus_manatus_rbh_intron","Cricetulus_griseus_rbh_intron",  "Cavia_porcellus_rbh_exon2", "Ceratotherium_simum_rbh_exon2", "Dasypus_novemcinctus_rbh_exon2", "Equus_ferus_caballus_rbh_exon2", "Felis_catus_rbh_exon2", "Mustela_putorius_furo_rbh_exon2", "Myotis_lucifugus_rbh_exon2", "Sus_scrofa_rbh_exon2", "Tursiops_truncatus_rbh_exon2", "Vicugna_pacos_rbh_exon2")

# remove long sequences before the split into peak pattern, quick and dirty
rbh.contig.final.set <- rbh.contig.final.set[!(names(rbh.contig.final.set) %in% needToRemove)]

for ( i in peaks$peakname ) {
  rbh.contig.all.seqs[[i]] <- rbh.contig.final.set[grep(pattern=peaks[peaks$peakname==i, "peakpattern"],names(rbh.contig.final.set))]
}
rbh.contig.all.seqs

tmp.stat <- lapply(rbh.contig.all.seqs, function(x) {
  y <- data.frame(seqLength=width(x), seqNames=names(x))
  y
})
sapply(tmp.stat, summary)


# Murine DHS --------------------------------------------------------------

# open chromatin (there are three expr 2 for E14.5 (SRX188655, SRX191055) and 1 for E18.5 (SRX191042))
mouse.dhs.files <- read.table(file="TRNP1/input/mm10.fetalbrain.dhs.files", header = F, stringsAsFactors = F)
mouse.dhs.files <- read.table(file="TRNP1/input/mm10.fetalbrain.dhs.files", header = F, stringsAsFactors = F)
mouse.dhs.files<- mouse.dhs.files$V1
epigenome.mm10.dhs <- lapply(mouse.dhs.files, function(x) read.table(paste("" , x, sep="/"), stringsAsFactors = F, as.is = T, header = F))
epigenome.mm10.dhs <- lapply(mouse.dhs.files, function(x) read.table(x, stringsAsFactors = F, as.is = T, header = F))
names(epigenome.mm10.dhs) <- c("SRX188655", "SRX191042", "SRX191055")
# convert to genomic ranges objects
epigenome.mm10.dhs.gr <- lapply(epigenome.mm10.dhs, function(x) {
  GRanges(seqnames = Rle(x[,1]),
          ranges=IRanges(start=x[,2], end=x[,3]),
          #strand=x[,4], # hotspot does not keep track of that, peakdeck does!
          #name=x[,6], # hotspot does not give names to peaks, peakdeck does!
          score=x[,5])
})

# extract nucleotide sequences of mouse DHS
epigenome.mm10.dhs.seq <- lapply(epigenome.mm10.dhs.gr, function(x) {
  getSeq(Mmusculus, x)
})

# calculate background frequency of a c g t in dhs peaks
mm10.dhs.bg <- colSums(alphabetFrequency(IRanges::unlist(DNAStringSetList(epigenome.mm10.dhs.seq))))[1:4]
mm10.dhs.bg  <- mm10.dhs.bg/sum(mm10.dhs.bg); mm10.dhs.bg

# mouse TRNP1 locus chr4:133,488,061-133,490,243
mm10.trnp1.locus <- GRanges(seqnames=Rle("chr4"), ranges=IRanges(start=133488000, end=133518000), strand=Rle("*"))

# DHS overlapping mouse TRNP1 locus
trnp1.mm10.dhs.gr <- lapply(epigenome.mm10.dhs.gr, function(x) {
  subsetByOverlaps(x, mm10.trnp1.locus) 
})
trnp1.mm10.dhs.gr  <- GenomicRanges::unlist(GRangesList(trnp1.mm10.dhs.gr))
# combine to one GRanges object
trnp1.mm10.dhs.gr  <- reduce(trnp1.mm10.dhs.gr ); trnp1.mm10.dhs.gr 
# pad seq to avoid edge effects
trnp1.mm10.dhs.gr <- GRanges(seqnames = seqnames(trnp1.mm10.dhs.gr), ranges=IRanges(start=start(trnp1.mm10.dhs.gr)-50, end=end(trnp1.mm10.dhs.gr)+50))
trnp1.mm10.dhs.seq <- getSeq(Mmusculus, trnp1.mm10.dhs.gr)
trnp1.mm10.dhsconcat.seq <- unlist(trnp1.mm10.dhs.seq)

# combine the two in one granges and one seq object
trnp1.mm10.all.seq <- c(DNAStringSet(trnp1.mm10.dhs.seq))
mousedhs <- c("intron", "unique1", "prom_exon1","unique2", "unique3","upstream3", "upstream2", "upstream1")
names(trnp1.mm10.all.seq) <- paste("Mus_musculus_dhs", mousedhs, sep="_")

# FINAL SET OF SEQUENCES --------------------------------------------------

# the final set contains targeted resequenced primate species, rbh results, ancestral sequences, real mouse dhs sequences.
rbh.contig.final.set
# substitute the rbh mouse seq with the real dhs sequences!
rbh.contig.final.all.set <- rbh.contig.final.set[!(grepl(pattern = "Mus_musculus",names(rbh.contig.final.set)))]
all.seqs.set <- c(DNAStringSet(trnp1.mm10.all.seq), DNAStringSet(rbh.contig.final.all.set),DNAStringSet(node.peak.set))
