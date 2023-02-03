library(Biostrings)

#ok 1 very promising sequence (end), 1 a bit promising (start) --> let's check

ferseq<-readDNAStringSet("/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trinity_kmer_bwa/trinity-final.Trinity.fasta")
fullest<-ferseq[names(ferseq)=="TRINITY_DN51_c0_g1_i1 len=501 path=[0:0-500]"]
writeXStringSet(reverseComplement(fullest),"/data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trinity_kmer_bwa/best.fa")

system("/home/zane/TRNP1/mafft7/mafft-linux64/mafft.bat --adjustdirection --reorder --add /data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trinity_kmer_bwa/best.fa /data/share/htp/TRNP1/paper_data/protein/fastas/prank_output/prank_TRNP1_coding_30species_forCoevol_longer.best.nuc.fas > /data/share/htp/TRNP1/dna_resequencing/Trnp1_Ferret/trinity_kmer_bwa/best.aligned.fa")



#concatenate first contig and original seq

full_aln<-readDNAStringSet("/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding_seqs_30species_forCoevol_longer.fa")
fer1<-subseq(full_aln[names(full_aln)=="Mustelaput"],1,464)
red_aln<-full_aln[names(full_aln)!="Mustelaput"]

fullest_cut<-subseq(reverseComplement(fullest), start=109,end=322+3)
ferfull<-xscat(fer1,fullest_cut)
names(ferfull)<-"Mustelaput"

writeXStringSet(ferfull,"/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/Mustela_putorius_reseq.fa")




#include the stopcodon for genbank

fullest_cut_stop<-subseq(reverseComplement(fullest), start=109,end=322+3)
ferfull_stop<-xscat(fer1,fullest_cut_stop)
names(ferfull_stop)<-"Mustelaput"

writeXStringSet(ferfull_stop,"/data/share/htp/TRNP1/paper_data/protein/fastas/TRNP1_coding/Mustela_putorius_reseq_stopCodon.fa")


