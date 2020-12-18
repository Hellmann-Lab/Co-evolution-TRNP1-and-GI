install.packages("PhyloOrchard", repos="http://R-Forge.R-project.org")
BiocManager::install("ggtree")
install.packages("phylobase")

library(PhyloOrchard)
library(phylobase)
library(ggtree)
library(ggplot2)


#load in phylogenetic trees
data(BinindaEmondsEtAl2007)
bestmammaltree <- BinindaEmondsEtAl2007[[1]]
bestmammaltree
complete_library_meta_information <- read.delim("regulation/data/MPRA/input/complete_library_meta_information.txt")
head(complete_library_meta_information)

species.all <- unique(complete_library_meta_information$species)
species.all <- gsub(pattern = " ", replacement = "_", species.all)

# correct tree tip labels to own sample names and choose nearest species whereever possible
mammaltreenames <- data.frame(mammaltree.tips=bestmammaltree$tip.label, number=1:length(bestmammaltree$tip.label))
mynames <- data.frame(samplenames=sort(species.all))
diffnames <-  data.frame(diffnames=sort(setdiff(species.all,bestmammaltree$tip.label)))
diffnames
mammaltreenames[grepl(pattern="Aotus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2187] <- "Aotus_azarae"
mammaltreenames[grepl(pattern="Cercocebus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2075] <- "Cercocebus_chrysogaster"
mammaltreenames[grepl(pattern="Chlorocebus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2072] <- "Chlorocebus_sabeus"
mammaltreenames[grepl(pattern="Cricetulus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[1012] <- "Cricetulus_griseus"
mammaltreenames[grepl(pattern="Equus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2612] <- "Equus_ferus"
mammaltreenames[grepl(pattern="Felis", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2793] <- "Felis_catus"
mammaltreenames[grepl(pattern="Lagothrix", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2165] <- "Lagothrix_cana"
mammaltreenames[grepl(pattern="Hylobates", mammaltreenames$mammaltree.tips),] # all Nomascus species belong to the gneus hylobates in this tree, so change them in the tree file!
bestmammaltree$tip.label[2148] <- "Nomascus_gabriellae"
bestmammaltree$tip.label[2147] <- "Nomascus_leucogenys"
mammaltreenames[grepl(pattern="Pongo", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2137] <- "Pongo_abelii"
mammaltreenames[grepl(pattern="Sarcophilus", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[4385] <- "Sarcophilus_harrisii"
mammaltreenames[grepl(pattern="Vicugna", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2595] <- "Vicugna_pacos"
mammaltreenames[grepl(pattern="Macaca", mammaltreenames$mammaltree.tips),]
bestmammaltree$tip.label[2097] <- "Macaca_leonia"

mammaltreenames <- data.frame(mammaltree.tips=bestmammaltree$tip.label, number=1:length(bestmammaltree$tip.label))
mynames <- data.frame(samplenames=sort(species.all))
diffnames <-  data.frame(diffnames=sort(setdiff(species.all,bestmammaltree$tip.label)))
diffnames

str(bestmammaltree)

# add papio anubis as sister to papio hamadryas
papio <- grep(pattern="*Papio*", bestmammaltree$tip.label); papio

tipgroup <- c((papio-2):(papio+2))
p <- ggtree(bestmammaltree)
node<- phylobase::MRCA(bestmammaltree, tip=tipgroup)
p2 <- viewClade(p+geom_tiplab(), node=node) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + theme_tree2()
p <- ggtree(bestmammaltree) + theme_tree2()
p2 <- viewClade(p+geom_tiplab(), node=node) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + scale_x_continuous(limits=c(100,200))
p2

tmp.phylo <- bestmammaltree %>% tidytree::as_tibble() 
head(tmp.phylo)
papio_length <- unlist(tmp.phylo[grep(pattern = "Papio*", tmp.phylo$label),"branch.length"])

tip.papio<-list(edge=matrix(c(2,1),1,2),
                tip.label="Papio_anubis",
                edge.length=papio_length,
                Nnode=1)
anc.node <- 5289
class(tip.papio)<-"phylo"
mammaltree<-bind.tree(bestmammaltree, tip.papio, where=anc.node)
p <- ggtree(mammaltree)
testnode <- phylobase::MRCA(mammaltree, tip=c(2070:2090))
viewClade(p+geom_tiplab(), node=testnode) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

# add Symphalangus_syndactylus as sister to node ancestor of hylobates
symphalangus <- grep(pattern="Hylobates", mammaltree$tip.label); symphalangus
tmp.phylo <- mammaltree %>% tidytree::as_tibble() 
tmp.phylo[grep(pattern = "5328", tmp.phylo$parent),]
tmp.phylo[grep(pattern = "5330", tmp.phylo$parent),]
tip.symphalangus<-list(edge=matrix(c(2,1),1,2),
                       tip.label="Symphalangus_syndactylus",
                       edge.length=8,
                       Nnode=1)
anc.node <- 5330
class(tip.symphalangus)<-"phylo"
mammaltree<-bind.tree(mammaltree,tip.symphalangus, where=anc.node)
p <- ggtree(mammaltree)
testnode <- 5325
viewClade(p+geom_tiplab(), node=testnode) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

# add Chlorocebus_aethiops as sister to Chlorocebus_sabeus
chlorocebus <- grep(pattern="*Chlorocebus*", mammaltree$tip.label); chlorocebus
tmp.phylo <- mammaltree %>% tidytree::as_tibble() 
tmp.phylo[grep(pattern = "2072", tmp.phylo$node),]

tip.chlorocebus<-list(edge=matrix(c(2,1),1,2),
                      tip.label="Chlorocebus_aethiops",
                      edge.length=5.3,
                      Nnode=1)
anc.node <- 5282
class(tip.chlorocebus)<-"phylo"
mammaltree<-bind.tree(mammaltree,tip.chlorocebus, where=anc.node)
p <- ggtree(mammaltree)
testnode <- 5282
viewClade(p+geom_tiplab(), node=testnode) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

