#Plots for paper

#libraries
#load libraries 
libs<-c("tidyverse","cowplot","GenomicRanges","data.table","reshape2","ggrepel", "broom", "ape", "ggtree", "readr","geiger","nlme","phytools","grid","gtable","xtable", "wesanderson","stringi","aplot","ggdendro","flashClust")
sapply(libs, require, character.only=T)


#set working directory
setwd("/data/share/htp/TRNP1/paper_data/")

#a helper fuction for plotting
source("/data/share/htp/TRNP1/paper_data/scripts/new_aes.R")

#color schemes####

clade_colors<-c('#297677','#666666','#78B5BC','#78B5BC','#CE8740')
names(clade_colors)<-c("Carnivore", "Rodent","Other", "Cetacean", "Primate")

pheno_colors <- c("#04516B",'#C54D4D')
names(pheno_colors) <- c("GI","EQ")

monkey_clade_colors<-c('#501537','#9B405F','#EF798A','#827191')
names(monkey_clade_colors)<-c("Great ape", "Old World monkey","New World monkey", "Other")

monkey_clade_colors2<-c('#501537','#9B405F','#EF798A','#827191')
names(monkey_clade_colors2)<-c("Great ape", "Old World monkey","New World monkey", "Other")

#combined clade colors 
combined_clade_colors<-c('#501537','#9B405F','#EF798A','#827191','#297677','#666666','#78B5BC')
names(combined_clade_colors)<-c("Great ape", "OWM","NWM", "Other primate", "Carnivore", "Rodent", "Other")


#clade colors with Human and Mouse having a distinct, individual color (DNase-seq species)
clade_colors_dnase<-c('#297677','#666666','#78B5BC','#78B5BC','#CE8740',"black")
names(clade_colors_dnase)<-c("Carnivore", "Rodent","Other", "Cetacean", "Primate", "DHS species")


#amino acid colors
AA_cols <- c("-" = "grey95", "A" = "#CCFF00", "E" = "#FF0066", "G" = "darkorange", "I"="#33FF33", "L"="#00FF00", "P"="#FFCC33", "Q"="#FF00CC", "R"="#0000CC", "S"="#FF3300", "T"="#FF6600", "V"="#99FF33", "F"="#33FF99")

#region colors
region_colors <- wes_palette("Zissou1", n=7, type = "continuous")
names(region_colors)<-c("upstream1","upstream2","upstream3","exon1","intron","exon2","downstream")






#PANEL 1: protein-evolution ####

#Figure A: brain mass and GI of mammalian species with Trnp1 coding sequence ####
# plot 31 phenotypes 
pheno_data_31species<-readRDS("protein/coevol/pheno_data/pheno_data_31species.rds") %>%
  mutate(`Brain size`=brain_mass/1000) #turn it to kgs

tree.exon1.coding.31sp<-read.tree("protein/trees/tree_TRNP1_coding_31sp.txt") 
length(tree.exon1.coding.31sp$tip.label)

pheno_tree<-ggtree(tree.exon1.coding.31sp, size=0.2)+theme(plot.margin = unit(c(0, 0, 0,0), "cm"))

#order the pheno data in the same order as the tree tips
pheno_tree_ord<-pheno_tree$data[which(pheno_tree$data$isTip),]
desiredOrder  <- rev(pheno_tree_ord[order(pheno_tree_ord$y),]$label)

ggtree(tree.exon1.coding.31sp)+theme(plot.margin = unit(c(0, 0, 0,0), "cm")) +geom_tiplab()+xlim(c(0,150))

pheno_data_31species$order<-match(pheno_data_31species$species, desiredOrder)
pheno_data_31sp_longtable<-pheno_data_31species %>% gather(phenotype,value,`Brain size`,GI)
pheno_data_31sp_longtable$species<-gsub("_"," ", pheno_data_31sp_longtable$species)
pheno_data_31sp_longtable$value2<-pheno_data_31sp_longtable$value


source("figures/new_aes.R")

both_heat<- ggplot(pheno_data_31sp_longtable[pheno_data_31sp_longtable$phenotype=="Brain size",], aes(x=phenotype, y=reorder(species,-order))) +
  geom_tile(aes(fill  =value)) +
  #rename legend to brain mass so that the order of legends is correct
  scale_fill_gradient(name="Brain mass", low = "grey90", high = "#C54D4D", guide = guide_colorbar(barwidth = 1.8, barheight = 0.5)) +
  new_scale("fill") +
  geom_tile(aes(fill = value), data = subset(pheno_data_31sp_longtable, phenotype=="GI")) +
  scale_fill_gradient("GI", low = "grey90", high = "#04516B",guide = guide_colorbar(barwidth = 1.7,barheight = 0.5))+
  scale_y_discrete(position = "right")+
  theme_classic()+
  theme(axis.title= element_blank(), 
        legend.position = "bottom", 
        legend.spacing.x = unit(1.5,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size=5,margin = margin(t = -3)),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-15,0,-10,0),
        axis.text.y.right = element_text(size=6,face="italic",margin= margin(l=-2)),
        axis.ticks=element_blank(), 
        axis.line = element_blank(),
        axis.text.x = element_text(size=7),
        plot.margin = unit(c(0.5, 0.5, 0, -0.2), "cm"))+
  coord_cartesian(xlim=c(1,2), ylim=c(1,31), clip="off") +
  annotate("segment", x = 5.2, xend = 5.2, y = 31.5, yend = 13.7, size=0.2)+
  annotate("text",x=5.4,  y=23,size=2.6, label="Primates", angle=270) +
  annotate("segment", x = 5.2, xend = 5.2, y = 13.3, yend = 10.7, size=0.2)+
  annotate("text",x=5.4,  y=13.5, size=2.6, label="Rodents", angle=270)+
  annotate("segment", x = 5.2, xend = 5.2, y = 10.3, yend = 6.7, size=0.2)+
  annotate("text",x=5.4,  y=7, size=2.6, label="Carnivores", angle=270)



pheno_legend<-get_legend(both_heat)
both_heat_noLeg<-both_heat + theme(legend.position = "none")

fig1A<-ggdraw(plot_grid(plot_grid(pheno_tree, both_heat_noLeg, ncol=2, align='h',axis="tb",scale=c(1,0.97), rel_widths = c(0.25,0.3)),
                        plot_grid(NULL, pheno_legend, rel_widths = c(0.03,0.3)),
                        rel_widths=c(1, 0.1), rel_heights=c(1,0.1), ncol=1))

ggsave("figures/out/panel1/fig1A.pdf",fig1A,width=10.8, height=8, units = 'cm')




#Figures B and C: Coevol results ####
#plot Coevol 30 species 
coevol_estimates<-readRDS("protein/coevol/results/for_figures/coevol_3phenos_31sp_summarized.rds")
coevol_estimates$clade<-factor(coevol_estimates$clade, levels=c("Primate", "Rodent", "Carnivore", "Other"))

#plot both axis in log10
options(scipen=10000) #turning off scientific notation
fig1B<-ggplot(coevol_estimates, aes(x=omega, y=signif(brain/1000,1)))+
  geom_smooth(method="lm", color="darkgrey", alpha=0.23)+
  xlab("Omega of TRNP1")+
  ylab("Brain size (kg)")+
  geom_errorbarh(aes(xmin=CI.low, xmax=CI.high), color="grey70",size=0.3)+
  theme_classic()+
  geom_point(aes(color=clade),size=1.2)+
  theme(legend.position = "none",
        axis.title.y=element_text( size=8.5), 
        axis.text.y=element_text(size=7, margin=margin(l=-4)),
        axis.title.x=element_text(size=8.5),
        axis.text.x=element_text(size=7),
        axis.line = element_line(size=0.3))+
  geom_text_repel(aes(label=interesting_species),  size=2.5, point.padding=0.05) +
  annotate("text",x=0.07,  y=7,size=2.5, label="0.64") +
  scale_color_manual(values=clade_colors)+
  scale_y_log10(breaks=c(0.01,1,10))+
  scale_x_log10(breaks=c(0.1,0.3,0.9))

ggsave(file="figures/out/panel1/fig1B.pdf",fig1B, width=6.5, height=6.5, units="cm")



fig1C_unfinished<-ggplot(coevol_estimates, aes(x=omega, y=GI))+
  geom_smooth(method="lm", color="darkgrey", alpha=0.23)+
  xlab("Omega of TRNP1")+
  ylab("Gyrification (GI)")+
  geom_errorbarh(aes(xmin=CI.low, xmax=CI.high), color="grey70", size=0.3)+
  theme_classic()+
  geom_point(aes(color=clade), size=1.2)+
  theme(legend.position = "top",
        axis.title.y=element_text( size=8.5), 
        axis.text.y=element_text(size=7, margin=margin(l=-2)),
        axis.title.x=element_text(size=8.5),
        axis.text.x=element_text(size=7),
        axis.line = element_line(size=0.3),
        legend.spacing.x = unit(1,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,0,-7,0),
        legend.box.margin=margin(0,0,-7,0))+
  geom_text_repel(aes(label=interesting_species), size=2.5, point.padding = 0.1)+
  annotate("text",x=0.07,  y=4.5,size=2.5, label="0.69") +
  scale_color_manual(values=clade_colors)+
  scale_y_log10()+
  scale_x_log10(breaks=c(0.1,0.3,0.9))


fig1C<-fig1C_unfinished+theme(legend.position = "none")
fig.legend.1 <- cowplot::get_legend(fig1C_unfinished)

ggsave(file="figures/out/panel1/fig1C.pdf", fig1C, width=6.5, height=6.5, units="cm")




#Figure D: Trnp1 structure ####

#sites under selection - 35,41,62,63,86,193 (relative to the 45-species aligment. it is 248 positions long in total)
full_prot_aligned_45sp<-readAAStringSet("protein/fastas/prank_output/prank_TRNP1_coding_45species.best.pep.fas") 
Mouse_seq_45<-AAString(as.character(full_prot_aligned_45sp[names(full_prot_aligned_45sp) %in% "Mus_musculus"]))
fun1 <- function(x) stri_length(x) - stri_count_fixed(x, "-")

#I have the IDR annotations from Miriam Esgleas on the mouse sequence and the sites under positive selection relative to the alignment --> need to identify these positions in the mouse sequence to accurately show the selected sites vs IDRs on the mouse sequence schematic
s_45 <- strsplit(x=as.character(full_prot_aligned_45sp[names(full_prot_aligned_45sp) %in% "Mus_musculus"]), split=character(0))
mouse_seq_df<-data.frame(seq_45=s_45, pos_alignment=1:length(s_45[[1]]))

for (i in 1:nrow(mouse_seq_df)){
  mouse_seq_df$pos_mouse_seq[i]<-fun1(Mouse_seq_45[1:i])
}

sites_under_sel_pp95<-readRDS("protein/PAML/PAML_M8_NEB_sites_pp95.rds")
mouse_seq_df$sel<-NA
for (i in sites_under_sel_pp95$`Alignment position`){
  mouse_seq_df$sel[mouse_seq_df$pos_alignment==i]<-mouse_seq_df$pos_mouse_seq[i]
}
#extract mouse-sites under positive selection
sel<-mouse_seq_df$sel[!is.na(mouse_seq_df$sel)]

#murine seq structure (Miriam's paper): IDR1: 1-103, (alpha-helix 87-144), IDR2: 165-178 & 196-223
#calculate the distances between region borders
df<-tibble(x=c(0,0,0,0,0), 
           y=c(103, 165-103, 178-165, 196-178, 223-196))

fig1D_R<-ggplot(df, aes(x=x, y=y))+
  geom_bar(stat="identity", color="grey30",size=0.2,fill=c("#F0D2D2","grey","#F0D2D2","grey","#F0D2D2"), width=0.2)+ 
  coord_flip() +
  theme_void()+
  annotate("text", y=-10, x=0, label="N", size=4, color="grey30")+
  annotate("text", y=223+10, x=0, label="C", size=4, color="grey30")+
  xlim(c(-0.3,0.3))+
  annotate("segment", y=sel, yend=sel, x=-0.1, xend=0.1, color="#04516B")+
ylim(c(-30,250))

fig1D<-plot_grid(NULL,fig1D_R,NULL, ncol=1, rel_heights = c(0.3,0.2,0.3))
ggsave("figures/out/panel1/fig1D_R.pdf",fig1D_R,height=1, width=5)

fig1D_upgraded<- ggdraw() + 
  cowplot::draw_image(magick::image_read("figures/out/panel1/fig1D_upgraded2.png"))





#Figure E: Overview of the experimental setup ####
fig1E<- ggdraw() + 
  cowplot::draw_image(magick::image_read("figures/out/panel1/mouse_exp_small2.png"))


#Figure F: The effect of Trnp1 orthologues on proliferation in NSCs ####
combined_glm_all<-readRDS("protein/data/proliferation/proliferation_LR_res_orthologues.rds") 
combined_glm_all$clade <-stringr::str_to_title(combined_glm_all$clade)

fig1F<-ggplot(combined_glm_all, aes(x=term, y=prolif_prob, color=clade))+
  geom_point(size=2.7)+
  geom_errorbar(aes(ymin=prolif_prob-prolif_stderr, ymax=prolif_prob+prolif_stderr), width=0.2)+
  theme_bw()+
  ylab("Proliferation rate")+
  scale_color_manual(values=clade_colors)+
  #human vs mouse
  geom_segment(aes(x=4, y=0.65, xend=5, yend=0.65), color="grey40", size=0.2)+
  geom_segment(aes(x=4, y=0.653, xend=4, yend=0.647), color="grey40", size=0.15)+
  geom_segment(aes(x=5, y=0.653, xend=5, yend=0.647), color="grey40", size=0.15)+
  annotate("text", x=4.5, y=0.67, label=".", size=3.5, color="grey30")+
  #human vs dolphin
  geom_segment(aes(x=6, y=0.72, xend=5, yend=0.72), color="grey40", size=0.2)+
  geom_segment(aes(x=6, y=0.723, xend=6, yend=0.717), color="grey40", size=0.15)+
  geom_segment(aes(x=5, y=0.723, xend=5, yend=0.717), color="grey40", size=0.15)+
  annotate("text", x=5.5, y=0.725, label="*", size=3.5, color="grey30")+
  #human vs macaque
  geom_segment(aes(x=1, y=0.69, xend=5, yend=0.69), color="grey40", size=0.2)+
  geom_segment(aes(x=1, y=0.693, xend=1, yend=0.687), color="grey40", size=0.15)+
  geom_segment(aes(x=5, y=0.693, xend=5, yend=0.687), color="grey40", size=0.15)+
  annotate("text", x=3, y=0.695, label="**", size=3.5, color="grey30")+
  #human vs galago
  geom_segment(aes(x=2, y=0.67, xend=5, yend=0.67), color="grey40", size=0.2)+
  geom_segment(aes(x=2, y=0.673, xend=2, yend=0.667), color="grey40", size=0.15)+
  geom_segment(aes(x=5, y=0.673, xend=5, yend=0.667), color="grey40", size=0.15)+
  annotate("text", x=3.5, y=0.675, label="*", size=3.5, color="grey30")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text( size=8.5), 
        axis.text.y=element_text(size=7, margin=margin(l=-3), color="grey30"),
        axis.text.x=element_text(size=8,color="black", angle=45, hjust=1, vjust=1),
        #axis.line = element_line(size=0.01),
        legend.position = "none")

ggsave("figures/out/panel1/fig1F.pdf", fig1F, width=3.5, height=3)



#Figure G: Prediction of GI using proliferation rates ####
#vs GI
prolif_vs_GI<-readRDS("protein/data/proliferation/prolif_vs_GI.rds")


fig1G<-ggplot(prolif_vs_GI, aes(x=prolif_prob, y=GI))+
  geom_smooth(method="lm", color="grey80", alpha=0.13)+
  geom_errorbarh(aes(xmin=prolif_stderr_min, xmax=prolif_stderr_max),color="darkgrey", height=0.1)+
  geom_point(aes(color=clade),size=2.7)+
  scale_color_manual(values=clade_colors)+
  theme_bw()+
  theme(axis.title = element_text(size=8.5),
        #axis.text = element_text(size=8, color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y=element_text(size=7, margin=margin(l=-3), color="grey30"),
        axis.text.x=element_text(size=7,color="grey30", margin=margin(b=1.5, t=1.5)))+
        #axis.line = element_line(size=0.05))+
  xlab("Proliferation rate")+
  ylab("Gyrification (GI)")+
  scale_x_log10()+
  scale_y_log10()+
  annotate("text", label='0.68',x=0.48,y=8, size=2.5)

ggsave("figures/out/panel1/fig1G.pdf", fig1G, width=5, height=4)




#combining all ####
plot4<-cowplot::plot_grid(fig1D_upgraded,fig1E,  ncol = 2, rel_widths = c(0.3,0.28), labels = c("b","e"), label_size = 11, align="h", axis="t", scale =c(0.9,1))

plot1<-cowplot::plot_grid(fig1A, plot4,  ncol = 1, rel_heights = c(0.8,0.4), labels = c("a",NULL), label_size = 11, label_y=0.95)

plot2_preliminary<-cowplot::plot_grid(fig1B, fig1C, ncol=2, labels = c("c","d"), rel_widths = c(0.32,0.3), label_size = 11)

plot2<-cowplot::plot_grid(fig.legend.1,plot2_preliminary, ncol=1, nrow=2, rel_heights = c(0.1,0.4), rel_widths = c(0.3,0.8))

plot3<-cowplot::plot_grid( fig1F, fig1G, ncol = 2, labels = c("f", "g"), rel_widths = c(0.4,0.4), label_size = 11, align = "h", axis="bt")

lower_panel<-plot_grid(plot3, ncol = 2, rel_widths = c(0.3,0.03, 0.85), labels=c("e",NULL, NULL), label_size = 11)

lower_and_upper<-cowplot::plot_grid(plot2,NULL,plot3, ncol=1, rel_heights = c(0.4,0.01,0.32))

final_plot<-cowplot::plot_grid(plot1,NULL,lower_and_upper, ncol=3, rel_widths = c(0.5,0.01,0.49)) #ignores rel_heights

ggsave("figures/out/panel1/final_panel1_2.pdf", final_plot, width=183, height=115, units = "mm")
ggsave("figures/out/panel1/final_panel1.png", final_plot, width=183, height=115, units = "mm")








#..................................................................#####
# Supplementary panel 1 ####

#Extended Fig A: Alignment of Trnp1 coding-sequence####
tree.exon1.coding.45sp<-read.tree("protein/trees/tree_TRNP1_coding_45sp.txt") 
basic.tree<-ggtree(tree.exon1.coding.45sp,size=0.2)+theme(plot.margin = unit(c(0, 0, 0,0), "cm"))

ext1A<-msaplot(basic.tree, "protein/fastas/prank_output/prank_TRNP1_coding_45species.best.pep.fas", width=7, bg_line=FALSE) +theme(legend.position = "none")

ggsave("figures/out/extended1/ext1A.pdf", ext1A, width=183, height=110,units = "mm")





#Extended Fig B: Positively selected sites####
#plot sites under selection (PAML)
full_prot_aligned_45sp<-readAAStringSet("protein/fastas/prank_output/prank_TRNP1_coding_45species.best.pep.fas") 
sites_under_sel_pp95<-readRDS("protein/PAML/PAML_M8_NEB_sites_pp95.rds")


sites_under_sel_45sp<-data.frame(species=names(full_prot_aligned_45sp))
for (i in sites_under_sel_pp95$`Alignment position`){
  stringz<-data.frame(species=names(full_prot_aligned_45sp),substr(full_prot_aligned_45sp, i, i))
  sites_under_sel_45sp<-inner_join(sites_under_sel_45sp, stringz, by="species")
}
colnames(sites_under_sel_45sp)[2:7]<-paste0("s",sites_under_sel_pp95$`Alignment position`)
sites_under_sel_45sp$s193[sites_under_sel_45sp$species=="Macaca_fascicularis"]<-NA

basic.tree_ord<-basic.tree$data[which(basic.tree$data$isTip),]
desiredOrder_sites  <- rev(basic.tree_ord[order(basic.tree_ord$y),]$label)
sites_under_sel_45sp$order<-match(sites_under_sel_45sp$species, desiredOrder_sites)
sites_under_sel_45sp_long<-tidyr::gather(sites_under_sel_45sp, key="site", value="AA",2:7) %>%
  mutate(AA=as.factor(AA),
         site=factor(site, levels=c("s35", "s41",  "s62", "s63", "s86","s193")),
         species=gsub("_"," ", species))


sites_plot<-ggplot(sites_under_sel_45sp_long, 
                   aes(x=site, y=reorder(species,-order )))+
  geom_tile(aes(fill=AA), color="white", size=0.3, alpha=0.7)+
  geom_text(aes(label=AA), color="grey20",size=1.6)+
  theme_light()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(size=6, color="grey30"), 
        legend.position = "none",
        axis.text.y=element_text(size=6, face="italic", family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_discrete(position="right")+
  scale_fill_manual(values=AA_cols)+ 
  theme(plot.margin = unit(c(0, 0, 0, -0.3), "cm"))                                           

ext1B<-plot_grid(basic.tree, sites_plot, rel_widths = c(0.15,0.3), align = "h", scale=c(1,0.94))
ggsave("figures/out/extended1/ext1B.pdf", ext1B, width=183, height=100,units = "mm")



#AA counts for 1D####
aa_counts<-sites_under_sel_45sp_long %>%
  group_by(site, AA) %>%
  summarize(n=length(AA)) %>%
  #exclude gap sign cause it looks weird
  filter(AA!="-") %>%
  arrange(site, -n) %>%
  group_by(site)%>% 
  mutate(rank=row_number()) %>%
  ungroup() %>%
  mutate(site2=gsub("s","",site),
         site2=as.numeric(site2))


fig1D_R2<-ggplot(df, aes(x=x, y=y))+
  geom_bar(stat="identity", color="grey30",size=0.2,fill=c("#F9D8A8","grey85","#F9D8A8","grey85","#F9D8A8"), width=0.2, alpha=0.9)+ 
  coord_flip() +
  theme_void()+
  annotate("text", y=-10, x=0, label="N", size=4, color="black")+
  annotate("text", y=223+10, x=0, label="C", size=4, color="black")+
  xlim(c(-0.2,0.2))+
  annotate("segment", y=sel, yend=sel, x=-0.1, xend=0.1, color="#BF4342")+
  ylim(c(-30,250))


#separate out the last site because it is located further away
fig1D_AAs1<-ggplot(aa_counts %>% filter(site2!=193), aes(x=site, y=-rank, size=n, color=AA))+
  geom_text(aes(label=AA), fontface="bold")+
  theme_void()+
  scale_color_manual(values=AA_cols)+
  scale_size_continuous(range=c(3,7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(-1, 0,0,0), "cm"))+
  ylim(c(-6,0.9))


fig1D_AAs2<-ggplot(aa_counts %>% filter(site2==193), aes(x=site, y=-rank, size=n, color=AA))+
  geom_text(aes(label=AA), fontface="bold")+
  theme_void()+
  scale_color_manual(values=AA_cols)+
  scale_size_continuous(range=c(3,7)) +
  theme(legend.position = "none",
        plot.margin = unit(c(-1, 0,0,0), "cm"))+
  ylim(c(-6,0.9))

plot_grid(fig1D_R2,
          NULL,
          plot_grid(NULL,fig1D_AAs1, NULL, fig1D_AAs2,NULL, ncol=5, rel_widths = c(0.07,0.4,0.05,0.1,0.1)),
          ncol=1, rel_heights = c(0.17,0.05,0.5))

ggsave("figures/out/panel1/fig1D_R_sites.pdf",height=6, width=10, units = "cm")





# Extended Fig C: Coevol estimates body mass vs GI ####
coevol_estimates$BodyM_kg<-coevol_estimates$BodyM*0.001

ext1C<-ggplot(coevol_estimates, aes(x=omega, y=BodyM_kg))+
  geom_smooth(method="lm", color="darkgrey", alpha=0.23)+
  xlab("omega of TRNP1")+
  geom_errorbarh(aes(xmin=CI.low, xmax=CI.high), color="grey70", size=0.3)+
  theme_classic()+
  geom_point(aes(color=clade))+
  theme(legend.position = "none",
        axis.title.y=element_text( size=8.5), 
        axis.text.y=element_text(size=7), # margin=margin(l=-4)
        axis.title.x=element_text(size=8.5),
        axis.text.x=element_text(size=7),
        axis.line = element_line(size=0.3))+
  geom_text_repel(aes(label=interesting_species), size=2.5)+
  scale_color_manual(values=clade_colors)+
  scale_y_log10(breaks=c(1,10,100,1000))+
  scale_x_log10(breaks=c(0.1,0.3,0.9))+
  ylab("Body mass (kg)")

ggsave("figures/out/extended1/ext1C.pdf", ext1C, width=6, height=6, units="cm")







#Extended Fig D: Effect of Trnp1 on proliferation ####
combined_glm_Trnp1<-readRDS("protein/data/proliferation/proliferation_LR_res_TRNP1.rds")

combined_glm_Trnp1$term<-gsub("no","control", combined_glm_Trnp1$term)
combined_glm_Trnp1$term<-gsub("yes","TRNP1", combined_glm_Trnp1$term)

ext1D<-ggplot(combined_glm_Trnp1, aes(x=term, y=prolif_prob, color=term))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=prolif_stderr_min, ymax=prolif_stderr_max), width=0.1)+
  theme_bw()+scale_color_manual(values=c("grey60","black"))+
  ylab("Proliferation rate")+
  geom_segment(aes(x=1, y=0.58, xend=2, yend=0.58), color="grey40", size=0.2)+
  geom_segment(aes(x=1, y=0.583, xend=1, yend=0.577), color="grey40", size=0.15)+
  geom_segment(aes(x=2, y=0.583, xend=2, yend=0.577), color="grey40", size=0.15)+
  annotate("text", x=1.5, y=0.59, label="****", size=3.5, color="grey30")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_text( size=8.5), 
        axis.text.y=element_text(size=7, margin=margin(l=-3), color="grey30"),
        axis.text.x=element_text(size=7,color="black"),
        #axis.line = element_line(size=0.01),
        legend.position = "none")

ggsave("figures/out/extended1/ext1D.pdf", ext1D, width=2.5, height=2.5)





# Extended Figs E & G: Coevol correlation matrices ####

#combining all####
lower_right<-plot_grid(ext1C, ext1D, ncol=1, rel_heights = c(0.65,0.5),
                       align="v",axis = "l", labels=c("c","d"), label_size = 11)

lower_panel_ext1<-plot_grid(ext1B, lower_right, NULL, rel_widths = c(0.52,0.21,0.03), ncol=3)
ext1<-plot_grid(ext1A,lower_panel_ext1, ncol=1, rel_heights = c(0.51,0.63), labels=c("a","b"), label_size = 11)

ggsave("figures/out/extended1/ext_fig1.pdf", ext1, scale=0.9, width=183, height=179,units = "mm")





#..................................................................#####

# PANEL 2: MPRA activity ####

# Figure A: MPRA overview ####
fig2A<- ggdraw() + 
  cowplot::draw_image(magick::image_read("figures/out/panel2/MPRA_overview_larger.png"))


# Figure B: Summarized activity per region per species ####
source("figures/new_aes.R")
activity_overlap_summary<-readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds")

#plot only the ones with GI and brain weight present +human1 activity --> 45 species btw
pheno_data<-readRDS("pheno_data/pheno_data.rds")


#data
activity_summary_hum1_withPheno<- activity_overlap_summary %>% 
  filter(cell_line=="human1") %>%
  left_join(pheno_data[,c("brain_mass","GI","species")],
            by=c("species")) %>%
  filter(!is.na(GI) & !is.na(brain_mass))

activity_summary_hum1_withPheno$region<-factor(activity_summary_hum1_withPheno$region, levels=c("upstream1", "upstream2", "upstream3","exon1", "intron", "exon2", "downstream"))


#tree 
mammaltree<-read.tree("protein/trees/mammaltree.txt") 
mpra_tree_full<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% activity_summary_hum1_withPheno$species])

mpra_tree_plot<-ggtree(mpra_tree_full, size=0.2)+theme(plot.margin = unit(c(0, 0, 0,0), "cm"))
mpra_tree_plot+geom_tiplab()+xlim(0,200)

#mark primates and OWMs+great apes
mpra_tree_plot2<-ggtree(mpra_tree_full, size=0.2)+
  ggimage::theme_transparent()+
  theme(plot.margin = unit(c(0, -0.1, 0,0), "cm"))+
  geom_hilight(node=46, fill="white")+
  geom_hilight(node=54, fill="grey85")+
  geom_hilight(node=57, fill="grey65")
length(mpra_tree_full$tip.label) #45

#order the data in the same order as the tree tips
mpra_tree_ord<-mpra_tree_plot$data[which(mpra_tree_plot$data$isTip),]
desiredOrder_mpra  <- rev(mpra_tree_ord[order(mpra_tree_ord$y),]$label)

activity_summary_hum1_withPheno$order<-match(activity_summary_hum1_withPheno$species, desiredOrder_mpra)
activity_summary_hum1_withPheno$species2<-gsub("_"," ",activity_summary_hum1_withPheno$species)

activity_summary_hum1_withPheno_long<-activity_summary_hum1_withPheno %>% 
  mutate(brain_mass=log2(brain_mass/100), GI_zzz=log2(GI)) %>%
  gather(phenotype,value, brain_mass,GI)

xval<-6.1
xtextval<-xval+0.3

both_heat_MPRA2<- ggplot(activity_summary_hum1_withPheno_long[activity_summary_hum1_withPheno_long$phenotype=="brain_mass",], aes(x=phenotype, y=reorder(species2,-order))) +
  geom_tile(aes(fill  =value), alpha=0.25) +
  scale_fill_gradient("brain_mass", low = "grey90", high = "#C54D4D", 
                      guide = guide_colorbar(barwidth = 1.1, barheight = 0.4)) +
  scale_x_discrete(label=c("Brain\nsize","GI"))+
  new_scale("fill") +
  geom_tile(aes(fill = value),alpha=0.25, data = subset(activity_summary_hum1_withPheno_long, phenotype=="GI")) +
  scale_fill_gradient("GI_zzz", low = "grey90", high = "#04516B",
                      guide = guide_colorbar(barwidth = 1.1,barheight = 0.4))+
  scale_y_discrete(position = "right")+
  theme_classic()+
  theme(axis.title= element_blank(), 
        legend.position = "top",
        legend.justification = "left",
        legend.spacing.x = unit(1,"mm"),
        legend.title = element_blank(),
        legend.text = element_text(size=5,margin = margin(t = -3, b=-2)),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-4,0,-4,0),
        legend.box.margin=margin(-5,0,-5,0),
        axis.text.y.right = element_text(size=6,face="italic",margin= margin(l=-2)),
        axis.ticks=element_blank(), 
        axis.line = element_blank(),
        axis.text.x = element_text(size=6.5, color="black", margin = margin(t=-1)),
        plot.margin = unit(c(0, 0.3, 0, -0.1), "cm"))+
  coord_cartesian(xlim=c(1,2), ylim=c(1,45), clip="off") +
  annotate("segment", x = xval, xend = xval, y = 45.3, yend = 21.7, size=0.2)+
  annotate("text",x=xtextval,  y=33,size=2.6, label="Primates", angle=270) +
  annotate("segment", x = xval, xend = xval, y = 20.3, yend = 18.7, size=0.2)+
  annotate("text",x=xtextval,  y=20.7, size=2.6, label="Rodents", angle=270)+
  annotate("segment", x = xval, xend = xval, y = 17.3, yend = 14.7, size=0.2)+
  annotate("text",x=xtextval,  y=14.2, size=2.6, label="Carnivores", angle=270)




#with region names below ####
total_activity_hum1_names<-ggplot(activity_summary_hum1_withPheno, aes(x=region, y=reorder(species2,-order)))+
  geom_tile(aes(fill=log2_total_activity))+ 
  scale_x_discrete(expand=c(0,0), labels=c("U1","U2","U3","E1","I","E2","D"))+
  scale_fill_gradient2(low="#F7C1BB",mid="grey70", high="#5F0F40", midpoint = 0)+ 
  theme_bw()+
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x = element_text( size=6.5, color="black", margin = margin(t=-1)),
        axis.ticks = element_blank(),
        legend.position = "top",
        legend.justification = "right",
        legend.spacing.x = unit(1.5,"mm"),
        legend.text = element_text(size=5,margin = margin(t = -3, b=-2)),
        legend.title=element_blank(),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-2,0,-4,0),
        legend.box.margin=margin(-5,0,-5,0),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey80"),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0, 0, -0.1), "cm"))+
  guides(fill = guide_colorbar( barwidth = 2,barheight = 0.4))+
  xlab("log2 (Regulatory activity)")

fig2B_2<- plot_grid(mpra_tree_plot2,total_activity_hum1_names, both_heat_MPRA2, rel_widths = c(0.11,0.18, 0.19), ncol=3, align = "h", axis="tb", scale=c(0.97,0.95,0.95))
ggsave("figures/out/panel2/fig2B.pdf", fig2B_2, height=97, width=99, units="mm")





# Figure C: Motif binding site score variation across OWMs & apes ####

scores_top22TFs_long2<-readRDS("regulation/data/TFs/ClusterBuster/results/scores_22TFs_long.rds")
tree.intron.tfs22.2<-readRDS("regulation/data/TFs/ClusterBuster/results/tree_22TFs.rds")
intron_motifs_1_OWMs<-readRDS("regulation/data/TFs/ClusterBuster/results/intron_motifs_1_OWMs.rds")
intron_22TFs_summarized<-readRDS("regulation/data/TFs/expression/summarized_info_top22TFs_ctrain.rds")



# bottom- plot motif distribution across the human intron ####
intron_motifs_hum<-read.table("regulation/data/TFs/ClusterBuster/results/OWMs_apes_intron_weighted/subset_cbust/motifs/motif_tables/Homo_sapiens.txt")
colnames(intron_motifs_hum)<-c("motif","start","end","strand","motif_score","MB_sequence")

intron_motifs_hum_max<-intron_motifs_hum %>%
  inner_join(intron_motifs_1_OWMs) %>%
  group_by(SYMBOL_clean) %>%
  slice_max(motif_score) %>%
  ungroup() %>%
  arrange(start) %>%
  mutate(plot_level=rep(1:3, length.out=length(start))) 

#as a line - at the middle spot
ggplot(intron_motifs_hum_max)+
  geom_segment(aes(x=189,xend=397, y=0, yend=0),color="grey90", size=5)+
  geom_segment(aes(x=start+(end-start)/2,xend=start+(end-start)/2+0.5, y=0, yend=0),color="black", size=5)+
  ggrepel::geom_text_repel(aes(x=(end-(end-start)/2), y=0, label=SYMBOL_clean),  size=1.5, nudge_y = 1) +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/out/panel2/motifs/ctrain_human_motif_pos_6_relClust_mid_text.pdf", height = 1, width=6)


ggplot(intron_motifs_hum_max)+
  geom_segment(aes(x=189,xend=397, y=0, yend=0),color="grey90", size=5)+
  geom_segment(aes(x=start+(end-start)/2,xend=start+(end-start)/2+0.5, y=0, yend=0),color="black", size=5)+
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/out/panel2/motifs/ctrain_human_motif_pos_6_relClust_mid.pdf", height = 1, width=6)



#middle- plot cbust scores####
source("figures/new_aes.R")
scores_top22TFs_long2_2<- scores_top22TFs_long2 %>% 
  left_join(pheno_data %>% 
              dplyr::select(species, GI)) %>%
  left_join(activity_overlap_summary %>% filter(region=="intron" & cell_line=="human1") %>%
              dplyr::select(species,log2_total_activity))  %>% 
  mutate(log2_GI=log2(GI)) %>%
  #include activity 
  gather(phenotype,value, log2_total_activity,log2_GI) %>%
  #reorder heatmap according to the TF motif positions
  left_join(intron_motifs_hum_max[,c("motif","start","end")]) %>%
  mutate(fontface=case_when(SYMBOL_clean %in% c("CTCF","SOX8", "ZBTB26") ~ "bold",
                            T ~"plain"))


p <- ggtree(tree.intron.tfs22.2, size=0.3) + 
  geom_tiplab(size=2)+xlim(0,120)+
  ggimage::theme_transparent()+
  theme(plot.margin=margin(0,-15,0,0))

p_scores <- ggplot(scores_top22TFs_long2_2 %>% 
                     mutate(SYMBOL_clean=gsub("_VAR2","",SYMBOL_clean)), 
                   aes(x=reorder(SYMBOL_clean,end), y=species2)) + 
  geom_tile(aes(fill=stand_score), color="white",alpha=0.8) + 
  theme_tree2()+
  scale_fill_gradient2(high="black",low="white", mid="grey50", 
                       #high="#8B786D",low="#EBF5EE", mid="#BFA89E", 
                       midpoint = 0.5,
                       guide = guide_colorbar(barwidth = 2.3,barheight = 0.5, title.position = "top"), 
                       name="standardized\nbinding score")+
  theme(legend.title = element_text(size=6),
        legend.text = element_text(size=5),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,0,0,0),
        plot.margin=margin(0,0,0,-15),
        legend.position = 'bottom',
        panel.background = element_rect(fill="transparent",colour=NA),
        legend.box.margin=margin(0,-3,0,-15),
        axis.text.x=element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title.align = 0.5)+ggimage::theme_transparent()
p_scores_legend<-get_legend(p_scores)
p_scores<-p_scores +theme(legend.position = "none")


both_heat_scores<- ggplot(scores_top22TFs_long2_2[scores_top22TFs_long2_2$phenotype=="log2_total_activity",], aes(x=phenotype, y=reorder(species2,-order))) +
  geom_tile(aes(fill  =value), color="white") +
  scale_fill_gradient("log2(intron\nactivity)", low = "grey90", high = "#FCC25F", 
                      guide = guide_colorbar(barwidth = 1.8, barheight = 0.5, title.position = "top")) +
  new_scale("fill") +
  geom_tile(aes(fill = value), color="white", data = subset(scores_top22TFs_long2_2,
                                                            phenotype=="log2_GI")) +
  scale_fill_gradient("log2(GI)", low = "grey90", high = "#04516B",
                      guide = guide_colorbar(barwidth = 1.8,barheight = 0.5, title.position = "top"))+
  theme_tree2()+
  scale_x_discrete(position="top", labels=c("GI","intron"))+
  theme(legend.title = element_text(size=6),
        legend.text = element_text(size=5),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,-3,0,-15),
        legend.position = "bottom",
        axis.text.x.top = element_text(angle=-90, hjust=1, vjust=0.2, size=7, margin = margin(b=-2),
                                       color="black"),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.title.align = 0.5)+ggimage::theme_transparent()

both_heat_scores_legend<-get_legend(both_heat_scores)
both_heat_scores<-both_heat_scores +theme(legend.position = "none")

library(aplot)
fig<-p_scores %>% insert_left(p, width=0.36) %>% insert_right(both_heat_scores, width=.12)
ggsave("figures/out/panel2/motifs/cbust_scores_heatmap_angle90.pdf",fig, height=7, width=12, units="cm")



#top- plot GO terms####
score_top<-ggplot(scores_top22TFs_long2_2 %>% 
                    left_join(intron_22TFs_summarized[,c("SYMBOL_clean","GO_term")]) %>%
                    mutate(SYMBOL_clean=gsub("_VAR2","",SYMBOL_clean)) %>% distinct(SYMBOL_clean, end,GO_term), 
                  aes(x=reorder(SYMBOL_clean,end), y=1)) + 
  geom_tile(aes(fill=GO_term), alpha=0.4)+
  geom_text(aes(label=SYMBOL_clean, color=GO_term), angle=-90,hjust=1, vjust=0.5, 
            size=4.5, color="black", nudge_y = -0.45)+
  theme_void()+
  scale_fill_manual(name="GO term",
                    values=c("#84a98c","#52796f","#DFE4DC"),
                    limits=c("cell proliferation","hormone levels and\ncell proliferation"))+
  theme(legend.title = element_text(size=6, color="black"),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,-2,0,-2),
        legend.key.size = unit(0.4, "cm"),
        legend.position = "bottom",
        legend.title.align = 0.5, legend.direction = "vertical")+
  guides(fill = guide_legend(title.position = "top"))


score_top_legend<-get_legend(score_top)
score_top<-score_top +theme(legend.position = "none")
ggsave("figures/out/panel2/motifs/GO_annot_top.pdf",score_top, height=2.2, width=13, units="cm")

plot_grid(score_top_legend,p_scores_legend, both_heat_scores_legend, ncol=3, align = "v", axis="tb")
ggsave("figures/out/panel2/motifs/legends.pdf", height=2, width=8, units="cm")






# Figure D: PGLS of intron activity vs GI ####
activity_summary_hum1_withPheno_intron<-activity_summary_hum1_withPheno %>% 
  filter(region=="intron")

intron_tree_full<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% activity_summary_hum1_withPheno_intron$species])
ggtree(intron_tree_full)+geom_tiplab()+xlim(0,150) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)
length(intron_tree_full$tip.label) #37


# fig D1: pgls across all species####
rownames(activity_summary_hum1_withPheno_intron)<-activity_summary_hum1_withPheno_intron$species
name.check(intron_tree_full,activity_summary_hum1_withPheno_intron)

mod.totalA.GI.intron<-gls(log2(GI)~log2_total_activity, data=activity_summary_hum1_withPheno_intron, correlation=corBrownian(value=1,phy=intron_tree_full, form=~species), method="ML")
mod.totalA.GI.zero<-gls(log2(GI)~1, data=activity_summary_hum1_withPheno_intron, correlation=corBrownian(value=1,phy=intron_tree_full, form=~species), method="ML")

#diagnostics
qqnorm(mod.totalA.GI.intron$residuals, pch = 1, frame = FALSE)
qqline(mod.totalA.GI.intron$residuals, col = "steelblue", lwd = 2)
plot(mod.totalA.GI.intron$fitted, mod.totalA.GI.intron$residuals)


#plot activity vs GI 
#join with the table
activity_summary_hum1_withPheno_intron<-activity_summary_hum1_withPheno_intron %>%
  left_join(as.data.frame(mod.totalA.GI.intron$fitted) %>%
              rownames_to_column("species") %>%
              dplyr::rename(fitted_log2_GI_all="mod.totalA.GI.intron$fitted"))

activity_summary_hum1_withPheno_intron$clade<-factor(activity_summary_hum1_withPheno_intron$clade, levels=c("Primate","Rodent","Carnivore", "Other"))

fig2C_1<-ggplot(activity_summary_hum1_withPheno_intron, aes(x=log2_total_activity, y=log2(GI)))+
  geom_point(aes(color=clade), size=1.5, alpha=0.7)+
  theme_bw()+
  xlab("log2 (Intron activity)")+
  geom_line(aes(y=fitted_log2_GI_all), color="grey50")+  
  scale_color_manual(values=clade_colors)+  
  scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  #scale_y_continuous(limits = c(-0.01, 2.7))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position = "right",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=8),
        axis.text=element_text(size=7),
        axis.text.y=element_text(margin=margin(l=-2)),
        axis.text.x=element_text(margin=margin(b=-1.5)),
        #legend.spacing.x = unit(0,"mm"),
        legend.spacing.y = unit(0, "mm"),
        legend.text = element_text(size=6,margin = margin(l = -5, r=-3, t=-7, b=-7)),
        legend.title=element_blank(),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-10,-10,-10,-10),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.key = element_rect(fill = NA),
        legend.key.width = unit(0.35, "cm"),
        legend.key.height = unit(0.35, "cm"),
        panel.border = element_rect(fill = NA, colour = "grey60"),
        plot.margin = unit(c(0.15,0,0,0), "cm"),
        plot.title = element_text(family = "Helvetica", face = "plain", size = (6), hjust=0.5, margin = margin()))+
  ylab("log2 (Gyrification)")+
  ggtitle("All mammals")+
  annotate("text", label=paste0(round(rr2::R2.resid(mod.totalA.GI.intron,mod.totalA.GI.zero), digits = 2)),  -Inf, Inf,  hjust = -0.3, vjust = 1.7, size=2.5)

figC1_legend<-get_legend(fig2C_1)
fig2C_1_final<-fig2C_1 + ggplot2::theme(legend.position = "none")




# fig D3: only OWMs and apes ####
intron_tree_apes_OWMs<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% activity_summary_hum1_withPheno_intron$species[activity_summary_hum1_withPheno_intron$primate_clade %in% c("Great ape","Old World monkey")]])

rownames(activity_summary_hum1_withPheno_intron)<-activity_summary_hum1_withPheno_intron$species

mod.totalA.GI.intron.apesOWMs<-gls(log2(GI)~log2_total_activity, 
                                   data=activity_summary_hum1_withPheno_intron,
                                   correlation=corBrownian(value=1,phy=intron_tree_apes_OWMs, form=~species),
                                   subset=c(primate_clade %in% c("Great ape","Old World monkey")), 
                                   method="ML")
mod.totalA.GI.null.apesOWMs<-gls(log2(GI)~1, 
                                 data=activity_summary_hum1_withPheno_intron,
                                 correlation=corBrownian(value=1,phy=intron_tree_apes_OWMs, form=~species),
                                 subset=c(primate_clade %in% c("Great ape","Old World monkey")), 
                                 method="ML")


#diagnostics
qqnorm(mod.totalA.GI.intron.apesOWMs$residuals, pch = 1, frame = FALSE)
qqline(mod.totalA.GI.intron.apesOWMs$residuals, col = "steelblue", lwd = 2)
plot(mod.totalA.GI.intron.apesOWMs$fitted, mod.totalA.GI.intron.apesOWMs$residuals)


#join with the table
activity_summary_hum1_withPheno_intron<-activity_summary_hum1_withPheno_intron %>%
  left_join(as.data.frame(mod.totalA.GI.intron.apesOWMs$fitted) %>%
              rownames_to_column("species") %>%
              dplyr::rename(fitted_log2_GI_apesOWMs="mod.totalA.GI.intron.apesOWMs$fitted"))


fig2C_3<-ggplot(activity_summary_hum1_withPheno_intron %>% filter(primate_clade %in% c("Great ape","Old World monkey")), 
                aes(x=log2_total_activity, y=log2(GI)))+
  geom_point(aes(color=clade), size=1.5, alpha=0.9)+theme_classic()+
  xlab("log2 (Intron activity)")+
  geom_line(aes(y=fitted_log2_GI_apesOWMs), color="grey50")+  
  scale_color_manual(values=clade_colors)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.text=element_text(size=7),
        axis.text.y=element_text(margin=margin(l=-5)),
        axis.text.x=element_text(margin=margin(b=-1.5)),
        legend.position="none",
        plot.margin = unit(c(0.15,-0.5,0.25,0), "cm"),
        panel.border = element_rect(fill = NA, colour = "grey60"),
        plot.title = element_text(family = "Helvetica", face = "plain", size = (6), hjust=0.5, margin=margin()))+
  ggtitle("Great apes &\n Old World monkeys")+
  scale_x_continuous(breaks=scales::breaks_extended(2))+
  scale_y_continuous(breaks=scales::breaks_extended(3))+
  annotate("text", label=paste0(round(rr2::R2.resid(mod.totalA.GI.intron.apesOWMs,mod.totalA.GI.null.apesOWMs), digits = 2)),  -Inf, Inf,  hjust = -0.3, vjust = 1.7, size=2.5)

fig2C_plots<-plot_grid(fig2C_1_final,fig2C_3, align="hv", ncol = 2)
fig2C_plots<-ggdraw(add_sub(fig2C_plots, "log2 (Intron activity)", vpadding=grid::unit(0,"lines"),y=4.6, x=0.52, vjust=4, size=8))
fig2C_final<-plot_grid(fig2C_plots, figC1_legend, ncol = 2, rel_widths = c(0.4,0.12))
fig2C_final





# Figure E: combined model ####
coevol_estimates_9sp_PGLS<-readRDS("regulation/data/MPRA/output/combined_model_PGLS.rds")

fig2E<-ggplot(coevol_estimates_9sp_PGLS, aes(x=fittedGI_dnds_intron, y=log2(GI)))+
  geom_abline(intercept = 0, slope=1, color="grey50")+
  geom_point(aes(color=clade), size=1.5, alpha=0.9)+
  theme_bw()+
  xlab("predicted log2(GI)")+
  scale_color_manual(values=clade_colors)+  
  scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x=element_text( size=8), 
        axis.title.y=element_text( size=8), 
        axis.text.y=element_text(size=7, margin=margin(l=-3)),
        axis.text.x=element_text(size=7),
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "grey60"))+
  annotate("text", label=paste0(round(coevol_estimates_9sp_PGLS$r2_dnds_intron[1], digits = 2)),  -Inf, Inf, hjust = -0.3, vjust = 1.7, size=2.5)


plot_grid(fig2C_final, fig2E, ncol=2, rel_widths = c(1,0.46), labels=c("d","e"), label_size = 11)
ggsave("figures/out/panel2/fig2DE.pdf", height=40*0.9, width=120*0.9, units="mm")









#..................................................................#####
# Supplementary panel 2 ####


# Extended Figs A-B: Detected tile distribution ####
chunk_distr<-readRDS("regulation/data/MPRA/output/tile_coverage_distr.rds")

chunk_distr$clade<-factor(chunk_distr$clade, levels=c("Primate","Rodent","Carnivore","Other"))
chunk_distr<-chunk_distr %>%
  mutate(order=case_when(clade=="Primate" ~ 1,
                         clade=="Rodent" ~ 2,
                         clade=="Carnivore" ~ 3,
                         T ~ 4))

ext2A<-ggplot(chunk_distr, aes(x=reorder(species,order), y=fract_present, fill=clade))+
  geom_boxplot(outlier.size = 0.4, lwd=0.3, alpha=0.9)+
  theme_bw()+
  scale_fill_manual(values=clade_colors)+
  ylab("Fraction of tiles\n present across regions")+
  theme(axis.title.x = element_text(size=9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=7),
        legend.position = "bottom",
        legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.spacing.y = unit(0, "cm"),
        legend.key.size = unit(0.35, "cm"),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-2,0,0,0),
        legend.box.margin=margin(-2,0,0,0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color="grey70"))+
  scale_y_continuous( minor_breaks =seq(0,1,0.249))+
  xlab("Species")




ext2B<-ggplot(chunk_distr, aes(x=region, y=fract_present, fill=region))+
  geom_boxplot(outlier.size = 0.4, lwd=0.3,alpha=0.9)+
  theme_bw()+
  scale_fill_manual(values=region_colors)+
  ylab("Fraction of tiles\n present across species")+
  theme(axis.title.x = element_blank(),
        axis.text.x=element_text(hjust=1, angle=45, size=8),
        axis.title.y = element_text(size=8),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_line(color="grey70"))+
  scale_y_continuous( minor_breaks =seq(0,1,0.249))

#1 outlier
d_ordii_cov<-chunk_distr %>% filter(species %in% "Dipodomys_ordii")

ext2A_B<-plot_grid(ext2A, ext2B, rel_widths = c(1,0.5), align = "hv", axis="b", scale = 0.95, labels = c("a","b"), label_size = 12)

ggsave("figures/out/extended2/ext2A_B.pdf", ext2A_B, width=183, height=75,units = "mm")




#Extended Fig C: tile activity correlation ####
median_activity_MPRA<-readRDS("regulation/data/MPRA/output/median_activity_MPRA.rds")

median_activity_MPRA_wide <- reshape2::dcast(median_activity_MPRA, SeqID ~ cell_line, value.var = "log2_activity", fill=0) %>%
  column_to_rownames("SeqID")

#human1 vs 2
cor.test(median_activity_MPRA_wide$human1, median_activity_MPRA_wide$human2)
#6D6875
activity_hum1_vs_hum2<-ggplot(median_activity_MPRA_wide, aes(x=human1, y=human2))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  xlab("Tile activity in human 1")+
  ylab("Tile activity in human 2")+
  ggtitle("Human1 vs Human2 (" ~rho~"0.5)")+
  ggpointdensity::geom_pointdensity()+
  scale_color_gradient2(low="#6D6875", mid="#7C7783", high="#FFCDB2")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-2,-5,0,0),
        legend.box.margin=margin(-2,-5,0,0),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(colour=guide_colourbar(barwidth=3, barheight=0.5))


#human1 vs macaque
cor.test(median_activity_MPRA_wide$human1, median_activity_MPRA_wide$macaque)
activity_hum1_vs_mac<-ggplot(median_activity_MPRA_wide, aes(x=human1, y=macaque))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  xlab("Tile activity in human 1")+
  ylab("Tile activity in macaque")+
  ggtitle("Human1 vs Macaque (" ~rho~"0.43)")+
  ggpointdensity::geom_pointdensity()+
  scale_color_gradient2(low="#6D6875", mid="#7C7783", high="#FFCDB2")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-2,-5,0,0),
        legend.box.margin=margin(-2,-5,0,0),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(colour=guide_colourbar(barwidth=3, barheight=0.5))


#human2 vs macaque
cor.test(median_activity_MPRA_wide$human2, median_activity_MPRA_wide$macaque)
activity_hum2_vs_mac<-ggplot(median_activity_MPRA_wide, aes(x=human1, y=macaque))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0)+
  theme_bw()+
  xlab("Tile activity in human 2")+
  ylab("Tile activity in macaque")+
  ggtitle("Human2 vs Macaque (" ~rho~"0.39)")+
  ggpointdensity::geom_pointdensity()+
  scale_color_gradient2(low="#6D6875", mid="#7C7783", high="#FFCDB2")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-2,-5,0,0),
        legend.box.margin=margin(-2,-5,0,0),
        legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(colour=guide_colourbar(barwidth=3, barheight=0.5))


#combine
ext2C<-plot_grid(activity_hum1_vs_hum2, activity_hum1_vs_mac, activity_hum2_vs_mac, ncol=3, scale=0.95)
ggsave("figures/out/extended2/ext2C.pdf", ext2C, width=183, height=70,units = "mm")






# Extended Fig D: Activity per region ####

activity_overlap_summary<-readRDS("regulation/data/MPRA/output/activity_overlap_summary.rds")

ggplot(activity_overlap_summary, aes(x=log2_total_activity, color=cell_line))+geom_density()

activity_overlap_summary_wide <-reshape2::dcast(activity_overlap_summary, cell_line ~ region +species , value.var = "log2_total_activity", fill=0) %>%
  column_to_rownames("cell_line")

activity_overlap_summary_wide<-as.data.frame(t(activity_overlap_summary_wide)) %>%
  rownames_to_column("region_species") %>%
  mutate(region=word(region_species,1,1, sep="_"))

activity_overlap_summary_wide$region<-factor(activity_overlap_summary_wide$region, levels=c("upstream1","upstream2","upstream3","exon1","intron","exon2","downstream"))


#pairwise
#human1 vs human2
cor.test(activity_overlap_summary_wide$human1, activity_overlap_summary_wide$human2)
hum1_vs_hum2_region<-ggplot(activity_overlap_summary_wide, aes(x=human1, y=human2, color=region))+
  geom_abline(slope = 1, intercept = 0, color="grey60")+
  geom_point(alpha=0.7, size=1)+
  scale_color_manual(values=region_colors)+
  theme_bw()+
  xlab("Region activity in human 1")+
  ylab("Region activity in human 2")+
  ggtitle("Human1 vs Human2 (" ~rho~"0.86)")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

#human1 vs macaque
cor.test(activity_overlap_summary_wide$human1, activity_overlap_summary_wide$macaque)
hum1_vs_mac_region<-ggplot(activity_overlap_summary_wide, aes(x=human1, y=macaque, color=region))+
  geom_abline(slope = 1, intercept = 0, color="grey60")+
  geom_point(alpha=0.7, size=1)+
  scale_color_manual(values=region_colors)+
  theme_bw()+
  xlab("Region activity in human 1")+
  ylab("Region activity in macaque")+
  ggtitle("Human1 vs Macaque (" ~rho~"0.85)")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

#human2 vs macaque
cor.test(activity_overlap_summary_wide$human2, activity_overlap_summary_wide$macaque)
hum2_vs_mac_region<-ggplot(activity_overlap_summary_wide, aes(x=human2, y=macaque, color=region))+
  geom_abline(slope = 1, intercept = 0, color="grey60")+
  geom_point(alpha=0.7, size=1)+
  scale_color_manual(values=region_colors)+
  theme_bw()+
  xlab("Region activity in human 2")+
  ylab("Region activity in macaque")+
  ggtitle("Human2 vs Macaque (" ~rho~"0.88)")+
  theme(axis.text=element_text(size=7),
        axis.title=element_text(size=8),
        plot.title = element_text(size=8),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,-2,0,0),
        legend.box.margin=margin(0,-2,0,0),
        #legend.position = "bottom",
        legend.title = element_text(size=7),
        legend.text = element_text(size=5),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.size = unit(0.35, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


ext2D<-plot_grid(hum1_vs_hum2_region, hum1_vs_mac_region, hum2_vs_mac_region, ncol=3, rel_widths = c(0.3,0.3,0.385), scale=0.95, labels = c("d",NULL,NULL), label_size = 12)
ggsave("figures/out/extended2/ext2D.pdf", ext2D, width=183, height=55,units = "mm")

ext2C_2<-plot_grid(ext2C,NULL, rel_widths = c(1,0.1), labels = c("c",NULL), label_size = 12)
ext2<-plot_grid(ext2A_B, ext2C_2, ext2D, ncol=1, rel_heights = c(1,0.95,0.83))

ggsave("figures/out/extended2/ext2.pdf", ext2_v2, width=183, height=185,units = "mm")







#..................................................................#####
# Supplementary panel 3 ####

# Region length across the tree ####
dnase_species<-activity_overlap_summary  %>%
  dplyr::select(species, clade) %>%
  dplyr::distinct(species, .keep_all=T) %>%
  mutate(species=gsub("_"," ",species),
        DNase_species=case_when(species=="Mus musculus" ~ "DHS species",
                                 species=="Homo sapiens" ~ "DHS species",
                                 T ~ clade),
         face=case_when(species %in% c("Mus musculus","Homo sapiens") ~ "bold.italic",  T ~ "italic"))

dnase_species$DNase_species<-factor(dnase_species$DNase_species, levels=c("DHS species","Primate","Rodent","Carnivore", "Other"))

mpra_tree_full<-ape::drop.tip(mammaltree, mammaltree$tip.label[!mammaltree$tip.label %in% activity_overlap_summary$species])
mpra_tree_full$tip.label<-gsub("_"," ", mpra_tree_full$tip.label)


p_full_tree0<-ggtree(mpra_tree_full, size=0.2) 
p_full_tree<-p_full_tree0 %<+% dnase_species +
  geom_tiplab(aes(color=factor(DNase_species), fontface=face), size=2,
              geom = "text")+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-10,0,0,0),
        legend.box.margin=margin(-10,0,0,0),
        legend.text = element_text(size=8),
        #legend.spacing.x = unit(0.4, "cm"),
        legend.key.size = unit(0.8, "cm"),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  scale_color_manual(values=c(clade_colors_dnase))




cov_region_length_wide <- reshape2::dcast(activity_overlap_summary %>% 
                                  filter(cell_line=="human1"), 
                                species ~ region , value.var = "length", fill=NA) %>%
  mutate(species=gsub("_"," ", species)) %>%
  left_join(dnase_species) %>%
  mutate(dnase=case_when(species=="Mus musculus" ~ "DHS species",
                                 species=="Homo sapiens" ~ "DHS species",
                                 T ~ clade)) 


#how does it know the right order actually?
p1<-facet_plot(p_full_tree, panel="Upstream 1", data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=upstream1, y=y, yend=y, colour=as.factor(dnase))) +
  #theme_tree2()+
  xlim_expand(c(0,300),"Tree")+
  theme(strip.text = element_text(size=7.5))
p2<-facet_plot(p1, panel="Upstream 2",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=upstream2, y=y, yend=y, colour=as.factor(dnase)))
p3<-facet_plot(p2, panel="Upstream 3",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=upstream3, y=y, yend=y, colour=as.factor(dnase)))
p4<-facet_plot(p3, panel="Exon 1",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=exon1, y=y, yend=y, colour=as.factor(dnase)))
p5<-facet_plot(p4, panel="Intron",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=intron, y=y, yend=y, colour=as.factor(dnase)))
p6<-facet_plot(p5, panel="Exon 2",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=exon2, y=y, yend=y, colour=as.factor(dnase)))
p7<-facet_plot(p6, panel="Downstream",  data=cov_region_length_wide, geom=geom_segment, aes(x=0, xend=downstream, y=y, yend=y, colour=as.factor(dnase)))


# increase the relative width of the tree compared to others
gt = ggplot_gtable(ggplot_build(p7))
gtable_show_layout(gt) # will show you the layout
gt # see plot layout in table format
gt$layout$l[grep('panel-1', gt$layout$name)] # find the column specific to panel-1
gt$widths[5] = 4*gt$widths[5] # it was column 5 - increase the width 4x

pdf("figures/out/extended3/ext3.pdf", paper="a4", width=0, height=6.7)
par(mar=c(0,0,0,0))
grid.draw(gt) # plot with grid draw
dev.off()





#..................................................................#####
# Supplementary panel 4 ####

#Extended Fig A: a boxplot showing activity distribution ####
median_length_df<-activity_summary_hum1_withPheno %>%
  group_by(region) %>%
  summarize(median_length=median(length)/1000)

ext4a<-activity_summary_hum1_withPheno %>%
  left_join(median_length_df) %>%
  ggplot(aes(x=region, y=log2_total_activity)) + 
  geom_boxplot(aes(fill=median_length), lwd=0.01, outlier.size=0.1, color="grey50")+
  theme_bw()+
  scale_fill_gradient(low="grey80",high="#FAA275")+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=7, colour="grey15"),
        axis.text.y = element_text(hjust = 0.95, size=7),
        axis.ticks.x = element_blank(),
        legend.spacing.x = unit(1.5,"mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=5,margin = margin(r = -3, l=-3)),
        legend.title=element_text(size=6, margin=margin(r=-3, l=-3), color="grey25"),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(0,-3,0,-1),
        legend.box.margin=margin(0,-3,0,-1),
        panel.border = element_rect(fill = NA, 
                                    colour = "grey60"),
        plot.margin = unit(c(-0.1, 0, 0.3, 0), "cm"))+
  guides(fill = guide_colorbar( barwidth = 0.4,barheight = 1.5))+
  scale_y_continuous(breaks=scales::breaks_extended(3), limits=c(-1,4.6))+
  scale_x_discrete(labels=c("U1","U2","U3","E1","I","E2","D"))+
  labs(fill="Median\nlength")+
  ylab("log2 (Regulatory\n activity)")





# Extended Fig B: TF expression ####
intron_tf_expr_NPCs_long<-readRDS("regulation/data/TFs/expression/vsdMat_22TFs.rds")
ext4b<-ggplot(intron_tf_expr_NPCs_long, aes(x=reorder(SYMBOL, norm_expr), y=norm_expr))+
  geom_boxplot(color="grey40", size=0.3)+geom_point(aes(color=clone))+
  theme_bw()+
  coord_flip()+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=8),
        axis.text.y=element_text(size=7.2,face=c("plain","plain","plain","plain","bold","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain","plain")),
        legend.position = "bottom",
        legend.spacing.x = unit(0.5,"mm"),
        legend.text = element_text(size=6,margin = margin()),
        legend.title=element_text(size=6),
        legend.background = element_rect(fill="transparent",colour=NA),
        legend.margin=margin(-5,0,0,-6),
        legend.box.margin=margin(-5,0,0,-6))+
  ylab("Normalized expression")+
  labs(color="Cell line")+
  scale_color_manual(labels=c("human1","human2","macaque"),values=c("#1B998B","#04516B","#8B786D"))




# Extended Fig C: TF motif distances ####
motif_distances_22_symb<-readRDS("regulation/data/TFs/ClusterBuster/results/motif_distances_22_symb.rds")

average_dist<-motif_distances_22_symb %>% 
  rownames_to_column("SYMBOL1") %>%
  gather(SYMBOL2,distance,-SYMBOL1) %>%
  filter(SYMBOL1!=SYMBOL2) %>%
  summarise(mean_dist=mean(distance))

hc <- flashClust::hclust(as.dist(motif_distances_22_symb))
ext4c<-ggdendro::ggdendrogram(hc, rotate = FALSE)+ 
  geom_hline(yintercept = average_dist$mean_dist, linetype="dashed", color="grey50")+  
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
  theme(axis.line.y = element_line(size=0.2),
        axis.text.x = element_text(size=7.2),
        axis.text.y=element_text(size=6.5))
ext4c



# Extended Fig D&E: Predicted activity ####
combined_coeffs<-readRDS("regulation/data/TFs/ClusterBuster/results/PGLS_3TF_coeffs.rds")

top<-ggplot(combined_coeffs %>% filter(motif!="(Intercept)"), 
            aes(x=motif, y=Value, color=y, group=y))+
  geom_point(size=1.8,position=position_dodge(width=0.9))+
  geom_errorbar(aes(ymin=Value-Std.Error, ymax=Value+Std.Error), 
                width=0.3, position=position_dodge(width=0.9))+
  coord_flip()+geom_abline(slope=0,intercept=0, linetype="dashed")+
  ylab("Coefficient")+
  theme_bw()+
  scale_color_manual(values=c("#40798E","#FCC25F"), labels=c("GI","intron\nactivity"))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(), legend.position = "right",
    axis.title = element_text(size=7.5),
    axis.text=element_text(size=7, color="black"),
    axis.text.y=element_text(margin=margin(l=-3)),
    axis.text.x=element_text(margin=margin(b=-1.5)),
    legend.spacing.y = unit(-3, "mm"),
    legend.text = element_text(size=6,margin = margin(l=-3, b=-3,t=-3)),
    legend.background = element_rect(fill="transparent",colour=NA),
    legend.margin=margin(0,-3,-3,-3),
    legend.box.margin=margin(0,-3,-3,-3),
    legend.title=element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"))
top

scores_top22TFs_mat<-readRDS("regulation/data/TFs/ClusterBuster/results/scores_top22TFs_mat.rds")
act_pp2<-ggplot(scores_top22TFs_mat, aes(x=fitted_intron_3TFs, y=log2_total_activity,color=clade))+
  geom_abline(intercept = 0,slope=1, color="grey50")+
  geom_point(size=1.8)+
  scale_color_manual(values=clade_colors)+ 
  theme_bw()+
  xlab("Predicted log2(intron activity)")+
  ylab("log2(intron activity)")+
  scale_x_continuous(breaks=scales::pretty_breaks(n=3))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  annotate("text", label=0.78,  -Inf, Inf,  hjust = -0.3, vjust = 1.7, size=2.5)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title = element_text(size=8),
    axis.text=element_text(size=7),
    axis.text.y=element_text(margin=margin(l=-3)),
    axis.text.x=element_text(margin=margin(b=-1.5)),
    plot.margin = unit(c(0,0,0,0), "cm"))

ext4e<-plot_grid(top, act_pp2, ncol = 1, align = "hv", axis="lbr", rel_heights = c(0.65,1))


#Extended Fig F: CTCF Chip-Seq in mouse ####
ext4f<- ggdraw() + 
  cowplot::draw_image(magick::image_read("figures/out/extended4/ctcf_chipseq2.png"))




#combine all####
ggdraw() +
  draw_plot(ext4a, x = 0.01, y = .83, width = .5, height = .15) +
  draw_plot(ext4b, x = 0.03, y = .25, width = .45, height = .58) +
  draw_plot(ext4e, x = 0.59, y = 0.29, width = .33, height = .35) +
  draw_plot(ext4c, x = 0.55, y = 0.64, width = 0.45, height = 0.35) +
  draw_plot(ext4f, x = 0.06, y = 0.01, width = 0.6+0.3, height = .25)+
  draw_plot_label(label = c("a", "c", "b", "d","e","f"), size = 11,
                  x = c(0.05, 0.54, 0.05, 0.55, 0.55, 0.05 ), 
                  y = c(1, 1, 0.86,0.69,0.53,0.28 ))

ggsave("figures/out/extended4/ext4.pdf",  width=140, height=150,units = "mm")

