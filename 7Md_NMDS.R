color.DaFILF<-c("#ff00aa","#ac00e6","#3333ff","#0ba29a","#ffff00","#ff9933","#F00000","#800055")
color.DaFInt<-c("#ff00aa","#ac00e6","#3333ff","#0ba29a","#ff9933","#F00000","#800055")


########## ILF/AnimalSource NMDS ##########
dsNMDS<-subset_samples(phy.mdILF,Da1F!="2") # define data for analysis
dsNMDS<-phy.mdILF # define data for analysis

mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "PCoA", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.Geo<-plot_ordination(dsNMDS, ord, color = "AnimalSource") + 
  geom_point(size=6,mapping = aes(color=AnimalSource)) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1"))
nmds.Geo
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.Geo), size = "last")), filename="NMDS/plotNMDS_ILFgeo.png", width=8,height=8)

########## Sample_Type NMDS ##########
### Find outlier samples ###
phy.mdIUB<-subset_samples(phy.Md,Sample_Type!="Ovary")
dsNMDS2<-subset_samples(phy.mdIUB,!Da1F%in%c("2","4")) #define data for analysis
dsNMDS2<-prune_taxa(ls.coreTot,dsNMDS2)
mNMDS2<-"wunifrac" # define metric for analysis
# calculate distances
distOrd2 = phyloseq::distance(dsNMDS2, method = c(mNMDS2)) # calculate distances
ord2 = ordinate(dsNMDS2, method = "PCoA", distance = distOrd2) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsType<-plot_ordination(dsNMDS2, ord2, color = "Sample_Type") + 
  geom_point(size=6,mapping = aes(color=Sample_Type)) +
  ggtitle(c(mNMDS2)) + 
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsType
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmdsType), size = "last")), filename="NMDS/plotNMDS_type.png", width=8,height=8)
##### Figure #####

########## DaF ILF NMDS ##########
sample_data(phy.mdDaFILF)$Da1F = factor(sample_data(phy.mdDaFILF)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F
dsNMDSdaf<-subset_samples(phy.mdDaFILF,Da1F!="2")
mNMDSdaf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDSdaf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFILF<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1F") + 
  geom_point(size=6,mapping = aes(color=Da1F,shape=AnimalSource)) +
  ggtitle(c(mNMDSdaf)) + 
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFILF) 
nmds.DaFILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.DaFILF), size = "last")), filename="NMDS/plotNMDS_Da1FILF.png", width=8,height=8)

########## DaF Intestinum NMDS ##########
sample_data(phy.mdDaFInt)$Da1F = factor(sample_data(phy.mdDaFInt)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F
dsNMDSdaf<-subset_samples(phy.mdDaFInt,Da1F!="2")
mNMDSdaf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDSdaf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFInt<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1F") + 
  geom_point(size=6,mapping = aes(color=Da1F,shape=AnimalSource)) +
  ggtitle(c(mNMDSdaf)) + 
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFInt) 
nmds.DaFInt
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.DaFInt), size = "last")), filename="NMDS/plotNMDS_Da1FInt.png", width=8,height=8)

########## Macrobdella Sample_Type/AnimalSource NMDS ##########
macroBaseI<-phy.mdDaFGI
DistBC = phyloseq::distance(macroBaseI, method = "unifrac") # calculate distances
ordBC = ordinate(macroBaseI, method = "PCoA", distance = DistBC) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsMacro<-plot_ordination(macroBaseI, ordBC, color = "Sample_Type",shape="AnimalSource") + 
  geom_point(size=8,mapping = aes(color=Sample_Type,shape=AnimalSource)) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  ggtitle("PCoA: Bray-Curtis") +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsMacro
ggsave(grid.draw(rbind(ggplotGrob(nmdsMacro), size = "last")), filename="NMDS/plotNMDS_ILFInt.png", width=8,height=8)




#### May need in future ####
sdat.mdILF<-sample_data(phy.mdILF)
sdat.mdILF$fedState<-with(sdat.mdILF,
                          ifelse(Da1F=="0","starve",
                                 ifelse(Da1F=="1","fed",
                                        ifelse(Da1F=="2","fed",
                                               ifelse(Da1F=="4","fed",
                                                      ifelse(Da1F=="7","fed",
                                                             ifelse(Da1F=="31","fed",
                                                                    ifelse(Da1F=="35","old",
                                                                           ifelse(Da1F=="100","old",
                                                                                  ifelse(Da1F=="113","old",
                                                                                         ifelse(Da1F=="215","starve",
                                                                                                as.character(Da1F))))))))))))
phy.mdILF.mod<-phyloseq(otu_table(phy.mdILF),tax_table(phy.mdILF),sdat.mdILF,phy_tree(phy.mdILF))
