color.DaFILF<-c("#e50000","#ac00e6","#3333ff","#0ba29a","#00cc00","#ff9933","#ff00aa","#F00000","#ffff00","#800055")
color.DaFInt<-c("#e50000","#ac00e6","#3333ff","#0ba29a","#ff9933","#ff00aa","#F00000","#800055")

########## ILF/AnimalSource NMDS ##########
dsNMDS<-subset_samples(phyT.mdILF,!Da1F%in%c("4")) # define data for analysis
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.Geo<-plot_ordination(dsNMDS, ord, color = "AnimalSource") + 
  geom_point(size=6,mapping = aes(color=AnimalSource)) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1"))
nmds.Geo
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.Geo), size = "last")), filename="NMDS/plotNMDS_ILFgeo.eps", device="eps", width=8,height=8)

########## Sample_Type NMDS ##########
### Find outlier samples ###
phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora")
phyT.mdMain<-subset_samples(phyT.md,!AnimalSource%in%c("MtSnowVT","CarogaNY"))
phyT.mdMainBlad<-subset_samples(phyT.mdMain,Sample_Type=="Bladder")
phyT.mdMainGI<-subset_samples(phyT.mdMain,Sample_Type!="Bladder")
phyT.mdUnfed<-subset_samples(phyT.mdMainGI,!Da1F%in%c("2","4","7"))
phyT.mdMerge<-merge_phyloseq(phyT.mdMainBlad,phyT.mdUnfed)

phyT.NMDStissue<-phyT.mdMerge #define data for analysis
m.NMDStissue<-"unifrac" # define metric for analysis
# calculate distances
dist.NMDStissue = phyloseq::distance(phyT.NMDStissue, method = c(m.NMDStissue)) # calculate distances
ord.NMDStissue = ordinate(phyT.NMDStissue, method = "NMDS", distance = dist.NMDStissue) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsType<-plot_ordination(phyT.NMDStissue, ord.NMDStissue, color = "Sample_Type") + 
  geom_point(size=5,mapping = aes(color=Sample_Type)) +
  ggtitle(c(m.NMDStissue)) + 
  stat_ellipse(type = "t", level = 0.7, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsType
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmdsType), size = "last")), filename="NMDS/plotNMDS_type.eps", device="eps", width=10,height=8)
##### Figure #####

########## DaF ILF NMDS ##########
phyT.mdDaFILF<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF")) # Keep only Intestinum samples
sample_data(phyT.mdDaFILF)$Da1Fb = factor(sample_data(phyT.mdDaFILF)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1F
#dsNMDSdaf<-subset_samples(phyT.mdDaFILF,AnimalSource=="Wlot")
dsNMDSdaf<-phyT.mdDaFILF
mNMDS.daf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDS.daf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "NMDS", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFILF<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1Fb") + 
  geom_point(size=6,mapping = aes(color=Da1Fb)) +
  ggtitle(bquote("ILF" + .(c(mNMDS.daf)))) +
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFILF) 
nmds.DaFILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.DaFILF), size = "last")), filename="NMDS/plotNMDS_Da1FILF.eps", device="eps", width=8,height=8)

########## DaF Intestinum NMDS ##########
phyT.mdDaFInt<-subset_samples(phyT.mdDaF,Sample_Type%in%c("Intestinum")) # Keep only Intestinum samples
sample_data(phyT.mdDaFInt)$Da1Fb = factor(sample_data(phyT.mdDaFInt)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1F
dsNMDSdaf<-phyT.mdDaFInt
mNMDS.daf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDS.daf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "NMDS", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFInt<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1Fb") + 
  geom_point(size=6,mapping = aes(color=Da1Fb)) +
  ggtitle(bquote("Intestinum" + .(c(mNMDS.daf)))) +
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFInt) 
nmds.DaFInt
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.DaFInt), size = "last")), filename="NMDS/plotNMDS_Da1FInt.eps", device="eps", width=8,height=8)

########## Macrobdella Sample_Type/AnimalSource NMDS ##########
phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
#phyT.mdCM<-subset_samples(phyT.md,AnimalSource%in%c("Wlot","GrotonMA")) # subset to include only MA and CT samples (phyloseq transformed . macrobdella decora CT MA)
phyT.md0<-subset_samples(phyT.md,!Da1Fb%in%c("2","4","7"))
phyT.mdGI<-subset_samples(phyT.md0,Sample_Type!="Bladder") # subset to exclude samples (phyloseq transformed . macrobdella decora gastrointestinal)

phyT.NMDSmdGI<-phyT.mdGI
m.NMDSmdGI<-"unifrac" # define metric for analysis
dist.NMDSmdGI = phyloseq::distance(phyT.NMDSmdGI, method = c(m.NMDSmdGI)) # calculate distances
ord.NMDSmdGI = ordinate(phyT.NMDSmdGI, method = "NMDS", distance = dist.NMDSmdGI) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsMacro<-plot_ordination(phyT.NMDSmdGI, ord.NMDSmdGI, color = "Header",shape="Sample_Type") + 
  geom_point(size=5,aes(color=Header,shape=Sample_Type)) + 
  scale_shape_manual(values = c(16,1)) +  
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  ggtitle(c(mNMDS)) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsMacro
ggsave(grid.draw(rbind(ggplotGrob(nmdsMacro), size = "last")), filename="NMDS/plotNMDS_ILFInt.eps", device="eps", width=10,height=8)


########## Species NMDS ##########
### ILF ###
phyT.total<-merge_phyloseq(phyT.Md,phyT.Hv) # merge Macrobdella and Hirudo phyloseq objects (physeq transformed . total)
phyT.total0<-subset_samples(phyT.total,!Da1Fb%in%c("2","4","7"))
phyT.totalILF<-subset_samples(phyT.total0,Sample_Type=="ILF")
dsNMDS.ILF<-subset_samples(phyT.totalILF,!Da1F%in%c("4")) #define data for analysis
mNMDS.ILF<-"wunifrac" # define metric for analysis

distOrd.ILF = phyloseq::distance(dsNMDS.ILF, method = c(mNMDS.ILF)) # calculate distances
ord.ILF = ordinate(dsNMDS.ILF, method = "NMDS", distance = distOrd.ILF) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
pNMDS.ILF<-plot_ordination(dsNMDS.ILF, ord.ILF, color = "Taxonomic_ID") + 
  geom_point(size=6,mapping = aes(color=Taxonomic_ID)) +
  ggtitle(bquote("ILF" + .(c(mNMDS.ILF)))) +
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
pNMDS.ILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pNMDS.ILF), size = "last")), filename="NMDS/plotNMDS_speciesILF.eps", device="eps", width=8,height=8)

### Intestinum ###
phyT.total<-merge_phyloseq(phyT.Md,phyT.Hv) # merge Macrobdella and Hirudo phyloseq objects (physeq transformed . total)
phyT.total0<-subset_samples(phyT.total,!Da1Fb%in%c("2","4","7"))
phyT.totalInt<-subset_samples(phyT.total0,Sample_Type=="Intestinum")
dsNMDS.Int<-subset_samples(phyT.totalInt,!Da1F%in%c("4")) #define data for analysis

mNMDS.Int<-"wunifrac" # define metric for analysis
distOrd.Int = phyloseq::distance(dsNMDS.Int, method = c(mNMDS.Int)) # calculate distances
ord.Int = ordinate(dsNMDS.Int, method = "NMDS", distance = distOrd.Int) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
pNMDS.Int<-plot_ordination(dsNMDS.Int, ord.Int, color = "Taxonomic_ID") + 
  geom_point(size=6,mapping = aes(color=Taxonomic_ID)) +
  ggtitle(bquote("Intestinum" + .(c(mNMDS.Int)))) +
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
pNMDS.Int
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pNMDS.Int), size = "last")), filename="NMDS/plotNMDS_speciesInt.eps", device="eps", width=8,height=8)
##### Figure #####






###############################################
################ Practice area ################ 
############################################### 

### Md Blader ###
dsNMDS.Int<-phyT.baseBmd #define data for analysis

mNMDS.Int<-"wunifrac" # define metric for analysis
distOrd.Int = phyloseq::distance(dsNMDS.Int, method = c(mNMDS.Int)) # calculate distances
ord.Int = ordinate(dsNMDS.Int, method = "NMDS", distance = distOrd.Int) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
pNMDS.Int<-plot_ordination(dsNMDS.Int, ord.Int, label="EmilyID") + 
  geom_point(size=1,mapping = aes()) +
  ggtitle(bquote("Bladder" + .(c(mNMDS.Int)))) +
  stat_ellipse(type = "t", level = 0.9, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
pNMDS.Int
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pNMDS.Int), size = "last")), filename="NMDS/plotNMDS_speciesInt.eps", device="eps", width=8,height=8)
##### Figure #####