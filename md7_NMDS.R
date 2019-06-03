color.DaFILF<-c("#e50000","#ac00e6","#3333ff","#0ba29a","#00cc00","#ff9933","#ff00aa","#F00000","#ffff00","#800055")
color.DaFInt<-c("#e50000","#ac00e6","#3333ff","#0ba29a","#ff9933","#ff00aa","#F00000","#800055")

###################################################################################
########## FIGURE 1 ##########
#### Host Species and Sample Type NMDS ###
dsNMDS<-phyT.leech # define data for analysis
sample_data(dsNMDS)$Sample_Type = factor(sample_data(dsNMDS)$Sample_Type, levels = c("ILF","Intestinum","Bladder"))
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination

nmds.SpeciesType<-plot_ordination(dsNMDS, ord, shape="Taxonomic_ID",color="Sample_Type") + 
  geom_point(aes(fill=Sample_Type),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position="none")
nmds.SpeciesType # uncomment to print plot

### Host Species and Sample Type by Order NMDS###
dsNMDSorder<-tax_glom(dsNMDS,taxrank="Order")
distOrdO = phyloseq::distance(dsNMDSorder, method = c(mNMDS)) # calculate distances
ordO = ordinate(dsNMDSorder, method = "NMDS", distance = distOrdO) # calculate ordination
nmds.SpeciesTypeOrder<-plot_ordination(dsNMDSorder, ordO, shape="Taxonomic_ID",color="Sample_Type") + 
  geom_point(aes(fill=Sample_Type),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position="none")
#nmds.SpeciesTypeOrder # uncomment to print plot
ggsave(plot_grid(nmds.SpeciesType, nmds.SpeciesTypeOrder, labels = "AUTO"), filename="NMDS/plotNMDS_fig1.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 1 to .eps file

###################################################################################

########## FIGURE 4 ##########
### ILF/WildTime NMDS ###
dsNMDS<-subset_samples(phyT.md,Sample_Type=="ILF") # define data for analysis
sample_data(dsNMDS)$WildMonth = factor(sample_data(dsNMDS)$WildMonth, levels = c("April","June","July","August","September","October")) # force Months to be in chronological order (instead of alphabetical default)
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination
nmds.collILF<-plot_ordination(dsNMDS, ord, color = "WildMonth") + 
  geom_point(aes(fill=WildMonth),size=3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#9bcdff","#ffb4b4","#ff0101","#810000","#810000","#005ab4")) +
  theme_bw() +
  theme(text=element_text(size=10), legend.position="none")
nmds.collILF # uncomment to print plot

### Intestinum/WildTime NMDS ###
dsNMDSint<-subset_samples(phyT.md,Sample_Type=="Intestinum") # define data for analysis
sample_data(dsNMDSint)$WildMonth = factor(sample_data(dsNMDSint)$WildMonth, levels = c("April","June","July","August","September","October")) # force Months to be in chronological order (instead of alphabetical default)
distOrdInt = phyloseq::distance(dsNMDSint, method = c(mNMDS)) # calculate distances
ordInt = ordinate(dsNMDSint, method = "NMDS", distance = distOrdInt) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.collInt<-plot_ordination(dsNMDSint, ordInt, color = "WildMonth") + 
  geom_point(aes(fill=WildMonth),size=3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#9bcdff","#ffb4b4","#ff0101","#810000","#810000","#005ab4")) +
  theme_bw() +
  theme(text=element_text(size=10), legend.position="none")
#nmds.collInt # uncomment to print plot
ggsave(plot_grid(nmds.collILF, nmds.collInt, labels = "AUTO"), filename="NMDS/plotNMDS_fig4.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 4 to .eps file
###################################################################################




########## Sample Type NMDS ##########
dsNMDS<-subset_samples(phyT.md) # define data for analysis
sample_data(dsNMDS)$Sample_Type = factor(sample_data(dsNMDS)$Sample_Type, levels = c("ILF","Intestinum","Bladder"))
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.Type<-plot_ordination(dsNMDS, ord, shape="Sample_Type") + 
  geom_point(size=6,alpha=0.75) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.99, linetype = 2) +
  scale_shape_manual(values=c(15,16,4),name="Organ Sampled") +
  scale_color_grey() +
  theme_bw()
nmds.Type
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.Type), size = "last")), filename="NMDS/plotNMDS_sampleType.eps", device="eps", width=10,height=8)




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

########## Macrobdella DaF ILF NMDS ##########
phyT.mdDaFILF<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF")) # Keep only Intestinum samples
sample_data(phyT.mdDaFILF)$Da1Fb = factor(sample_data(phyT.mdDaFILF)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1F
#dsNMDSdaf<-subset_samples(phyT.mdDaFILF,AnimalSource=="Wlot")
dsNMDSdaf<-phyT.mdDaFILF
mNMDS.daf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDS.daf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFILF<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1Fb") + 
  geom_point(size=6,mapping = aes(color=Da1Fb)) +
  ggtitle(bquote("ILF" + .(c(mNMDS.daf)))) +
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFILF) +
  theme_bw()
nmds.DaFILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.DaFILF), size = "last")), filename="NMDS/plotNMDS_Da1FILF.eps", device="eps", width=8,height=8)

########## Hirudo DaF ILF NMDS ##########
phyT.hvILF<-subset_samples(phyT.Hv,Sample_Type%in%c("ILF")) # Keep only Intestinum samples
sample_data(phyT.hvILF)$Da1Fb = factor(sample_data(phyT.hvILF)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1F
#dsNMDSdaf<-subset_samples(phyT.mdDaFILF,AnimalSource=="Wlot")
dsNMDSdaf<-phyT.hvILF
mNMDS.daf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDS.daf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.hvDaFILF<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1Fb") + 
  geom_point(size=6,mapping = aes(color=Da1Fb)) +
  ggtitle(bquote("Hirudo verbana ILF" + .(c(mNMDS.daf)))) +
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFILF) +
  theme_bw()
nmds.hvDaFILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmds.hvDaFILF), size = "last")), filename="NMDS/plotNMDS_hvDa1FILF.eps", device="eps", width=8,height=8)

########## DaF Intestinum NMDS ##########
phyT.mdDaFInt<-subset_samples(phyT.mdDaF,Sample_Type%in%c("Intestinum")) # Keep only Intestinum samples
sample_data(phyT.mdDaFInt)$Da1Fb = factor(sample_data(phyT.mdDaFInt)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1F
dsNMDSdaf<-phyT.mdDaFInt
mNMDS.daf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDS.daf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmds.DaFInt<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1Fb") + 
  geom_point(size=6,mapping = aes(color=Da1Fb)) +
  ggtitle(bquote("Intestinum" + .(c(mNMDS.daf)))) +
  stat_ellipse(type = "t", level = 0.75, linetype = 2) +
  scale_color_manual(values=color.DaFInt) +
  theme_bw()
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