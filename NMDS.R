########## ILF/AnimalSource NMDS ##########
dir.create("NMDS") # create folder for data output(s)

macroA<-subset_samples(prAdult,Taxonomic_ID=="Mdecora") # keep only Macrobdella samples
#macroA<-subset_samples(macroAa,sample_names(macroAa)!="MAa0dF110814EMb.iD449") # remove obvious outlier that seems to be pure bladder
dsNMDS<-subset_samples(macroA,Sample_Type=="ILF") # keep only ILF samples. define data for analysis
mNMDS<-"unifrac" # define metric for analysis
# calculate distances
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "PCoA", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsGeo<-plot_ordination(dsNMDS, ord, color = "AnimalSource") + 
  geom_point(size=6,mapping = aes(color=AnimalSource)) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=brewer.pal(7,"Set1")) 
nmdsGeo
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmdsGeo), size = "last")), filename="NMDS/plotNMDSgeo.png", width=8,height=8)
##### Figure #####

########## Sample_Type NMDS ##########
### Find outlier samples ###
macroA<-subset_samples(corePhy,Taxonomic_ID=="Mdecora") # keep only Macrobdella samples
dsNMDS2<-macroA #define data for analysis
mNMDS2<-"wunifrac" # define metric for analysis
# calculate distances
distOrd2 = phyloseq::distance(dsNMDS2, method = c(mNMDS2)) # calculate distances
ord2 = ordinate(dsNMDS2, method = "PCoA", distance = distOrd2) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsType<-plot_ordination(dsNMDS2, ord2, color = "Sample_Type") + 
  geom_point(size=6,mapping = aes(color=Sample_Type)) +
  ggtitle(c(mNMDS2)) + 
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsType


########## DaF NMDS ##########
#dsNMDSdaf<-subset_samples(macDaFpr,Sample_Type=="Intestinum") #define data for analysis
sample_data(macDaF)$Da1F = factor(sample_data(macDaF)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F
dsNMDSdaf<-macDaF
mNMDSdaf<-"unifrac" # define metric for analysis
# calculate distances
distOrdDaf = phyloseq::distance(dsNMDSdaf, method = c(mNMDSdaf)) # calculate distances
OrdDaf = ordinate(dsNMDSdaf, method = "PCoA", distance = distOrdDaf) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsDaF<-plot_ordination(dsNMDSdaf, OrdDaf, color = "Da1F") + 
  geom_point(size=6,mapping = aes(color=Da1F)) +
  ggtitle(c(mNMDSdaf)) + 
  stat_ellipse(type = "t", level = 0.5, linetype = 2) +
  scale_color_manual(values=rainbow) 
nmdsDaF

##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(nmdsDaF), size = "last")), filename="plotNMDSdaf.png", width=8,height=8)
##### Figure #####


########## Macrobdella Sample_Type/AnimalSource NMDS ##########
macroBase<-subset_samples(macDaF,Taxonomic_ID=="Mdecora")
#macroBaseD<-subset_samples(macroBase,!Da1F=="0")
macroBaseI<-subset_samples(macroBase,Sample_Type%in%c("ILF","Intestinum"))
DistBC = phyloseq::distance(macroBaseI, method = "unifrac") # calculate distances
ordBC = ordinate(macroBaseI, method = "PCoA", distance = DistBC) # calculate ordination

# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsMacro<-plot_ordination(macroBaseI, ordBC, color = "Sample_Type") + 
  geom_point(size=8,mapping = aes(color=Sample_Type)) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  ggtitle("PCoA: Bray-Curtis") +
  scale_color_manual(values=brewer.pal(6,"Set1")) 
nmdsMacro
ggsave(grid.draw(rbind(ggplotGrob(nmdsMacro), size = "last")), filename="plotNMDSmacro.png", width=8,height=8)




