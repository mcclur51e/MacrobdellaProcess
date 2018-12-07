############## Ordination/Beta Diversity #################
#library("phyloseq") # should already be loaded
#library("ggplot2") # should already be loaded
library("vegan") # should already be loaded
library("plyr") # should already be loaded

phyR.adult<-subset_samples(phy.out,Age=="A") # keep adult samples (prune Adult)
phyR.md<-subset_samples(phyR.adult,Taxonomic_ID%in%c("Mdecora","Unk")) # subset containing only Macrobdella samples

sdat.md<-sample_data(phyR.md)
sdat.md$Da1F<-with(sdat.md,
                      ifelse(Da1F=="31","30",
                      ifelse(Da1F=="35","30",
                      ifelse(Da1F=="100","100",
                      ifelse(Da1F=="113","100",
                      ifelse(Da1F=="215","100",
                      as.character(Da1F)))))))  
phyR.md<-phyloseq(otu_table(phyR.md),tax_table(phyR.md),sdat.md,phy_tree(phyR.md))
phyR.mdILF<-subset_samples(phyR.md,Sample_Type=="ILF")
sample_data(phyR.mdILF)$Da1F = factor(sample_data(phyR.mdILF)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F

########## Alpha Diversity ##########
dtAD<-phyR.mdILF # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.ILF<-plot_richness(dtAD, x = "Da1F",measures=c("Chao1","InvSimpson")) + 
  geom_boxplot() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_brewer(palette="Set1") 
pDiv.ILF
##### Figure #####
ggsave(grid.draw(pDiv.ILF),filename="Plots/plotAlpha_DaFILF.png",width=8,height=4)

 

phyR.mdInt<-subset_samples(phyR.md,Sample_Type=="Intestinum")
sample_data(phyR.mdInt)$Da1F = factor(sample_data(phyR.mdInt)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F

dtAD<-phyR.mdInt # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.Int<-plot_richness(dtAD, x = "Da1F",measures=c("Chao1","InvSimpson")) + 
  geom_boxplot() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_brewer(palette="Set1") 
pDiv.Int
##### Figure #####
ggsave(grid.draw(pDiv.Int),filename="Plots/plotAlpha_DaFInt.png",width=8,height=4)


