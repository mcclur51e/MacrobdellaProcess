############## Ordination/Beta Diversity #################
#library("phyloseq") # should already be loaded
#library("ggplot2") # should already be loaded
#library("vegan") # should already be loaded
#library("plyr") # should already be loaded

phyR.adult<-subset_samples(phyR.out,Age=="A") # keep adult samples (prune Adult)
phyR.md<-subset_samples(phyR.adult,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
phyR.hv<-subset_samples(phyR.adult,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A042317EMe.u","A042317EMe.i","A042317EMf.u","A042317EMh.u","A042117JGbb","A042317EMi","A042317EMj","A042317EMj.b","A042317EMj.i","A042317EMj.u","A042317EMk","A042317EMk.b","A042317EMk.i","A043017EMl","A043017EMm.u","A42hF050317EMb","A050217EMo","A050217EMo.u","A050217EMr.i"))
phyR.AD<-merge_phyloseq(phyR.md,phyR.hv)

sdat.AD<-sample_data(phyR.AD)
sdat.AD$Da1F<-with(sdat.AD,
                      ifelse(Da1F=="31","30",
                      ifelse(Da1F=="35","30",
                      ifelse(Da1F=="100","100",
                      ifelse(Da1F=="113","100",
                      ifelse(Da1F=="215","100",
                      ifelse(Taxonomic_ID=="Hverbana","Hverbana",
                      as.character(Da1F))))))))  
phyR.AD<-phyloseq(otu_table(phyR.AD),tax_table(phyR.AD),sdat.AD,phy_tree(phyR.AD)) # return modified MAP to phyR.md
phyR.adILF<-subset_samples(phyR.AD,Sample_Type=="ILF") # subset phyR.md to include only ILF samples
sample_data(phyR.adILF)$Da1F = factor(sample_data(phyR.adILF)$Da1F, levels = c(0,1,4,7,30,31,35,100,113,215,"Hverbana")) # Reorder Da1F

########## Alpha Diversity ##########
dtAD<-phyR.adILF # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.ILF<-plot_richness(dtAD, x = "Da1F", measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_brewer(palette="Set1") 
pDiv.ILF
##### Figure #####
ggsave(grid.draw(pDiv.ILF),filename="Plots/plotAlpha_DaFILF.png",width=8,height=4)

 

phyR.adInt<-subset_samples(phyR.AD,Sample_Type=="Intestinum")
sample_data(phyR.adInt)$Da1F = factor(sample_data(phyR.adInt)$Da1F, levels = c(0,1,4,7,30,31,35,100,113,215,"Hverbana")) # Reorder Da1F

dtAD<-phyR.adInt # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.Int<-plot_richness(dtAD, x="Da1F", measures=c("Observed","Chao1","ACE","Shannon","Simpson","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_brewer(palette="Set1") 
pDiv.Int
##### Figure #####
ggsave(grid.draw(pDiv.Int),filename="Plots/plotAlpha_DaFInt.png",width=8,height=4)


