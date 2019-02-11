############## Ordination/Beta Diversity #################
#library("phyloseq") # should already be loaded
#library("ggplot2") # should already be loaded
#library("vegan") # should already be loaded
#library("plyr") # should already be loaded
color.grey<-c("#808080","#000000")

phyR.adult<-subset_samples(phyR.out,Age=="A") # keep adult samples (prune Adult)
phyR.md<-subset_samples(phyR.adult,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
phyR.hv<-subset_samples(phyR.adult,sample_names(phyR.adult)%in%c(sample_names(phyT.Hv)))
phyR.AD<-merge_phyloseq(phyR.md,phyR.hv)
phyR.adILF<-subset_samples(phyR.AD,Sample_Type=="ILF") # subset phyR.md to include only ILF samples
sample_data(phyR.adILF)$Da1Fb = factor(sample_data(phyR.adILF)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215,"Hverbana")) # Reorder Da1F

########## Alpha Diversity ##########
dtAD<-phyR.adILF # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.ILF<-plot_richness(dtAD, x = "Da1Fb",measures=c("Observed","Shannon","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  ggtitle("ILF") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_manual(values=color.grey)
pDiv.ILF
##### Figure #####
ggsave(grid.draw(pDiv.ILF),filename="Plots/plotAlpha_DaFILF.eps", device="eps",width=8,height=4)

 
phyR.adInt<-subset_samples(phyR.AD,Sample_Type=="Intestinum")
sample_data(phyR.adInt)$Da1Fb = factor(sample_data(phyR.adInt)$Da1F, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215,"Hverbana")) # Reorder Da1F
dtAD<-phyR.adInt # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.Int<-plot_richness(dtAD, x="Da1Fb", measures=c("Observed","Shannon","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  ggtitle("Intestinum") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_manual(values=color.grey) 
pDiv.Int
##### Figure #####
ggsave(grid.draw(pDiv.Int),filename="Plots/plotAlpha_DaFInt.eps", device="eps",width=8,height=4)

##### ILF vs Intestinum #####
phyT.hvGI <- subset_samples(phyT.Hv,Sample_Type%in%c("ILF","Intestinum"))
phyR.adhvGI <- subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.hvGI)))
dtAD<-phyR.adhvGI # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.hvGI<-plot_richness(dtAD, x="Sample_Type", measures=c("Observed","Shannon","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  ggtitle("Hirudo verbana GI") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_manual(values=color.grey) 
pDiv.hvGI
##### Figure #####
ggsave(grid.draw(pDiv.hvGI),filename="Plots/plotAlpha_hvGI.eps", device="eps",width=8,height=4)

##### ILF vs Intestinum #####
phyT.mdGI <- subset_samples(phyT.Md,Sample_Type%in%c("ILF","Intestinum"))
phyR.admdGI <- subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.mdGI)))
dtAD<-phyR.admdGI # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.mdGI<-plot_richness(dtAD, x="Sample_Type", measures=c("Observed","Shannon","InvSimpson")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  ggtitle("Macrobdella decora GI") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_manual(values=color.grey) 
pDiv.mdGI
##### Figure #####
ggsave(grid.draw(pDiv.mdGI),filename="Plots/plotAlpha_mdGI.eps", device="eps",width=8,height=4)


##### ILF vs Intestinum #####
phyT.GI <- subset_samples(phyT.start,Sample_Type%in%c("ILF","Intestinum"))
phyR.adGI <- subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.GI)))
dtAD<-phyR.adGI # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.GI<-plot_richness(dtAD, x="Sample_Type", color="Taxonomic_ID", measures=c("Observed","Shannon","InvSimpson")) + 
  geom_violin(aes(fill=Taxonomic_ID),color="black",alpha=0.7) +
  geom_jitter(width = 0.2) +
  ggtitle("GI") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(6,"Set1")) +
  scale_color_brewer(type = 'qual', palette = "Set1")
pDiv.GI
##### Figure #####
ggsave(grid.draw(pDiv.GI),filename="Plots/plotViolin_GI.eps", device="eps",width=8,height=4)




