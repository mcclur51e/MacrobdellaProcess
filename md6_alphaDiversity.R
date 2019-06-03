############## Ordination/Beta Diversity #################
#library("phyloseq") # should already be loaded
#library("ggplot2") # should already be loaded
#library("vegan") # should already be loaded
#library("plyr") # should already be loaded

phyR.adult<-subset_samples(phyR.sin,Age=="A") # keep adult samples (prune Adult)
phyR.md<-subset_samples(phyR.adult,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
phyR.hv<-subset_samples(phyR.adult,sample_names(phyR.adult)%in%c(sample_names(phyT.hv)) & !Da1F%in%c("1","2","4"))
phyR.AD<-merge_phyloseq(phyR.md,phyR.hv)
phyR.adILF<-subset_samples(phyR.AD,Sample_Type%in%c("ILF")) # subset phyR.md to include only ILF samples
sample_data(phyR.adILF)$Da1Fb = factor(sample_data(phyR.adILF)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,76,82,"90+",100,113,215,"Hverbana")) # Reorder Da1F

########## Alpha Diversity ##########
phyR.adILFunfed<-subset_samples(phyR.adILF,!Da2F%in%c("1","2","4","7","8","28")) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed
dtADILF<-subset_samples(phyR.adILFunfed,!Da1F%in%c("1","2","4","7")) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.ILF<-plot_richness(dtADILF, x = "Taxonomic_ID",measures=c("Shannon")) +
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,3)) +
  theme(text=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  scale_fill_manual(values=c("#808080","#000000"))
pDiv.ILF
##### Figure #####

phyR.adInt<-subset_samples(phyR.AD,Sample_Type=="Intestinum")
phyR.adIntunfed<-subset_samples(phyR.adInt,!Da2F%in%c("1","2","4","7","8","28")) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed
dtADInt<-subset_samples(phyR.adIntunfed,!Da1F%in%c("1","2","4","7")) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed

pDiv.Int<-plot_richness(dtADInt, x="Taxonomic_ID", measures=c("Shannon")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,3)) +
  theme(text=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank(),legend.position = "none") +
  scale_fill_manual(values=c("#808080","#000000"))
pDiv.Int
##### Figure #####
ggsave(plot_grid(pDiv.ILF, pDiv.Int, labels = "AUTO"), filename="Plots/plotAlpha_fig2.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 1 to .eps file

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
phyT.mdGI <- subset_samples(phyT.md,Sample_Type%in%c("ILF","Intestinum"))
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
phyT.GI <- subset_samples(phyT.leech,Sample_Type%in%c("ILF","Intestinum"))
phyR.adGI <- subset_samples(phyR.sin,sample_names(phyR.sin)%in%c(sample_names(phyT.GI)))
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




