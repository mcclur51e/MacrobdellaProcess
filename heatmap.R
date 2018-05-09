# Basic heatmap
coreILF<-subset_samples(coreTab,Sample_Type=="ILF") # keep only ILF samples
pCoreILF<-plot_heatmap(coreILF,sample.label="AnimalSource",taxa.label="Genus",title="ILF") # plot heatmap with samples labeled by 'source' and OTUs labeled by 'Genus'
pCoreILF


coreInt<-subset_samples(coreTab,Sample_Type=="Intestinum") # keep only ILF samples
pCoreInt<-plot_heatmap(coreInt,sample.label="AnimalSource",taxa.label="Genus",title="Intestinum") # plot heatmap with samples labeled by 'source' and OTUs labeled by 'Genus'
pCoreInt