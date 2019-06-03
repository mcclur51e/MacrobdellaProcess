library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. #
library("structSSI")
#library("ggplot2") # should already be loaded
#theme_set(theme_bw())

##### Shifts in abundance v. species #####
phyT.mdILF<-subset_samples(phyT.md,Sample_Type=="ILF" & Da1Fb%in%c(0,2))
phyR.dsq<-subset_samples(phyR.sin,sample_names(phyR.sin)%in%c(sample_names(phyT.mdILF)))
dsq.start = phyloseq_to_deseq2(phyR.dsq, ~ Da1Fb)

geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 1, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP")]
ls.02<-levels(posigtab$Number)

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=Genus2, x=log2FoldChange, color=Order)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(9,"Set1")) 
###
phyT.mdILF<-subset_samples(phyT.md,Sample_Type=="ILF" & Da1Fb%in%c(0,4))
phyR.dsq<-subset_samples(phyR.sin,sample_names(phyR.sin)%in%c(sample_names(phyT.mdILF)))
dsq.start = phyloseq_to_deseq2(phyR.dsq, ~ Da1Fb)

geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 1, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP")]
ls.04<-levels(posigtab$Number)

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=Genus2, x=log2FoldChange, color=Order)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(9,"Set1")) 


ls.change<-union(ls.02,ls.04)
phyT.BWchange<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA),Sample_Type=="ILF" & Da2F%in%c("none","Unk","0"))
phyT.BWchange<-prune_taxa(c(ls.change),phyT.BWchange) # keep only core ILF and intestinum taxa (phyloseq Transformed . Days after Feeding ILF pruned)

sample_data(phyT.BWchange)$Da1Fb = factor(sample_data(phyT.BWchange)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+")) # Reorder Da1Fb

dt.BWmerge <- psmelt(phyT.BWchange) # create data.table from phyloseq object (data table . Box+Whisker merged)
pBW.DaF <- ggplot(dt.BWmerge) +
  stat_boxplot(data=dt.BWmerge, aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  scale_y_log10(limits=c(.001,1)) + 
  facet_grid(Taxonomic_ID+Sample_Type~Order+Genus2, scales="free_x",space="free") +
  theme_bw() +
  xlab("Days after feeding") +
  theme(text=element_text(size=10),strip.text=element_text(size=6), axis.text.x=element_text(angle=90,hjust=1), legend.position="none") +
  scale_fill_manual(values=brewer.pal(8,"Greys")) 
pBW.DaF

##### Shifts in abundance v. sample type in Macrobdella decora #####
phyT.mdMain<-subset_samples(phyT.md,!AnimalSource%in%c("MtSnowVT","CarogaNY"))
phyT.mdGI<-subset_samples(phyT.mdMain,Sample_Type%in%c("Intestinum"))
phyR.dsq<-subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.mdGI)))
phyT.dsq<-subset_samples(phyT.Tform,sample_names(phyR.out)%in%c(sample_names(phyT.mdGI)))
dsq.mdGI = phyloseq_to_deseq2(phyR.dsq, ~ AnimalSource)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dsq.mdGI), 1, gm_mean)
esf.mdGI = estimateSizeFactors(dsq.mdGI, geoMeans = geoMeans)
dsq.mdGI = DESeq(esf.mdGI, fitType="local")

res = results(dsq.mdGI)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"), as(as.matrix(taxa_sums(phyT.dsq))[rownames(sigtab),],"numeric"))
names(sigtab)[21]<-"count"
head(sigtab)

posigtab = sigtab[(abs(sigtab[, "log2FoldChange"]) > 2 & sigtab[,"count"]>.001), ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP","Blastn","Identity","count")]

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=Genus2, x=log2FoldChange, color=Order)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(9,"Set1")) 

phyT.mdGIdelta<-subset_taxa(phyT.mdGI,taxa_names(phyT.mdGI)%in%c(rownames(posigtab)))
phyT.mdGIdelta<-subset_taxa(phyT.mdGIdelta,Genus2!="Aeromonas2")
### Stacked Bar Plot
pBar.mdGIdelta<-plot_bar(phyT.mdGIdelta,fill="Family") + 
  facet_grid(Sample_Type~Genus2, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(8,"Set1"))
pBar.mdGIdelta
ggsave(grid.draw(rbind(ggplotGrob(pBar.mdGIdelta), size = "last")), filename="Plots/plotBarStack_mdGIdelta.eps", device="eps", width=12,height=8)

max(sample_sums(subset_samples(phyT.mdGIdelta,Sample_Type=="ILF")))
min(sample_sums(subset_samples(phyT.mdGIdelta,Sample_Type=="ILF")))

max(sample_sums(subset_samples(phyT.mdGIdelta,Sample_Type=="Intestinum")))

##### Violin Plot #####
dt.mdGIdelta<-data.table(psmelt(phyT.mdGIdelta))
pVio.mdGIdelta<-ggplot(dt.mdGIdelta, aes(x=Sample_Type,y=Abundance, fill=Sample_Type)) + 
  geom_violin(aes(fill=Sample_Type),alpha=0.7) +
  geom_jitter(width = 0.2) +
  ggtitle("OTUs that drive change between ILF and intestinum microbiota in M.decora") + 
  theme(text=element_text(size=10),axis.title.x=element_blank(),legend.position="none") +
  facet_grid(AnimalSource~Genus2, scales="free_x",space="free") + 
  scale_y_log10() + 
  scale_fill_manual(values=brewer.pal(6,"Set1")) +
  scale_color_brewer(type = 'qual', palette = "Set1")
pVio.mdGIdelta
ggsave(grid.draw(rbind(ggplotGrob(pVio.mdGIdelta), size = "last")), filename="Plots/plotViolin_mdGIdelta.eps", device="eps", width=12,height=8)






##### Shifts in abundance v. 0Da1F #####
phyT.ilf04<-subset_samples(phyT.start,Sample_Type=="ILF" & Taxonomic_ID=="Mdecora" & Da1Fb%in%c(0,4))
phyR.dsq<-subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.ilf04)))
dsq.start = phyloseq_to_deseq2(phyR.dsq, ~ Da1Fb)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 2, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP")]

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=Number, x=log2FoldChange, color=Genus2)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(9,"Set1")) 





###### Bladder vs non
##### Shifts in abundance v. sample type in Macrobdella decora #####
sample_data(phyR.out)$System<-with(sample_data(phyR.out),
                                     ifelse(Sample_Type=="Bladder","Bladder",
                                     ifelse(Sample_Type=="ILF","GI",
                                     ifelse(Sample_Type=="Intestinum","GI",
                                     as.factor(Sample_Type)))))
phyT.dsq<-subset_samples(phyT.md, Sample_Type!="Ovary")
phyR.dsq<-subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.dsq)))
dsq.mdGI = phyloseq_to_deseq2(phyR.dsq, ~ System)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(dsq.mdGI), 1, gm_mean)
esf.mdGI = estimateSizeFactors(dsq.mdGI, geoMeans = geoMeans)
dsq.mdGI = DESeq(esf.mdGI, fitType="local")

res = results(dsq.mdGI)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"), as(as.matrix(taxa_sums(phyT.dsq))[rownames(sigtab),],"numeric"))
names(sigtab)[19]<-"count"
head(sigtab)

posigtab = sigtab[(abs(sigtab[, "log2FoldChange"]) > 2 & sigtab[,"count"]>.01), ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP","Blastn","Identity","count")]

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=RDP, x=log2FoldChange, color=Class)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(10,"Paired")) 

phyT.mdBladGI<-subset_taxa(phyT.dsq,taxa_names(phyT.dsq)%in%c(rownames(posigtab)))
### Stacked Bar Plot
pBar.mdBladGI<-plot_bar(phyT.mdBladGI,fill="Genus2") + 
  facet_grid(Sample_Type~Genus2, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(9,"Set1"))
pBar.mdBladGI
ggsave(grid.draw(rbind(ggplotGrob(pBar.mdBladGI), size = "last")), filename="Plots/plotBarStack_mdGIdelta.eps", device="eps", width=12,height=8)
### Box & Whisker Plot
dt.mdBladGI<-data.table(psmelt(phyT.mdBladGI))
pVio.mdBladGI<-ggplot(dt.mdBladGI, aes(x=System,y=Abundance, fill=System)) + 
  geom_violin(aes(fill=System),alpha=0.7) +
  geom_jitter(width = 0.2) +
  ggtitle("OTUs that drive change between GI and bladder microbiota in Macrobdella decora") + 
  theme(text=element_text(size=10),axis.title.x=element_blank(),legend.position="none") +
  facet_grid(.~Genus2, scales="free_x",space="free") + 
  scale_y_log10() + 
  scale_fill_manual(values=brewer.pal(6,"Set1")) +
  scale_color_brewer(type = 'qual', palette = "Set1")
pVio.mdBladGI
#pBar.mdBladGI
ggsave(grid.draw(rbind(ggplotGrob(pVio.mdBladGI), size = "last")), filename="Plots/plotViolin_mdGIdelta.eps", device="eps", width=12,height=8)









##### Shifts in abundance v. sample type in Hirudo #####
phyT.dsq<-subset_samples(phyT.base, Taxonomic_ID=="Mdecora")
phyT.dsq<-subset_samples(phyT.dsq,Sample_Type!="Bladder")

ls.hvBlad<-taxa_names(prune_taxa(taxa_sums(subset_samples(phyT.dsq,Sample_Type=="Intestinum")) > .001, phyT.dsq))
ls.hvGI<-taxa_names(prune_taxa(taxa_sums(subset_samples(phyT.dsq,Sample_Type=="ILF")) > .001, phyT.dsq))
ls.hvBladUnq<-setdiff(ls.hvBlad,ls.hvGI)
ls.hvGIUnq<-setdiff(ls.hvGI,ls.hvBlad)

phyR.dsq<-subset_samples(phyR.out,sample_names(phyR.out)%in%c(sample_names(phyT.dsq)))
dsq.hv = phyloseq_to_deseq2(phyR.dsq, ~ Sample_Type)
geoMeans = apply(counts(dsq.hv), 1, gm_mean)
esf.hv = estimateSizeFactors(dsq.hv, geoMeans = geoMeans)
dsq.hv = DESeq(esf.hv, fitType="local")

res = results(dsq.hv)[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"), as(as.matrix(taxa_sums(phyT.dsq))[rownames(sigtab),],"numeric"))
names(sigtab)[19]<-"count"
#head(sigtab)

posigtab = sigtab[(abs(sigtab[, "log2FoldChange"]) > 1), ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family","Genus","Genus2","Number","RDP","BLASTn","percentID")]

sigtabgen = subset(sigtab, !is.na(Genus2))
# Order order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x)) 
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x)) 
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Number, function(x) max(x)) 
x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Number), levels=names(x)) 
ggplot(posigtab, aes(y=Number, x=log2FoldChange, color=Class)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=brewer.pal(9,"Set1")) 

write.table(posigtab[,c("baseMean","log2FoldChange","RDP")],"table_foldChange.csv",sep=",") # Make list of core OTUs for Macrobdella



