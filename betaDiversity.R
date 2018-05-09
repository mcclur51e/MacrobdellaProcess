############## Ordination/Beta Diversity #################
#library("phyloseq") # should already be loaded
#library("ggplot2") # should already be loaded
library("vegan") # should already be loaded
library("plyr") # should already be loaded

pr10Adult<-subset_samples(physeq10,Age=="A") # keep adult samples (prune Adult)
hi10Adult<-prune_taxa(taxa_sums(prAdult)>.01,pr10Adult) # keep taxa with at least 1% of 1 sample
macro10<-subset_samples(hi10Adult,Taxonomic_ID%in%c("Mdecora","Unk")) # subset containing only Macrobdella samples
macro10N<-subset_samples(macro10,!Replicate%in%c("MN2","MN3","MN4"))
macro10S<-subset_samples(macro10N,!sample_names(macro10N)%in%c("MAa0dF110814EMb.iD691","MAa0dF110814EMc.iD55","MAa0dF110814EMb.uD551","Ma5wF082117EMc.uD695","MAa4dF110514EMa.uD571","MAa4dF110514EMc.uD693","MAa7dF110814EMa.uD555","MAa7dF110814EMc.uD553","MAa4dF110514EMa.iD24","MAa4dF110514EMa.iD18","MAa4dF110514EMb.iD25","MAa4dF110514EMb.uD536","MAa4dF110514EMb.iD24","MAa4dF110514EMb.iD488","MAa0dF110814EMb.iD34","MAa0dF110814EMc.iD450","MAa2dF110314EMc.iD54","Ma061817EMc.iD577","MAa4dF110514EMc.iD91","MAa2dF110314EMc.uD581")) # Remove duplicate samples
macro10daf0<-subset_samples(macro10S,Da1F=="0")
macro10BIU<-subset_samples(macro10daf0,Sample_Type!="Ovary")

########## Alpha Diversity ##########
dtAD<-macro10BIU # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed
#sample_data(addat)$test = factor(sample_data(addat)$test, levels = c("Adult","Hatch","Juvenile","+Water","+Rock","+Adults")) 

pAlpha<-plot_richness(dtAD, x = "Sample_Type",color="Sample_Type",measures=c("Chao1","InvSimpson")) + 
  geom_boxplot() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=45),axis.title.x=element_blank()) +
  scale_fill_brewer(palette="Set1") 
pAlpha
##### Fig #####
ggsave(grid.draw(pAlpha),filename="alpha.png",width=8,height=4)
##### Fig #####
 