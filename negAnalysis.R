neg<-subset_samples(outlier,Age=="G")
matNeg <- sort(taxa_sums(neg), TRUE)[1:30] # Identify 27 most abundant taxa
pruNeg <- prune_taxa(names(matNeg),neg) # Create a subset of data including only 10 most abundant taxa

otuNeg<-taxa_names(pruNeg)

intersect(otuNeg,coreTot)

pNeg<-plot_bar(pruNeg,fill="Genus") +
#  ylim(0,500) +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome) 
pNeg # print plot

##### Table #####
write.table(tax_table(pruNeg), "negTaxTable.csv", sep=",")
write.table(otu_table(pruNeg),"negOTUTable.csv",sep=",")
##### Table #####

