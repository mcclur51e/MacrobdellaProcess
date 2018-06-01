neg<-subset_samples(outlier,Age=="G")
negP<-subset_samples(neg,!sample_names(neg)=="NegEMbeatEM122116.431")
matNeg <- sort(taxa_sums(negP), TRUE)[1:30] # Identify 27 most abundant taxa
pruNeg <- prune_taxa(names(matNeg),negP) # Create a subset of data including only 10 most abundant taxa

otuNeg<-taxa_names(pruNeg)

intersect(otuNeg,coreTot)

pNeg<-plot_bar(pruNeg,fill="Genus") +
#  ylim(0,500) +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome) 
pNeg # print plot
ggsave(grid.draw(rbind(ggplotGrob(pNeg), size = "last")), filename="plotNeg.png", width=12,height=8)

##### Table #####
write.table(tax_table(pruNeg), "table_taxNeg.csv", sep=",")
write.table(otu_table(pruNeg),"table_OTUneg.csv",sep=",")
##### Table #####

