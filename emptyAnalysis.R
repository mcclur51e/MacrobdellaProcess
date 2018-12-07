dir.create("TroubleShoot")

phy.empty<-subset_samples(physeq,Age=="E")
phy.pres<-prune_taxa(taxa_sums(phy.empty)>1,phy.empty) # Remove OTUs with < 100 reads over all samples (1% of 1 sample)

phy.4plates<-subset_samples(phy.empty,!Sample_Plate=="v4.SDSD")
phy.4pru<-prune_taxa(taxa_sums(phy.4plates)>1,phy.4plates) # Rem
plot_bar(phy.4pru,fill="Number")

phy.SDSD<-subset_samples(phy.empty,Sample_Plate=="v4.SDSD")
phy.SDSDpru<-prune_taxa(taxa_sums(phy.SDSD)>1,phy.SDSD) # Rem

mat.SDSDpru <- sort(taxa_sums(phy.SDSDpru), TRUE)[1:50] # Identify 27 most abundant taxa
phy.SDSDhi <- prune_taxa(names(mat.SDSDpru),phy.SDSDpru) # Create a subset of data including only 27 most abundant taxa
pSDSDhi<-plot_bar(phy.SDSDhi,fill="Genus") + scale_fill_manual(values=pairBiome)
ggsave(grid.draw(rbind(ggplotGrob(pSDSDhi), size = "last")), filename="plotEmpty.png", width=12,height=8)
write.table(tax_table(phy.SDSDpru), "table_taxEmpty.csv", sep=",")
write.table(otu_table(phy.SDSDpru), "table_otuEmpty.csv", sep=",")
write.table(sample_data(phy.SDSDpru), "table_mapEmpty.csv", sep=",")



phy.Md<-subset_samples(physeq,Taxonomic_ID=="Mdecora")
phy.MdPres<-prune_taxa(taxa_sums(phy.Md)>1,phy.Md) # Remove OTUs with < 100 reads over all samples (1% of 1 sample)
phy.MdT<-transform_sample_counts(phy.MdPres, function(x) x/sum(x)) # transform raw counts to fraction

mat.MdPres <- sort(taxa_sums(phy.MdT), TRUE)[1:28] # Identify 27 most abundant taxa
phy.MdHi <- prune_taxa(names(mat.MdPres),phy.MdT) # Create a subset of data including only 27 most abundant taxa
plot_bar(phy.MdHi,x="Replicate",fill="Family") + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) + scale_fill_manual(values=pairBiome)

mat.phy <- sort(taxa_sums(physeq), TRUE)[1:28] # Identify 27 most abundant taxa
phy.hi <- prune_taxa(names(mat.phy),physeq) # Create a subset of data including only 27 most abundant taxa
plot_bar(phy.hi,x="Taxonomic_ID",fill="Family") + scale_fill_manual(values=pairBiome)

phy.strain<-subset_samples(physeq,Age=="S")
phy.strainPru<-prune_taxa(taxa_sums(phy.strain)>1,phy.strain) # Rem
mat.strain <- sort(taxa_sums(phy.strain), TRUE)[1:28] # Identify 27 most abundant taxa
phy.strainHi <- prune_taxa(names(mat.strain),phy.strain) # Create a subset of data including only 27 most abundant taxa
plot_bar(phy.strainHi,fill="Number") + scale_fill_manual(values=pairBiome)


########## 2018_0604 ##########
phy.empty<-subset_samples(physeq,Age=="E")
phy.SDSD<-subset_samples(phy.empty,Sample_Plate=="v4.SDSD")

#phy.SDSC<-subset_samples(physeq,Sample_Plate=="v4.SDSC")
#phy.10<-subset_samples(phy.pres,sample_sums(phy.pres)>=10000) # keep samples 


phy.SD704<-subset_samples(physeq,i5_Index_ID%in%c(levels(sample_data(phy.SDSD)$i5_Index_ID)))
#phy.AJE<-subset_samples(phy.SD704,InvertAge%in%c("A","J","E"))
phy.AJE<-phy.SD704
phy.pres<-prune_taxa(taxa_sums(phy.AJE)>10,phy.AJE) # Remove OTUs with < 100 reads over all samples (1% of 1 sample)
phy.hi<-subset_samples(phy.pres,sample_sums(phy.pres)>=1000) # keep samples with >=10000 reads (physeq.minimum 10,000)
phy.Tform.SD704<-transform_sample_counts(phy.hi, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)

mat.SD704 <- sort(taxa_sums(phy.Tform.SD704), TRUE)[1:30] # Identify 50 most abundant taxa
phy.trim.704 <- prune_taxa(names(mat.SD704),phy.Tform.SD704) # Create a subset of data including only most abundant taxa
phy.trim.704<-subset_samples(phy.trim.704,!sample_names(phy.trim.704)=="Empty20180522v4SDSB31")
pEmpty<-plot_bar(phy.trim.704,x="Sample_Well",fill="Genus") + facet_grid(InvertTax~., scales="free_x",space="free") + scale_fill_manual(values=pairBiome)
pEmpty

phy.adult<-subset_samples(phy.trim.704,InvertAge%in%c("E","A"))
pAdult<-plot_bar(phy.adult,x="Sample_Well",fill="Genus") + facet_grid(InvertTax~InvertType, scales="free_x",space="free") + scale_fill_manual(values=pairBiome)
pAdult

dsNMDS<-phy.trim.704 # keep only ILF samples. define data for analysis
mNMDS<-"bray" # define metric for analysis
# calculate distances
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "PCoA", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsGeo<-plot_ordination(dsNMDS, ord, color = "InvertType") + 
  geom_point(size=4,mapping = aes(color=InvertType)) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=brewer.pal(7,"Set1")) 
nmdsGeo





g04<-subset_samples(phy.trim.704,Sample_Well=="G04")

#sample_data(phy.SDSDpru)$i7_Index_ID
#sample_data(phy.SDSDpru)$i5_Index_ID

phy.SD704<-subset_samples(phy.0522,i5_Index_ID%in%c(levels(sample_data(phy.SDSDpru)$i5_Index_ID)))
phy.SD500s<-merge_phyloseq(phy.SD704,phy.SDSD)
phy.Tform.SD500s<-transform_sample_counts(phy.10min, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)

mat.SD500s <- sort(taxa_sums(phy.Tform.SD500s), TRUE)[1:30] # Identify 50 most abundant taxa
phy.trim.500 <- prune_taxa(names(mat.SD500s),phy.Tform.SD500s) # Create a subset of data including only most abundant taxa
plot_bar(phy.trim.500,x="i5_Index_ID",fill="Genus") + facet_grid(i7_Index_ID~., scales="free_x",space="free") + scale_fill_manual(values=pairBiome)



