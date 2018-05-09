########## Separate hatchlings and plot to evaluate if enough samples have been sequenced ##########

# Control
hatchT<-transform_sample_counts(hatch, function(x) x/sum(x)) # transform raw counts to fraction
most_abundant_taxa <- sort(taxa_sums(hatchT), TRUE)[1:15] # Identify 15 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),hatchT) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus") + scale_fill_manual(values=pairBiome)



phyYo<-merge_phyloseq(hatch,juvCon) # merge hatchling and untreated-juvenile phyloseq objects (physeq young)
phyYoT<-transform_sample_counts(phyYo, function(x) x/sum(x)) # transform raw counts to fraction

########## ILF/AnimalSource NMDS ##########
dsNMDS<-subset_samples(phyYoT,Sample_Type=="ILF") # keep only ILF samples. define data for analysis
mNMDS<-"unifrac" # define metric for analysis
# calculate distances
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "PCoA", distance = distOrd) # calculate ordination
# PCoA. May want to change color= , shape= , alpha= , and size=
nmdsYo<-plot_ordination(dsNMDS, ord, color = "Da1F") + 
  geom_point(size=6,mapping = aes(color=Da1F)) +
  ggtitle(c(mNMDS)) + 
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pairBiome)
nmdsYo