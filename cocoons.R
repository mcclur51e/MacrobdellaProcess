coc<-subset_samples(physeqAn,Age=="C")

### Add new column to map data. 'Header' column will be used for pCore figure
MAPcoc<-sample_data(coc) # pull map data from phyloseq object
MAPcoc$CocoonHealth<-with(MAPcoc,
  ifelse(AlbumenColor=="brown","sick",
  ifelse(AlbumenColor=="dk.green","sick",
  ifelse(AlbumenColor=="green","sick",
  ifelse(AlbumenColor=="lt.yellow","healthy",
  ifelse(AlbumenColor=="tan","healthy",
  ifelse(AlbumenColor=="white","healthy",
  as.character(AlbumenColor))))))))
coc = merge_phyloseq(coc,MAPcoc) # return map data to the phyloseq object

cocH<-subset_samples(coc, CocoonHealth=="healthy")
most_abundant_taxa <- sort(taxa_sums(cocH), TRUE)[1:40] # Identify 24 most abundant taxa
pruneCoc <- prune_taxa(names(most_abundant_taxa),cocH)
plot_bar(pruneCoc,fill="Genus",facet_grid=AnimalSource~Var.75) + scale_fill_manual(values=pairBiome) 




