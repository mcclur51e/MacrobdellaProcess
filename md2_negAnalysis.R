######################################################################################
########## Create plot of negative PCR controls from same dates as samples  ##########
######################################################################################
phyR.PCRneg<-subset_samples(phyR.sin,Sample_Type%in%c("PCRneg","Reagents")) # keep PCR negative control samples
phyR.Md<-subset_samples(phyR.sin,Taxonomic_ID%in%c("Mdecora")) # subset containing only Macrobdella samples
phyR.Hv<-subset_samples(phyR.sin,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A042317EMe.u","A042317EMe.i","A042317EMf.u","A042317EMh.u","A042117JGbb","A042317EMi","A042317EMj","A042317EMj.b","A042317EMj.i","A042317EMj.u","A042317EMk","A042317EMk.b","A042317EMk.i","A043017EMl","A043017EMm.u","A42hF050317EMb","A050217EMo","A050217EMo.u","A050217EMr.i","A0dF091814EMa","A0dF091814EMb"))
ls.PCRnegDates<-intersect(sample_data(phyR.Md)$csv_PCRDate,sample_data(phyR.PCRneg)$csv_PCRDate) # create list of dates on which PCR negative controls were run with samples of interest
phyR.trim.PCRneg<-subset_samples(phyR.PCRneg,csv_PCRDate%in%c(ls.PCRnegDates)) # restrict PCR negative samples to only those from ls.PCRnegDates
phyR.mat.PCRneg <- prune_taxa(names(sort(taxa_sums(phyR.trim.PCRneg),TRUE)[1:30]),phyR.trim.PCRneg) # Create a subset of data including only (30) most abundant taxa

### Plot ###
pbar.PCRneg<-plot_bar(phyR.mat.PCRneg,x="EmilyID",fill="Genus") +
  ylim(0,500) +
  theme(text=element_text(size=10), axis.title.x=element_blank()) + 
  facet_grid(Sample_Type~., scales="free_x",space="free") +
  scale_fill_manual(values=pairBiome) 
pbar.PCRneg # print stacked bar plot
ggsave(grid.draw(rbind(ggplotGrob(pbar.PCRneg), size = "last")), filename="NegControls/plot_stack_negPCR.png", width=12,height=8)
### Table ###
phyR.taxaNeg<-subset_taxa(phyR.trim.PCRneg,taxa_sums(phyR.trim.PCRneg)>1) # keep only OTUs with >1 read total
write.table(tax_table(phyR.taxaNeg), "NegControls/table_tax_negPCR.csv", sep=",")
write.table(otu_table(phyR.taxaNeg),"NegControls/table_OTU_negPCR.csv",sep=",")
var.maxNeg <- max(otu_table(phyR.taxaNeg)) # find the greatest count of any OTU in a single sample
var.meanNeg <- mean(otu_table(phyR.taxaNeg)) # find the median count of any OTU in a single sample

########## Identify 'contaminant' OTUs  ##########
##### Copied (and slightly modified) from deontam vignette #####
#library("decontam") # should already be loaded
dir.create("NegControls/Decontam",showWarnings = FALSE) # create folder for data output(s)
### Trim unuseable samples ###
phyR.MdnegRaw<-merge_phyloseq(phyR.Md,phyR.trim.PCRneg,phyR.Hv) # combine phyloseq objects with Macrobdella data and negative controls from those MiSeq runs (phyloseq.Macrobdella & negatives)
phyR.MdnegConc<-subset_samples(phyR.MdnegRaw,!PostPCRDNA_ng_uL%in%c("Unk")) # Remove any samples where Post PCR [DNA] was not calculated
phyR.MdnegCount<-subset_samples(phyR.MdnegConc,sample_sums(phyR.MdnegConc)>=1) # keep samples with >=1 reads
phyR.Mdneg<-prune_taxa(taxa_sums(phyR.MdnegCount)>1,phyR.MdnegCount) # keep only taxa with >=1 reads
### Map ###
map.Mdneg<-sample_data(phyR.Mdneg) # pull map data from phyloseq object
map.Mdneg$Conc<-with(map.Mdneg,
                 ifelse(PostPCRDNA_ng_uL%in%c("NA","zero","Unk",""), "0", 
                        as.character(PostPCRDNA_ng_uL)))
map.Mdneg$c <- as.numeric(as.character(map.Mdneg$Conc)) * 10 + 1 # ridiculous math to make decontam work
map.Mdneg$q2 <- as.numeric(as.character(map.Mdneg$c))
phyR.Mdneg = merge_phyloseq(phyR.Mdneg,map.Mdneg) # return map data to the phyloseq object
### Table ###
df.Mdneg <- as.data.frame(sample_data(phyR.Mdneg)) # Put sample_data into a ggplot-friendly data.frame
df.Mdneg$LibrarySize <- sample_sums(phyR.Mdneg)
df.Mdneg <- df.Mdneg[order(df.Mdneg$LibrarySize),]
df.Mdneg$Index <- seq(nrow(df.Mdneg))
### Plot ###
pScat.libSize<-ggplot(data=df.Mdneg, aes(x=Index, y=LibrarySize, color=Sample_Type)) + 
  geom_point() +
  theme(text=element_text(size=10), axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_y_log10() + 
  scale_fill_manual(values=brewer.pal(6,"Set1"))
pScat.libSize

### FREQUENCY ###
contamdf.freq <- isContaminant(phyR.Mdneg, method="frequency", conc="q2")
table(contamdf.freq$contaminant) # list how many OTUs are contaminants (TRUE/FALSE)
phyR.contamFreq <- prune_taxa(contamdf.freq$contaminant, phyR.Mdneg)
#ls.contamFreq <- contamdf.freq$contaminant==TRUE
ls.mat.contamFreq <- sort(taxa_sums(phyR.contamFreq), TRUE)[1:16] # Identify 16 most abundant contaminating taxa
### Plot ###
plot_frequency(phyR.Mdneg, names(ls.mat.contamFreq), conc="q2") # use this to only plot specific contaminants
### Table ###
write.table(contamdf.freq, "NegControls/Decontam/table_decontam_contamStatsFrequency.csv", row.names=TRUE, sep=",")

### PREVALENCE ###
sample_data(phyR.Mdneg)$is.neg <- sample_data(phyR.Mdneg)$Sample_or_Control=="NegControl"
contamdf.prev <- isContaminant(phyR.Mdneg, method="prevalence", neg="is.neg")
#ls.contamPrev <- contamdf.prev$contaminant

# Make phyloseq object of presence-absence in negative controls
phyR.neg <- prune_samples(sample_data(phyR.Mdneg)$Sample_or_Control=="NegControl", phyR.Mdneg)
phyR.neg.presence <- transform_sample_counts(phyR.neg, function( abund) 1*(abund>0))
# Make phyloseq object of presence-absence in samples
phyR.pos <- prune_samples(sample_data(phyR.Mdneg)$Sample_or_Control%in%c("Sample","PosControl"), phyR.Mdneg)
phyR.pos.presence <- transform_sample_counts(phyR.pos, function( abund) 1*(abund>0))
# Make data.frame of prevalence in positive and negative samples
df.pres <- data.frame(prevPos=taxa_sums(phyR.pos.presence), prevNeg=taxa_sums(phyR.neg.presence),contaminant=contamdf.prev$contaminant)
pScat.contamPrev<-ggplot(data=df.pres, aes(x=prevNeg, y=prevPos,color=contaminant)) + 
  geom_point() +
  ggtitle("OTU prevalence in positive and negative samples") + 
  theme(text=element_text(size=10), axis.text.x=element_text(angle=90)) +
  scale_fill_manual(values=brewer.pal(6,"Set1"))
pScat.contamPrev
write.table(contamdf.prev, "NegControls/Decontam/table_decontam_contamStatsPrevalence.csv", row.names=TRUE, sep=",")

### COMBINED ###
contamdf.com <- isContaminant(phyR.Mdneg, method="combined", neg="is.neg", conc="q2")
#ls.contamCom <- contamdf.com$contaminant
write.table(contamdf.com, "NegControls/Decontam/table_decontam_contamStatsCombined.csv", row.names=TRUE, sep=",")

### MINIMUM ###
contamdf.min <- isContaminant(phyR.Mdneg, method="minimum", neg="is.neg", conc="q2")
#ls.contamMin <- contamdf.min$contaminant
write.table(contamdf.min, "NegControls/Decontam/table_decontam_contamStatsMinimum.csv", row.names=TRUE, sep=",")

phyR.allc <- prune_taxa(taxa_names(phyR.contamFreq), phyR.Mdneg) # (physeq all contaminants)
phyR.nc <- subset_taxa(phyR.Mdneg,!taxa_names(phyR.Mdneg)%in%c(taxa_names(phyR.contamFreq)))

##### Plot of contaminants in samples #####
phyT.Mdneg<-transform_sample_counts(phyR.Mdneg, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)
phyT.allc<- prune_taxa(taxa_names(phyR.contamFreq), phyT.Mdneg) # (physeq all contaminants)
phyG.allc<-tax_glom(phyT.allc,taxrank="Genus")
phyG.cont<-subset_samples(phyG.allc,!Age%in%c("M","G"))

pbar.Contam<-plot_bar(phyG.cont,x="EmilyID",fill="Genus") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pbar.Contam # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbar.Contam), size = "last")), filename="Plots/plotBarStack_Contaminants.png", width=16,height=8)





##############################
### Temporary Testing Area ###
##############################
