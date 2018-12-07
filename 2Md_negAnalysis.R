######################################################################################
########## Create plot of negative PCR controls from same dates as samples  ##########
######################################################################################
phy.PCRneg<-subset_samples(phy.sin,Sample_Type%in%c("PCRneg","Reagents")) # keep PCR negative control samples
phy.Md<-subset_samples(phy.sin,Taxonomic_ID%in%c("Mdecora")) # subset containing only Macrobdella samples
phy.Hv<-subset_samples(phy.sin,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A042317EMe.u","A042317EMe.i","A042317EMf.u","A042317EMh.u","A042117JGbb","A042317EMi","A042317EMj","A042317EMj.b","A042317EMj.i","A042317EMj.u","A042317EMk","A042317EMk.b","A042317EMk.i","A043017EMl","A043017EMm.u","A42hF050317EMb","A050217EMo","A050217EMo.u","A050217EMr.i","A0dF091814EMa","A0dF091814EMb"))
ls.PCRnegDates<-intersect(sample_data(phy.Md)$csv_PCRDate,sample_data(phy.PCRneg)$csv_PCRDate) # create list of dates on which PCR negative controls were run with samples of interest
phy.trim.PCRneg<-subset_samples(phy.PCRneg,csv_PCRDate%in%c(ls.PCRnegDates)) # restrict PCR negative samples to only those from ls.PCRnegDates
list.mat.PCRneg <- sort(taxa_sums(phy.trim.PCRneg), TRUE)[1:30] # Identify 30 most abundant taxa
phy.mat.PCRneg <- prune_taxa(names(list.mat.PCRneg),phy.trim.PCRneg) # Create a subset of data including only most abundant taxa
maxNeg <- max(otu_table(phy.mat.PCRneg)) # find the greatest count of any OTU in a single sample

### Plot ###
stack.PCRneg<-plot_bar(phy.mat.PCRneg,x="EmilyID",fill="Genus") +
  ylim(0,500) +
  theme(text=element_text(size=10), axis.title.x=element_blank()) + 
  facet_grid(Sample_Type~., scales="free_x",space="free") +
  scale_fill_manual(values=pairBiome) 
stack.PCRneg # print stacked bar plot
ggsave(grid.draw(rbind(ggplotGrob(stack.PCRneg), size = "last")), filename="NegControls/plot_stack_negPCR.png", width=12,height=8)
### Table ###
phy.taxaNeg<-subset_taxa(phy.trim.PCRneg,taxa_sums(phy.trim.PCRneg)>1) # keep only OTUs with >1 read total
write.table(tax_table(phy.taxaNeg), "NegControls/table_tax_negPCR.csv", sep=",")
write.table(otu_table(phy.taxaNeg),"NegControls/table_OTU_negPCR.csv",sep=",")

########## Identify 'contaminant' OTUs  ##########
##### Copied from deontam vignette #####
library("decontam")
dir.create("NegControls/Decontam",showWarnings = FALSE) # create folder for data output(s)
### Trim unuseable samples ###
phy.MdnegRaw<-merge_phyloseq(phy.Md,phy.trim.PCRneg,phy.Hv) # combine phyloseq objects with Macrobdella data and negative controls from those MiSeq runs (phyloseq.Macrobdella & negatives)
ps2<-subset_samples(phy.MdnegRaw,!PostPCRDNA_ng_uL%in%c("Unk")) # Remove any samples where Post PCR [DNA] was not calculated
ps1<-subset_samples(ps2,sample_sums(ps2)>=1) # keep samples with >=1 reads
phy.Mdneg<-prune_taxa(taxa_sums(ps1)>1,ps1) # keep only taxa with >=1 reads
### Map ###
map.Mdneg<-sample_data(phy.Mdneg) # pull map data from phyloseq object
map.Mdneg$Conc<-with(map.Mdneg,
                 ifelse(PostPCRDNA_ng_uL%in%c("NA","zero","Unk",""), "0", 
                        as.character(PostPCRDNA_ng_uL)))
map.Mdneg$c <- as.numeric(as.character(map.Mdneg$Conc)) * 10 + 1 # ridiculous math to make decontam work
map.Mdneg$q2 <- as.numeric(as.character(map.Mdneg$c))
phy.Mdneg = merge_phyloseq(phy.Mdneg,map.Mdneg) # return map data to the phyloseq object
write.table(map.Mdneg, "NegControls/Decontam/table_decontam_MAP.csv", sep=",")
### Table ###
df.Mdneg <- as.data.frame(sample_data(phy.Mdneg)) # Put sample_data into a ggplot-friendly data.frame
df.Mdneg$LibrarySize <- sample_sums(phy.Mdneg)
df.Mdneg <- df.Mdneg[order(df.Mdneg$LibrarySize),]
df.Mdneg$Index <- seq(nrow(df.Mdneg))
write.table(df.Mdneg, "NegControls/Decontam/table_decontam_df.csv", sep=",")
### Plot ###
ggplot(data=df.Mdneg, aes(x=Index, y=LibrarySize, color=Sample_Type)) + geom_point()

### FREQUENCY ###
contamdf.freq <- isContaminant(phy.Mdneg, method="frequency", conc="q2")
table(contamdf.freq$contaminant) # list how many OTUs are contaminants (TRUE/FALSE)
phy.contamFreq <- prune_taxa(contamdf.freq$contaminant, phy.Mdneg)
list.mat.contamFreq <- sort(taxa_sums(phy.contamFreq), TRUE)[1:16] # Identify 16 most abundant contaminating taxa
phy.mat.PCRneg <- prune_taxa(names(list.mat.PCRneg),phy.trim.PCRneg) # Create a subset of data including only most abundant taxa
### Plot ###
plot_frequency(phy.Mdneg, names(list.mat.contamFreq), conc="q2") # use this to only plot specific contaminants
### Table ###
write.table(contamdf.freq, "NegControls/Decontam/table_decontam_contamStatsFreq.csv", row.names=TRUE, sep=",")
write.table(tax_table(phy.contamFreq), "NegControls/Decontam/table_decontam_contamTaxFreq.csv", sep=",")

### PREVALENCE ###
sample_data(phy.Mdneg)$is.neg <- sample_data(phy.Mdneg)$Sample_or_Control=="NegControl"
contamdf.prev <- isContaminant(phy.Mdneg, method="prevalence", neg="is.neg")

# Make phyloseq object of presence-absence in negative controls
phy.neg <- prune_samples(sample_data(phy.Mdneg)$Sample_or_Control=="NegControl", phy.Mdneg)
phy.neg.presence <- transform_sample_counts(phy.neg, function( abund) 1*(abund>0))
# Make phyloseq object of presence-absence in samples
phy.pos <- prune_samples(sample_data(phy.Mdneg)$Sample_or_Control%in%c("Sample","PosControl"), phy.Mdneg)
phy.pos.presence <- transform_sample_counts(phy.pos, function( abund) 1*(abund>0))
# Make data.frame of prevalence in positive and negative samples
df.pres <- data.frame(prevalence.pos=taxa_sums(phy.pos.presence), prevalence.neg=taxa_sums(phy.neg.presence),contam.prev=contamdf.prev$contaminant)
ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos,color=contam.prev)) + 
  geom_point()

phy.contamPrev <- prune_taxa(contamdf.prev$contaminant, phy.Mdneg)
### Table ###
write.table(tax_table(phy.contamPrev), "NegControls/Decontam/table_decontam_contamTaxPrev.csv", sep=",")

### COMBINED ###
contamdf.com <- isContaminant(phy.Mdneg, method="combined", neg="is.neg", conc="q2")
phy.contamCom <- prune_taxa(contamdf.com$contaminant, phy.Mdneg)
### Table ###
write.table(tax_table(phy.contamCom), "NegControls/Decontam/table_decontam_contamTaxCom.csv", sep=",")

### MINIMUM ###
contamdf.min <- isContaminant(phy.Mdneg, method="minimum", neg="is.neg", conc="q2")
phy.contamMin <- prune_taxa(contamdf.min$contaminant, phy.Mdneg)
### Table ###
write.table(tax_table(phy.contamMin), "NegControls/Decontam/table_decontam_contamTaxMin.csv", sep=",")

phyR.nc <- prune_taxa(taxa_names(phy.contamFreq), phy.Mdneg) # (physeq no contaminants)
phyR.allc<- prune_taxa(taxa_names(phy.contamFreq), phy.Mdneg) # (physeq all contaminants)


##############################
### Temporary Testing Area ###
##############################
maxA <- max(otu_table(phy.allc)) # find the greatest count of any OTU in a single sample
list.mat.A <- sort(taxa_sums(phy.allc), TRUE)[1:35] # Identify most abundant contaminating taxa
phy.mat.A <- prune_taxa(names(list.mat.A),phy.allc) # Create a subset of data including only most abundant taxa

plot_bar(phy.mat.A,x="EmilyID",fill="Genus2") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) 

phy.B<-subset_samples(phy.mat.A,sample_sums(phy.allc)>maxA) # keep samples with > maxA reads (physeq.minimum 10,000)


plot_bar(phy.B,x="EmilyID",fill="Genus2") +
#  ylim(0,1000) +
  theme(text=element_text(size=10), axis.title.x=element_blank()) 

phy.decontam<-subset_samples(phy.sin,!sample_names(phy.sin)%in%c(sample_names(phy.B)))

write.table(tax_table(phy.decontam), "NegControls/Decontam/table_decontam_Tax.csv", sep=",")

