### Copied from deontam vignette ###
library("decontam")
ps4<-subset_samples(toolow,Taxonomic_ID%in%c("Mdecora","NegControl")) # subset data so only looking at Mdecora and negative controls

ps2<-subset_samples(ps4,!PostPCRDNA_ng_uL%in%c("Unk")) # Remove any samples where Post PCR [DNA] was not calculated
ps1<-subset_samples(ps2,sample_sums(ps2)>=1) # keep samples with >=1 reads
ps<-prune_taxa(taxa_sums(ps1)>1,ps1) # keep only taxa with >=1 reads

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame

#write.table(sample_data(ps), "psTable.csv", sep=",") # write sample metadata table to csv. Used for quality control

#MAPps$q2<-with(MAPps,ifelse(is.numeric(MAPps$qTest)==TRUE,qTest,abs(as.numeric(as.character(df$qTest)))))
MAPps<-sample_data(ps) # pull map data from phyloseq object
MAPps$q2<-as.numeric(MAPps$PostPCRDNA_ng_uL) # force column to be numeric because R
MAPps$Qubit<-as.numeric(MAPps$Qubit)
MAPps$q2<-prod(10,MAPps$PostPCRDNA_ng_uL)
MAPps$q2<-sum(prod(10,MAPps$PostPCRDNA_ng_uL),1)
ps = merge_phyloseq(ps,MAPps) # return map data to the phyloseq object

df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
write.table(df, "dfTable.csv", sep=",")

ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

### Frequency ###
contamdf.freq <- isContaminant(ps, method="frequency", conc="q2")
table(contamdf.freq$contaminant) # list how many OTUs are contaminants (TRUE/FALSE)
head(contamdf.freq)
head(which(contamdf.freq$contaminant)) # list ranking of contaminant OTUs
ps.contamFreq <- prune_taxa(contamdf.freq$contaminant, ps)
taxa_names(ps.contamFreq) # list contaminant OTUs
tax_table(ps.contamFreq)[,6] # list contaminant genera
plot_frequency(ps, taxa_names(ps)[c(1,33,106,247,264,310)], conc="q2")

#write.table(contamdf.freq, "contamFreqTable.csv", sep=",")
#write.table(tax_table(ps.contamFreq), "contamTaxFreq.csv", sep=",")

### Prevalence ###
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

# Make phyloseq object of presence-absence in negative controls
ps.neg <- prune_samples(sample_data(ps)$Sample_or_Control == "Control", ps)
ps.neg.presence <- transform_sample_counts(ps.neg, function( abund) 1*(abund>0))
# Make phyloseq object of presence-absence in samples
ps.pos <- prune_samples(sample_data(ps)$Sample_or_Control == "Sample", ps)
ps.pos.presence <- transform_sample_counts(ps.pos, function( abund) 1*(abund>0))
# Make data.frame of prevalence in positive and negative samples
df.pres <- data.frame(prevalence.pos=taxa_sums(ps.pos.presence), prevalence.neg=taxa_sums(ps.neg.presence),contam.prev=contamdf.prev$contaminant)
ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos,color=contam.prev)) + 
  geom_point()

ps.contamPrev <- prune_taxa(contamdf.prev$contaminant, ps)
write.table(tax_table(ps.contamPrev), "contamTaxPrev.csv", sep=",")
#write.table(contamdf.freq, "contamFreqTable.csv", sep=",")


### Combined ###
contamdf.com <- isContaminant(ps, method="combined", neg="is.neg", conc="q2")
table(contamdf.com$contaminant)
head(which(contamdf.com$contaminant))
ps.contamCom <- prune_taxa(contamdf.com$contaminant, ps)
write.table(tax_table(ps.contamCom), "contamTaxCom.csv", sep=",")

### Minimum ###
contamdf.min <- isContaminant(ps, method="minimum", neg="is.neg", conc="q2")
table(contamdf.min$contaminant)
head(which(contamdf.min$contaminant))
ps.contamMin <- prune_taxa(contamdf.min$contaminant, ps)
write.table(tax_table(ps.contamMin), "contamTaxMin.csv", sep=",")

### Either ###
contamdf.eit <- isContaminant(ps, method="minimum", neg="is.neg", conc="q2")
table(contamdf.eit$contaminant)
head(which(contamdf.eit$contaminant))
ps.contamEit <- prune_taxa(contamdf.eit$contaminant, ps)
write.table(tax_table(ps.contamEit), "contamTaxEit.csv", sep=",")
