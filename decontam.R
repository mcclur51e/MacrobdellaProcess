### Copied from deontam vignette ###
library("decontam")
ps4<-subset_samples(physeq,Taxonomic_ID%in%c("Mdecora","NegControl")) # subset data so only looking at Mdecora and negative controls

ps2<-subset_samples(ps4,PostPCRDNA_ng_uL!="Unk") # Remove any samples where Post PCR [DNA] was not calculated
ps1<-subset_samples(ps2,sample_sums(ps2)>=1) # keep samples with >=1 reads
ps<-prune_taxa(taxa_sums(ps1)>1,ps1) # keep only taxa with >=1 reads

df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame

#write.table(sample_data(ps), "psTable.csv", sep=",") # write sample metadata table to csv. Used for quality control

#MAPps$q2<-with(MAPps,ifelse(is.numeric(MAPps$qTest)==TRUE,qTest,abs(as.numeric(as.character(df$qTest)))))
MAPps<-sample_data(ps) # pull map data from phyloseq object
MAPps$q2<-as.numeric(MAPps$qTest) # force column to be numeric because R
MAPps$Qubit<-as.numeric(MAPps$Qubit)
ps = merge_phyloseq(ps,MAPps) # return map data to the phyloseq object

df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
write.table(df, "dfTable.csv", sep=",")

ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

contamdf.freq <- isContaminant(ps, method="frequency", conc="q2")
table(contamdf.freq$contaminant)
head(contamdf.freq)

table(contamdf.freq$contaminant) # list how many OTUs are contaminants (TRUE/FALSE)
head(which(contamdf.freq$contaminant)) # list contaminant OTUs
plot_frequency(ps, taxa_names(ps)[c(114,477,479,147,150,538)], conc="q2")

### Prevalence ###
sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Control"
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
