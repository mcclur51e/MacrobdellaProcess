#############################################
##### Practice for evaluating core taxa #####
#############################################

ls.outliers<-as.character(read.csv("DataFiles/list_outliers.csv",header=FALSE,sep=",")$V1) # Import list of samples determined to be outliers. Last evaluated 2018-11-22. (list.outliers)
phy.red<-subset_samples(phy.sin, !sample_names(phy.sin)%in%c(ls.outliers))# Remove outliers from phy.sin (physeq.reduced)

phy.Adult<-subset_samples(phy.red,Age%in%c("A")) # Keep only adult samples (physeq.Adult)
phy.ILF<-subset_samples(phy.Adult,Sample_Type%in%c("ILF","Intestinum"))
phy.Mdeco<-subset_samples(phy.ILF,Taxonomic_ID=="Mdecora") # Keep only Mdecora samples (physeq.Mdecora)
phyT.Mdeco<-transform_sample_counts(phy.Mdeco, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform.Mdecora)
mat.Mdeco<-sort(taxa_sums(phyT.Mdeco), TRUE)[1:50] # Identiy top 28 OTUs from phy.Mdeco (most abundant taxa.Mdecora)
phyT.pruMdeco<-prune_taxa(names(mat.Mdeco),phyT.Mdeco) # prune phy.Mdeco to include only the most abundant taxa as defined in mat.Mdeco (physeq.prunedMdecora)

pBar.pruMdeco<-plot_bar(phyT.pruMdeco,x="Source_ID",fill="Genus2") + 
  facet_grid(Sample_Type~Da1F, scales="free_x",space="free") + 
  scale_fill_manual(values=pairBiome)
pBar.pruMdeco

phyT.noCT<-subset_samples(phyT.pruMdeco,AnimalSource!="Wlot") # keep only samples from GrotonMA
macDaF<-phyT.Gr # define to be put into NMDS analysis

fm.noCT = fast_melt(phyT.noCT) # (fast melt Vermont ILF)
dt.prev.noCT = fm.noCT[, list(Prevalence = sum(count >= .005), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MedCount = median(count),
                             MaxCount = max(count)),
                      by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont ILF)

ls.coreMdILF = dt.prev.noCT[(Prevalence >= .75*nsamples(phyT.noCT) & MaxCount >= .001), TaxaID] # Make list of core OTUs for Hirudo ILF
phy.coreMd<-prune_taxa(ls.coreMdILF,phyT.noCT)

macDaF<-phy.coreMd




#################################################
##### Boxplot comparing AnimalSource #####
##########################################

bwdat<-testILF

#sample_data(bwdat)$AnimalSource = factor(sample_data(bwdat)$AnimalSource, levels = c(0,1,2,4,7,14,28,30,35,113,215)) # Reorder Da1F
bwdatP<-prune_taxa(c(coreTot),bwdat) # keep only core ILF taxa (box+whisker data pruned)
bwdatP<-prune_taxa(c("denovo199378","denovo113251","denovo225500","denovo52983","denovo68214"),bwdatP)
bwdatPmer<-tax_glom(bwdatP,"Genus2")

lowA<-1e-3 # set ymin
bwMelt <- psmelt(bwdatPmer) # create data.table from phyloseq object, bwdatFam (box+whisker melt)
bwGen<-as.character(get_taxa_unique(bwdatPmer, "Genus2")) # compile list of genera to examine (box + whisker genera)

bwdatP0<-subset_samples(bwdatP,AnimalSource=="GrotonMA")
bwMelt0<-psmelt(bwdatP0)

bwMelt$Q25<-0
bwMelt$Q50<-0
bwMelt$Q75<-0
bwMelt0$Q25<-0
bwMelt0$Q50<-0
bwMelt0$Q75<-0
bwMelt0$box2<-0
bwMelt0$box3<-0
bwMelt0$box4<-0
bwMelt0$conf1<-0
bwMelt0$conf2<-0

for (f in c(bwGen)){
  bwMelt0[,f]<-with(bwMelt0,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Genus2==eval(f) | AnimalSource == GrotonMA)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Genus2) %>%
  summarise(avg_Q25 = mean(Q25),
            avg_Q50 = mean(Q50),
            avg_Q75 = mean(Q75),
            hinge1 = mean(box2),
            avg_box3 = mean(box3),
            hinge2 = mean(box4),
            conf1 = mean(conf1),
            conf2 = mean(conf2))

# Enter Abundance in new taxa-specific columns in data.table (bwMelt). 
for (f in c(bwGen)){
  bwMelt[,f]<-with(bwMelt,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt,Genus2==eval(f) | AnimalSource==GrotonMA)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt). #   
  bwMelt$Q25<-with(bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q25,bwMelt$Q25))
  bwMelt$Q50<-with(bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q50,bwMelt$Q50))
  bwMelt$Q75<-with(bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q75,bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pProt <- ggplot(bwMelt) +
  ggtitle("Macrobdella ILF Taxa") + 
  stat_boxplot(data=bwMelt, aes(x=AnimalSource, y=Abundance, fill=AnimalSource)) +    
  scale_y_log10(limits=c(lowA,1)) + 
  facet_grid(~Genus2, scales="free_x",space="free") +
  geom_line(data=bwMelt,aes(x=as.numeric(AnimalSource),y=Q50),linetype=2) + 
  geom_ribbon(data=bwMelt, aes(x=as.numeric(AnimalSource), ymin=Q25, ymax=Q75),fill="red", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=rainbow)   
pProt


