# library(pairwiseAdonis) # should already be loaded

phyT.adonis<-phyT.start # assign data for permanova testing (phyloseq transformed . adonis function)
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.md<-adonis(bray.md ~ Taxonomic_ID, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

##### PERMANOVA for taxonomic groups merged at "Order" #####
phyT.adonis<-phyT.start # assign data for permanova testing (phyloseq transformed . adonis function)
phyT.adonis<-tax_glom(phyT.adonis,taxrank="Order")
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.md<-adonis(bray.md ~ Taxonomic_ID, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

###### Macrobdella decora #####
phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
phyT.adonis<-phyT.md # assign data for permanova testing (phyloseq transformed . adonis function)
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.md<-adonis(bray.md ~ Sample_Type * WildMonth * Da1Fb * AnimalSource, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

##### Affect of extraction method #####
phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
phyT.adonis<-phyT.md # assign data for permanova testing (phyloseq transformed . adonis function)
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make
perm.md<-adonis(bray.md ~ Kit_Name, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

###### Hirudo verbana #####
phyT.hv<-subset_samples(phyT.start,Taxonomic_ID=="Hverbana") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
phyT.adonis<-phyT.hv # assign data for permanova testing (phyloseq transformed . adonis function)
bray.hv<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.hv <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.hv<-adonis(bray.hv ~ Sample_Type * Da1Fb * AnimalSource, data = sd.hv, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.hv



pairwise.adonis(bray.md,sd.md$Lot_Number)

pairwise.adonis(bray.md,paste(sd.md$Lot_Number,sd.md$Sample_Type))

phyT.mdCM<-subset_samples(phyT.md,AnimalSource%in%c("Wlot","GrotonMA")) # subset to include only MA and CT samples (phyloseq transformed . macrobdella decora CT MA)
bray.mdCM <- phyloseq::distance(phyT.mdCM, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora CT MA)
sd.mdCM <- data.frame(sample_data(phyT.mdCM)) # make a data frame from the sample_data (sample data . macrobdella decora CT MA)

pairwise.adonis(bray.mdCM,paste(sd.mdCM$Sample_Type,sd.mdCM$AnimalSource))









########## Compare statistical differences between sample groups using PERMANOVA ##########
### Define groups (based off phyloseq object from 'DaF.R') ### 
#phyT.start<-subset_samples(phyT.start,!AnimalSource%in%c("CarogaNY","MtSnowVT"))
phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
sample_data(phyT.md)$isBlad<-with(sample_data(phyT.md),
                                  ifelse(Sample_Type=="Bladder","yes",
                                         "no")) # create new column to indicate if sample is from a bladder

phyT.mdMain<-subset_samples(phyT.md,AnimalSource!="CarogaNY")
phyT.mdGI<-subset_samples(phyT.md,Sample_Type%in%c("ILF","Intestinum"))
phyT.mdILF<-subset_samples(phyT.md,Sample_Type%in%c("ILF")) # keep only ILF samples
phyT.mdInt<-subset_samples(phyT.md,Sample_Type%in%c("Intestinum")) # keep only Intestinum samples
phyT.mdGI0<-subset_samples(phyT.mdGI,Da1Fb=="0")

ls.mdInthi<-taxa_names(prune_taxa(names(sort(taxa_sums(phyT.mdInt), TRUE)[1:10]),phyT.mdInt))

phyT.adonis<-phyT.mdILF
bray.mdGI<- phyloseq::distance(phyT.adonis, method = "wunifrac") # Calculate bray curtis distance matrix
sd.mdGI <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data
pairwise.adonis(bray.mdGI,paste(sd.mdGI$Da1F))

perm.mdGI<-adonis(bray.mdGI ~ Da1Fb, data = sd.mdGI, method = "bray") # Adonis test
perm.mdGI

### Determine if AnimalSource affects ILF m### 
phyT.md <- subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
phyT.mdILF <- subset_samples(phyT.md,Sample_Type=="ILF")
phyT.mdILF0 <- subset_samples(phyT.mdILF,Da1Fb%in%c("0","100"))
ls.mdILFhi<-taxa_names(prune_taxa(names(sort(taxa_sums(phyT.mdILF0), TRUE)[1:10]),phyT.mdILF))

phyT.adonis <- phyT.mdILF0
bray.mdILF <- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix
sd.mdILF <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data
pairwise.adonis(bray.mdILF,paste(sd.mdILF$AnimalSource))

dt.mdILF0 = fast_melt(phyT.mdILF0) # make data table from physeq object (data table. physeq Pruned)
prev.mdILF0 = dt.mdILF0[, list(Prevalence = sum(count >= .0001), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.mdILF0 = prev.mdILF0[(Prevalence > 1 & MaxCount >= .001), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)

phyT.adonisR <- phyT.adonis
for(otu in c(ls.mdILF0)){
  phyT.adonis<-subset_taxa(phyT.adonisR,Number==c(otu))
  print(otu)
  bray.mdILF<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix
  sd.mdILF <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data
  print(pairwise.adonis(bray.mdILF,paste(sd.mdILF$AnimalSource)))
}

### Determine if AnimalSource affects Intestinum### 
phyT.mdInt<-subset_samples(phyT.md,Sample_Type=="Intestinum")
phyT.mdInt0 <- subset_samples(phyT.mdInt,Da1Fb%in%c("0","100"))

phyT.adonis<-phyT.mdInt0
bray.mdInt<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix
sd.mdInt <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data
pairwise.adonis(bray.mdInt,paste(sd.mdInt$AnimalSource))

dt.mdInt0 = fast_melt(phyT.mdInt0) # make data table from physeq object (data table. physeq Pruned)
prev.mdInt0 = dt.mdInt0[, list(Prevalence = sum(count >= .0001), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.mdInt0 = prev.mdInt0[(Prevalence > 1 & MaxCount >= .001), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)

phyT.adonisR <- phyT.adonis
for(otu in c(ls.mdInt0)){
  phyT.adonis<-subset_taxa(phyT.adonisR,Number==c(otu))
  print(otu)
  bray.mdInt<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix
  sd.mdInt <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data
  print(pairwise.adonis(bray.mdInt,paste(sd.mdInt$AnimalSource)))
}




########## DaF in a loop ##########
### Macrobdella decora ILF by Da1F ###
ls.Da1F<-levels(factor(sample_data(phyT.mdDaFILF)$Da1Fb))
mtx.sigDaFILF<-matrix(ncol=length(ls.Da1F),nrow=length(ls.Da1F))
dimnames(mtx.sigDaFILF) = list(c(ls.Da1F),c(ls.Da1F))
for(day1 in c(ls.Da1F)){
  for(day in c(ls.Da1F)){ 
    if(day==day1){next}
    loop.mdILF<-subset_samples(phyT.mdDaFILF,Da1Fb%in%c(day1,day)) # Examine only ILF data
    bray.mdILF<- phyloseq::distance(loop.mdILF, method = "bray") # Calculate bray curtis distance matrix
    df.mdILF <- data.frame(sample_data(loop.mdILF)) # make a data frame from the sample_data
    perm.mdILF<-adonis(bray.mdILF ~ Da1Fb * AnimalSource, data = df.mdILF, method = "bray") # Adonis test
    sig.mdILF<-as.data.frame(perm.mdILF$aov.tab)["Da1F", "Pr(>F)"]
    print(perm.mdILF$aov.tab)
    mtx.sigDaFILF[day1,day] <- as.numeric(sig.mdILF)
  }
}
write.table(mtx.sigDaFILF, "tableSig_DaFMd_ILF.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included



###############################################
################ Practice area ################ 
############################################### 

# by ILF
phyT.coreILF<-subset_samples(phyT.mdDaFGIpr,Sample_Type=="ILF") # Examine only ILF data
bray.coreILF<- phyloseq::distance(phyT.coreILF, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreILF <- data.frame(sample_data(phyT.coreILF)) # make a data frame from the sample_data
perm.coreILF<-adonis(bray.coreILF ~ Taxonomic_ID, data = df.sampleCoreILF, method = "bray") # Adonis test
sig.coreILF<-as.data.frame(perm.coreILF$aov.tab)["Taxonomic_ID", "Pr(>F)"]
beta <- betadisper(bray.coreILF, df.sampleCoreILF$Taxonomic_ID)
permutest(beta)
# by Intestinum
phyT.coreInt<-subset_samples(phyT.mdDaFGIpr,Sample_Type=="Intestinum") # Examine only Intestinum data
bray.coreInt<- phyloseq::distance(phyT.coreInt, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreInt <- data.frame(sample_data(phyT.coreInt)) # make a data frame from the sample_data
perm.coreInt<-adonis(bray.coreInt ~ Taxonomic_ID, data = df.sampleCoreInt, method = "bray") # Adonis test
sig.coreInt<-as.data.frame(perm.coreInt$aov.tab)["Taxonomic_ID", "Pr(>F)"]
# by Bladder
phyT.coreBlad<-subset_samples(phyT.start,Sample_Type=="Bladder") # Examine only Bladder data
bray.coreBlad<- phyloseq::distance(phyT.coreBlad, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreBlad <- data.frame(sample_data(phyT.coreBlad)) # make a data frame from the sample_data
perm.coreBlad<-adonis(bray.coreBlad ~ Taxonomic_ID, data = df.sampleCoreBlad, method = "bray") # Adonis test
sig.coreBlad<-as.data.frame(perm.coreBlad$aov.tab)["Taxonomic_ID", "Pr(>F)"]
# by Hirudo
phyT.coreHv<-subset_samples(phyT.mdDaFGIpr,Taxonomic_ID=="Hverbana") # Examine only Hverbana data
bray.coreHv<- phyloseq::distance(phyT.coreHv, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreHv <- data.frame(sample_data(phyT.coreHv)) # make a data frame from the sample_data
perm.coreHv<-adonis(bray.coreHv ~ Sample_Type, data = df.sampleCoreHv, method = "bray") # Adonis test
sig.coreHv<-as.data.frame(perm.coreHv$aov.tab)["Sample_Type", "Pr(>F)"]
# by Macrobdella
phyT.coreMd<-subset_samples(phyT.mdDaFGIpr,Taxonomic_ID=="Mdecora") # Examine only Mdecora data
bray.coreMd<- phyloseq::distance(phyT.coreMd, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreMd <- data.frame(sample_data(phyT.coreMd)) # make a data frame from the sample_data
perm.coreMd<-adonis(bray.coreMd ~ Sample_Type, data = df.sampleCoreMd, method = "bray") # Adonis test
sig.coreMd<-as.data.frame(perm.coreMd$aov.tab)["Sample_Type", "Pr(>F)"]
# Hirudo GI
phyT.coreHvGI<-subset_samples(phyT.coreHv,Sample_Type%in%c("ILF","Intestinum")) # Examine only Hverbana ILF and Intestinum data
bray.coreHvGI<- phyloseq::distance(phyT.coreHvGI, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreHvGI <- data.frame(sample_data(phyT.coreHvGI)) # make a data frame from the sample_data
perm.coreHvGI<-adonis(bray.coreHvGI ~ Sample_Type, data = df.sampleCoreHvGI, method = "bray") # Adonis test
sig.coreHvGI<-as.data.frame(perm.coreHvGI$aov.tab)["Sample_Type", "Pr(>F)"]
sig.coreHvGI
# Macrobdella GI
phyT.coreMdGI<-subset_samples(phyT.coreMd,Sample_Type%in%c("ILF","Intestinum")) # Examine only Mdecora ILF and Intestinum data
bray.coreMdGI<- phyloseq::distance(phyT.coreMdGI, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreMdGI <- data.frame(sample_data(phyT.coreMdGI)) # make a data frame from the sample_data
perm.coreMdGI<-adonis(bray.coreMdGI ~ Sample_Type, data = df.sampleCoreMdGI, method = "bray") # Adonis test
sig.coreMdGI<-as.data.frame(perm.coreMdGI$aov.tab)["Sample_Type", "Pr(>F)"]
sig.coreMdGI



sig.coreILF
sig.coreInt
sig.coreBlad
sig.coreHv
sig.coreMd
sig.coreHvGI
sig.coreMdGI







###############################################
################ Practice area ################ 
############################################### 

#### Are LinA samples ok to use? ####


phyT.md<-subset_samples(phyT.start,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq trasnformed . macrobdella decora)
phyT.adonis<-phyT.md # assign data for permanova testing (phyloseq transformed . adonis function)
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.md<-adonis(bray.md ~ Sample_Type * Da1Fb * DNA_extractor, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

pairwise.adonis(bray.md,sd.md$DNA_extractor)

pairwise.adonis(bray.md,paste(sd.md$DNA_extractor,sd.md$Sample_Type))
