########## Compare statistical differences between sample groups using PERMANOVA ##########
### Define groups (based off phyloseq object from 'DaF.R') ### 
macDafGI<-subset_samples(macDaf,O_Sample_Type%in%c("ILF","Intestinum")) # Keep only ILF and intestinum samples
macDaFp<-prune_taxa(taxa_sums(macDafGI)>.0,macDafGI) # keep taxa present in at least 1 sample
MdDaFILF<-subset_samples(macDaFp,O_Sample_Type%in%c("ILF")) # keep only ILF samples
MdDaFInt<-subset_samples(macDaFp,O_Sample_Type%in%c("Intestinum")) # keep only Intestinum samples


# by ILF
coreILF<-subset_samples(coreMdHv,O_Sample_Type=="ILF") # Examine only ILF data
bray.coreILF<- phyloseq::distance(coreILF, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreILF <- data.frame(sample_data(coreILF)) # make a data frame from the sample_data
perm.coreILF<-adonis(bray.coreILF ~ O_Taxonomic_ID, data = df.sampleCoreILF, method = "bray") # Adonis test
sig.coreILF<-as.data.frame(perm.coreILF$aov.tab)["O_Taxonomic_ID", "Pr(>F)"]
beta <- betadisper(bray.coreILF, df.sampleCoreILF$O_Taxonomic_ID)
permutest(beta)
# by Intestinum
coreInt<-subset_samples(coreMdHv,O_Sample_Type=="Intestinum") # Examine only Intestinum data
bray.coreInt<- phyloseq::distance(coreInt, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreInt <- data.frame(sample_data(coreInt)) # make a data frame from the sample_data
perm.coreInt<-adonis(bray.coreInt ~ O_Taxonomic_ID, data = df.sampleCoreInt, method = "bray") # Adonis test
sig.coreInt<-as.data.frame(perm.coreInt$aov.tab)["O_Taxonomic_ID", "Pr(>F)"]
# by Bladder
coreBlad<-subset_samples(coreMdHv,O_Sample_Type=="Bladder") # Examine only Bladder data
bray.coreBlad<- phyloseq::distance(coreBlad, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreBlad <- data.frame(sample_data(coreBlad)) # make a data frame from the sample_data
perm.coreBlad<-adonis(bray.coreBlad ~ O_Taxonomic_ID, data = df.sampleCoreBlad, method = "bray") # Adonis test
sig.coreBlad<-as.data.frame(perm.coreBlad$aov.tab)["O_Taxonomic_ID", "Pr(>F)"]
# by Hirudo
coreHv<-subset_samples(coreMdHv,O_Taxonomic_ID=="Hverbana") # Examine only Hverbana data
bray.coreHv<- phyloseq::distance(coreHv, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreHv <- data.frame(sample_data(coreHv)) # make a data frame from the sample_data
perm.coreHv<-adonis(bray.coreHv ~ O_Sample_Type, data = df.sampleCoreHv, method = "bray") # Adonis test
sig.coreHv<-as.data.frame(perm.coreHv$aov.tab)["O_Sample_Type", "Pr(>F)"]
# by Macrobdella
coreMd<-subset_samples(coreMdHv,O_Taxonomic_ID=="Mdecora") # Examine only Mdecora data
bray.coreMd<- phyloseq::distance(coreMd, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreMd <- data.frame(sample_data(coreMd)) # make a data frame from the sample_data
perm.coreMd<-adonis(bray.coreMd ~ O_Sample_Type, data = df.sampleCoreMd, method = "bray") # Adonis test
sig.coreMd<-as.data.frame(perm.coreMd$aov.tab)["O_Sample_Type", "Pr(>F)"]
# Hirudo GI
coreHvGI<-subset_samples(coreHv,O_Sample_Type%in%c("ILF","Intestinum")) # Examine only Hverbana ILF and Intestinum data
bray.coreHvGI<- phyloseq::distance(coreHvGI, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreHvGI <- data.frame(sample_data(coreHvGI)) # make a data frame from the sample_data
perm.coreHvGI<-adonis(bray.coreHvGI ~ O_Sample_Type, data = df.sampleCoreHvGI, method = "bray") # Adonis test
sig.coreHvGI<-as.data.frame(perm.coreHvGI$aov.tab)["O_Sample_Type", "Pr(>F)"]
sig.coreHvGI
# Macrobdella GI
coreMdGI<-subset_samples(coreMd,O_Sample_Type%in%c("ILF","Intestinum")) # Examine only Mdecora ILF and Intestinum data
bray.coreMdGI<- phyloseq::distance(coreMdGI, method = "bray") # Calculate bray curtis distance matrix
df.sampleCoreMdGI <- data.frame(sample_data(coreMdGI)) # make a data frame from the sample_data
perm.coreMdGI<-adonis(bray.coreMdGI ~ O_Sample_Type, data = df.sampleCoreMdGI, method = "bray") # Adonis test
sig.coreMdGI<-as.data.frame(perm.coreMdGI$aov.tab)["O_Sample_Type", "Pr(>F)"]
sig.coreMdGI



sig.coreILF
sig.coreInt
sig.coreBlad
sig.coreHv
sig.coreMd
sig.coreHvGI
sig.coreMdGI

########## DaF in a loop ##########
### Macrobdella decora ILF by Da1F ###
list.Da1F<-levels(factor(sample_data(MdDaFILF)$O_Da1F))
mat.sigDaFILF<-matrix(ncol=length(list.Da1F),nrow=length(list.Da1F))
dimnames(mat.sigDaFILF) = list(c(list.Da1F),c(list.Da1F))
for(day1 in c(list.Da1F)){
  print(day1)
  for(day in c(list.Da1F)){ 
    print(day)
    if(day==day1){next}
    MdLoopILF<-subset_samples(MdDaFILF,O_Da1F%in%c(day1,day)) # Examine only ILF data
    bray.MdLoopILF<- phyloseq::distance(MdLoopILF, method = "bray") # Calculate bray curtis distance matrix
    df.sampleMdLoopILF <- data.frame(sample_data(MdLoopILF)) # make a data frame from the sample_data
    perm.MdLoopILF<-adonis(bray.MdLoopILF ~ O_Da1F, data = df.sampleMdLoopILF, method = "bray") # Adonis test
    sig.MdLoopILF<-as.data.frame(perm.MdLoopILF$aov.tab)["O_Da1F", "Pr(>F)"]
    print(sig.MdLoopILF)
    mat.sigDaFILF[day1,day] <- as.numeric(sig.MdLoopILF)
  }
}
write.table(mat.sigDaFILF, "tableSig_DaFMd_ILF.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included

### Macrobdella decora Intestinum by Da1F ###
list.Da1F<-levels(factor(sample_data(MdDaFInt)$O_Da1F))
mat.sigDaFInt<-matrix(ncol=length(list.Da1F),nrow=length(list.Da1F))
dimnames(mat.sigDaFInt) = list(c(list.Da1F),c(list.Da1F))
for(day1 in c(list.Da1F)){
  print(day1)
  for(day in c(list.Da1F)){ 
    print(day)
    if(day==day1){next}
    MdLoopInt<-subset_samples(MdDaFInt,O_Da1F%in%c(day1,day)) # Examine only Int data
    bray.MdLoopInt<- phyloseq::distance(MdLoopInt, method = "bray") # Calculate bray curtis distance matrix
    df.sampleMdLoopInt <- data.frame(sample_data(MdLoopInt)) # make a data frame from the sample_data
    perm.MdLoopInt<-adonis(bray.MdLoopInt ~ O_Da1F, data = df.sampleMdLoopInt, method = "bray") # Adonis test
    sig.MdLoopInt<-as.data.frame(perm.MdLoopInt$aov.tab)["O_Da1F", "Pr(>F)"]
    print(sig.MdLoopInt)
    mat.sigDaFInt[day1,day] <- as.numeric(sig.MdLoopInt)
  }
}
write.table(mat.sigDaFInt, "tableSig_DaFMd_Intestinum.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included

### Macrobdella decora by sample type ###
phy.Md<-subset_samples(phyBasC,O_Taxonomic_ID=="Mdecora") # Keep only Macrobdella samples
list.Type<-levels(factor(sample_data(phy.Md)$O_Sample_Type))
mat.sigTypeMd<-matrix(ncol=length(list.Type),nrow=length(list.Type))
dimnames(mat.sigTypeMd) = list(c(list.Type),c(list.Type))
for(t1 in c(list.Type)){
  print(t1)
  for(type in c(list.Type)){ 
    print(type)
    if(type==t1){next}
    MdLoopType<-subset_samples(phy.Md,O_Sample_Type%in%c(t1,type)) # Examine only Type data
    bray.MdLoopType<- phyloseq::distance(MdLoopType, method = "bray") # Calculate bray curtis distance matrix
    df.sampleMdLoopType <- data.frame(sample_data(MdLoopType)) # make a data frame from the sample_data
    perm.MdLoopType<-adonis(bray.MdLoopType ~ O_Sample_Type, data = df.sampleMdLoopType, method = "bray") # Adonis test
    sig.MdLoopType<-as.data.frame(perm.MdLoopType$aov.tab)["O_Sample_Type", "Pr(>F)"]
    print(sig.MdLoopType)
    mat.sigTypeMd[t1,type] <- as.numeric(sig.MdLoopType)
  }
}
write.table(mat.sigTypeMd, "tableSig_typeMd.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included

### Hirudo verbana by sample type ###
phy.Hv<-subset_samples(phyBasC,O_Taxonomic_ID=="Hverbana") # Keep only Macrobdella samples
list.Type<-levels(factor(sample_data(phy.Hv)$O_Sample_Type))
mat.sigTypeHv<-matrix(ncol=length(list.Type),nrow=length(list.Type))
dimnames(mat.sigTypeHv) = list(c(list.Type),c(list.Type))
for(t1 in c(list.Type)){
  print(t1)
  for(type in c(list.Type)){ 
    print(type)
    if(type==t1){next}
    HvLoopType<-subset_samples(phy.Hv,O_Sample_Type%in%c(t1,type)) # Examine only Type data
    bray.HvLoopType<- phyloseq::distance(HvLoopType, method = "bray") # Calculate bray curtis distance matrix
    df.sampleHvLoopType <- data.frame(sample_data(HvLoopType)) # make a data frame from the sample_data
    perm.HvLoopType<-adonis(bray.HvLoopType ~ O_Sample_Type, data = df.sampleHvLoopType, method = "bray") # Adonis test
    sig.HvLoopType<-as.data.frame(perm.HvLoopType$aov.tab)["O_Sample_Type", "Pr(>F)"]
    print(sig.HvLoopType)
    mat.sigTypeHv[t1,type] <- as.numeric(sig.HvLoopType)
  }
}
write.table(mat.sigTypeHv, "tableSig_typeHv.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included


