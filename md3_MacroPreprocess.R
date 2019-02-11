#########################################################################################################
########## Pre-process data to remove low count OTUs, contaminant OTUs, and duplicate samples  ##########
#########################################################################################################
phyR.In<-subset_samples(phyR.sin,Age%in%c("A","M","G")) # Keep only adult, mock, and reagent/negative samples (physeq.Initial)
phyR.In1<-subset_samples(phyR.In,!Sample_Type%in%c("Ovary")) # remove ovary samples (physeq bladder, ILF, intestinum) 
#phyR.In2<-subset_samples(phyR.In1,Da1F!="2") # 2Da1F samples have too low of n and are problematic in the Macrobdella dataset, so removing them entirely for this study
phyR.10min<-subset_samples(phyR.In1,sample_sums(phyR.In1)>=5000) # keep samples with >=10000 reads (physeq.minimum 10,000)
dt.phyPru = fast_melt(phyR.10min) # make data table from physeq object (data table. physeq Pruned)
prev.phyPru = dt.phyPru[, list(Prevalence = sum(count >= var.maxNeg), 
                              TotalPer = sum(count),
                              MinCount = min(count),
                              MaxCount = max(count)),
                              by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.Pres = prev.phyPru[(Prevalence > 0 & MaxCount >= var.maxNeg), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)

phyR.pru2 <- prune_taxa(ls.Pres,phyR.10min) # remove identified taxa from phyloseq object (physeq.Pruned 2)
phyR.out<-subset_samples(phyR.pru2, !Date_Collected=="2015_0830")
phyT.Tform<-transform_sample_counts(phyR.out, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)

### Add new column to map data. 'Header' column will be used for many figures
sample_data(phyT.Tform)$Header<-with(sample_data(phyT.Tform),
                                     ifelse(AnimalSource=="USA","",
                                     ifelse(AnimalSource=="BBEZ","",
                                     ifelse(AnimalSource=="GrotonMA","MA",
                                     ifelse(AnimalSource=="CarogaNY","NY",
                                     ifelse(AnimalSource=="Wlot","CT",
                                     ifelse(AnimalSource=="MtSnowVT","VT",
                                     ifelse(AnimalSource=="Schoolhouse_Brook","CT",     
                                     as.character(AnimalSource)))))))))
### Add new column to map data. 'Da1fb' column will be used to merge some time points together
sample_data(phyT.Tform)$Da1Fb<-with(sample_data(phyT.Tform),
                                    ifelse(Da1F=="8","7",
                                    ifelse(Da1F=="Unk","4",
                                    ifelse(Da1F=="31","30",
                                    ifelse(Da1F=="35","30",
                                    ifelse(Da1F=="90","90+",
                                    ifelse(Da1F=="99","90+",
                                    ifelse(Da1F=="100","90+",
                                    ifelse(Da1F=="101","90+",
                                    ifelse(Da1F=="108","90+",
                                    ifelse(Da1F=="110","90+",
                                    ifelse(Da1F=="113","90+",
                                    ifelse(Da1F=="130","90+",
                                    ifelse(Da1F=="165","90+",
                                    ifelse(Da1F=="215","90+",
                                    as.character(Da1F)))))))))))))))) 
### Add new column to map data. 'Da1fb' column will be used to merge some time points together
sample_data(phyR.out)$Da1Fb<-with(sample_data(phyR.out),
                                  ifelse(Da1F=="8","7",
                                  ifelse(Da1F=="31","30",
                                  ifelse(Da1F=="35","30",
                                  ifelse(Da1F=="90","90+",
                                  ifelse(Da1F=="99","90+",
                                  ifelse(Da1F=="100","90+",
                                  ifelse(Da1F=="101","90+",
                                  ifelse(Da1F=="108","90+",
                                  ifelse(Da1F=="110","90+",
                                  ifelse(Da1F=="113","90+",
                                  ifelse(Da1F=="130","90+",
                                  ifelse(Da1F=="165","90+",
                                  ifelse(Da1F=="215","90+",
                                  as.character(Da1F)))))))))))))))


