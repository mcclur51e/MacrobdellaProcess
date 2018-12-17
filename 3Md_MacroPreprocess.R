#########################################################################################################
########## Pre-process data to remove low count OTUs, contaminant OTUs, and duplicate samples  ##########
#########################################################################################################
phyR.In<-subset_samples(phyR.sin,Age%in%c("A","M","G")) # Keep only adult, mock, and reagent/negative samples (physeq.Initial)
phyR.In2<-subset_samples(phyR.In,Da1F!="2") # 2Da1F samples have too low of n and are problematic in the Macrobdella dataset, so removing them entirely for this study
phyR.10min<-subset_samples(phyR.In2,sample_sums(phyR.In2)>=10000) # keep samples with >=10000 reads (physeq.minimum 10,000)
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
