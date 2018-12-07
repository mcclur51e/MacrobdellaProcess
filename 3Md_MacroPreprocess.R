#########################################################################################################
########## Pre-process data to remove low count OTUs, contaminant OTUs, and duplicate samples  ##########
#########################################################################################################
phy.In<-subset_samples(phy.sin,Age%in%c("A","M","G")) # Keep only adult, mock, and reagent/negative samples (physeq.Initial)
phy.10min<-subset_samples(phy.In,sample_sums(phy.In)>=10000) # keep samples with >=10000 reads (physeq.minimum 10,000)
dt.phyPru = fast_melt(phy.10min) # make data table from physeq object (data table. physeq Pruned)
prev.phyPru = dt.phyPru[, list(Prevalence = sum(count >= maxNeg), 
                              TotalPer = sum(count),
                              MinCount = min(count),
                              MaxCount = max(count)),
                              by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.hiPrev = prev.phyPru[(Prevalence > 2 & MaxCount >= maxNeg), TaxaID] # Make list of OTUs present in more than 2 samples and having at least maxNeg reads in at least one sample (list.high prevalence)
phy.pru2 <- prune_taxa(ls.hiPrev,phy.10min) # remove identified taxa from phyloseq object (physeq.Pruned 2)
phy.out<-subset_samples(phy.pru2, !Date_Collected=="2015_0830")
phyT.Tform<-transform_sample_counts(phy.out, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)
