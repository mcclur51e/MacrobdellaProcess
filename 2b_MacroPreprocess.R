########## Pre-process Data ##########
### Remove taxa in fewer than 3 (total) samples and with fewer than 100 reads (<=1%) in at least one sample
phy.In<-subset_samples(physeq,Age%in%c("A","M","G")) # Keep only adult, mock, and reagent/negative samples (physeq.Initial)
phy.Pru<-prune_taxa(taxa_sums(phy.In)>100,phy.In) # Remove OTUs with < 100 reads over all samples (1% of 1 sample) (physeq.Pruned)
dt.phyPru = fast_melt(phy.Pru) # make data table from physeq object (data table. physeq Pruned)
prev.phyPru = dt.phyPru[, list(Prevalence = sum(count >= cMin), 
                              TotalPer = sum(count),
                              MinCount = min(count),
                              MaxCount = max(count)),
                              by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
list.hiPrev = prev.phyPru[(Prevalence > 2 & MaxCount >= 10), TaxaID] # Make list of OTUs present in more than 2 samples and having at least 100 read (>= 1%) in at least one sample (list.high prevalence)
phy.Pru2 <- prune_taxa(list.hiPrev,phy.Pru) # remove identified taxa from phyloseq object (physeq.Pruned 2)

#outlier3<-subset_samples(physeq,!DateHatch%in%c("2014_0914","2014_0911"))
phy.out3<-subset_samples(phy.Pru2,!sample_names(phy.Pru2)%in%c("A1dF100913EM.iD39","A4dF101213EM.iD73","Ba0dF110614EMc.iD11","Ba0dF110614EMa.iD148"))
phy.out2<-subset_samples(phy.out3,!sample_names(phy.out3)%in%c("FAa1dF100913EMa.iD6","FAa1dF100913EMb.iD9","SAa1dF100913EMb.iD141","A042317EMj.iD441","A042317EMk.iD444","A042317EMe.uD542","A042317EMh.uD558","A042317EMg.uD569","A042317EMi.uD560","A042317EMg.iD690","A042317EMf.iD688","A042317EMe.iD689","A042317EMj.bD440","A042317EMk.bD445","A042117JGa.bD438")) # Remove duplicates and bad GrafJ samples
phy.out1<-subset_samples(phy.out2, !Date_Collected=="2015_0830")
phy.10min<-subset_samples(phy.out1,sample_sums(phy.out1)>=10000) # keep samples with >=10000 reads (physeq.minimum 10,000)
#physeqG<-tax_glom(physeq10,taxrank="Genus")
phy.Tform<-transform_sample_counts(phy.10min, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)

### Remove Halomonas-contaminated samples ###
phy.Hal<-subset_taxa(phy.Tform, Family=="Halomonadaceae") # make phyloseq object with only Halomonadaceae (physeq.Halomonadaceae)
list.Hal<-prune_samples(sample_sums(phy.Hal)<.012,phy.Hal) # list samples with <1% Halomonas (list.Halomonadaceae)
phy.noHal<-prune_samples(sample_names(list.Hal),phy.Tform) # keep samples with <= 1% Halomonas (physeq.no Halomonadaceae)
### Remove Anaplasma-contaminated samples ###
phy.Ana<-subset_taxa(phy.Tform,Family=="Anaplasmataceae") # make phyloseq object with only Anaplasmataceae (physeq.Anaplasmataceae)
list.Ana<-prune_samples(sample_sums(phy.Ana)<.01,phy.Ana) # list samples with <1% Anaplasma (list.Anaplasmataceae)
phy.noAn<-prune_samples(sample_names(list.Ana),phy.noHal) # keep samples with <= 1% Anaplasmataceae (physeq.no Anaplasmataceae)
#physeqMer<-tax_glom(physeqAn,"Genus")
#save(physeqMer,file=("DataFiles/physeq_currentMer.RData")) # Save the phyloseq data object in a .RData file 
########## Pre-process Data ##########

### Common groups ###
