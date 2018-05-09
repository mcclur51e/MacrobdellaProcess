########## Pre-process Data ##########
toolow<-prune_taxa(taxa_sums(physeq)>100,physeq) # Remove OTUs with < 100 reads over all samples (1% of 1 sample)
#outlier3<-subset_samples(physeq,!DateHatch%in%c("2014_0914","2014_0911"))
outlier3<-subset_samples(toolow,!sample_names(toolow)%in%c("A1dF100913EM.iD39","A4dF101213EM.iD73","j7dH092114EMb.iD76","j7dH091814EMc.iD53","j2d2F020715EMc.iD70","Ba0dF110614EMc.iD11","Ba0dF110614EMa.iD148"))
outlier2<-subset_samples(outlier3,!sample_names(outlier3)%in%c("FAa1dF100913EMa.iD6","FAa1dF100913EMb.iD9","SAa1dF100913EMb.iD141","A042317EMj.iD441","A042317EMk.iD444","A042317EMe.uD542","A042317EMh.uD558","A042317EMg.uD569","A042317EMi.uD560","A042317EMg.iD690","A042317EMf.iD688","A042317EMe.iD689","A042317EMj.bD440","A042317EMk.bD445","A042117JGa.bD438","URj4dF110514EMa.iD39","URj4dF110514EMb.iD40","URj4dF110514EMc.iD225","UWj4dF110514EMe.iD72")) # Remove duplicates and bad GrafJ samples
outlier<-subset_samples(outlier2, !Date_Collected=="2015_0830")
physeq10<-subset_samples(outlier,sample_sums(outlier)>=10000) # keep samples with >=10000 reads
#physeqG<-tax_glom(physeq10,taxrank="Genus")
physeqT<-transform_sample_counts(physeq10, function(x) x/sum(x)) # transform raw counts to fraction

### Remove Halomonas-contaminated samples ###
pHalo<-subset_taxa(physeqT, Family=="Halomonadaceae") # make phyloseq object with only Halomonadaceae
listH<-prune_samples(sample_sums(pHalo)<.012,pHalo) # list samples with <1% Halomonas
physeqH<-prune_samples(sample_names(listH),physeqT) # keep samples with <= 1% Halomonas
### Remove Anaplasma-contaminated samples ###
pAna<-subset_taxa(physeqT,Family=="Anaplasmataceae") # make phyloseq object with only Anaplasmataceae
listA<-prune_samples(sample_sums(pAna)<.012,pAna) # list samples with <1% Anaplasma
physeqAn<-prune_samples(sample_names(listA),physeqH) # keep samples with <= 1% Anaplasmataceae
#physeqMer<-tax_glom(physeqAn,"Genus")
#save(physeqMer,file=("DataFiles/physeq_currentMer.RData")) # Save the phyloseq data object in a .RData file 
########## Pre-process Data ##########

### Common groups ###
juvRaw<-subset_samples(physeqT,Age=="J") # to be used to check for samples that may have been eliminated due to Anaplasma/Halomonas contamination
juv<-subset_samples(physeqAn,Age=="J")
juvTreat<-subset_samples(juv,!Treatment%in%c("BBEZadults","USAadult","1feed","2feeds","Refused"))
juvCon<-subset_samples(juv,Treatment%in%c("1feed","2feeds"))
juv1F<-subset_samples(juv,Da2F=="none")

hatch<-subset_samples(physeqAn,Age=="H")
