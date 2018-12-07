########## Separate environmental interaction groups and plot to evaluate if enough samples have been sequenced ##########
# 
most_abundant_taxa <- sort(taxa_sums(juv), TRUE)[1:27] # Identify 15 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),juv) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus") + facet_grid(RunDate~., scales="free_x",space="free") + scale_fill_manual(values=pairBiome)




# Control
juvCon1F<-subset_samples(juvCon,O_Da1F%in%c("0","1","2","4","7"))
juvConT<-transform_sample_counts(juvCon1F, function(x) x/sum(x)) # transform raw counts to fraction
most_abundant_taxa <- sort(taxa_sums(juvConT), TRUE)[1:15] # Identify 15 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),juvConT) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus") + scale_fill_manual(values=pairBiome)

# 2 Feeds # weird results here
juvCon2F<-subset_samples(juvCon,Da2F%in%c("0","1","2","4","7"))
juvCon2T<-transform_sample_counts(juvCon2F, function(x) x/sum(x)) # transform raw counts to fraction
most_abundant_taxa <- sort(taxa_sums(juvCon2T), TRUE)[1:15] # Identify 15 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),juvCon2T) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus") + scale_fill_manual(values=pairBiome)


### Pure Culture treatments ###
# Pure culture feed
juvTreatPure<-subset_samples(juvTreat,Subproject=="PureCulture")
juvTpT<-transform_sample_counts(juvTreatPure, function(x) x/sum(x)) # transform raw counts to fraction
most_abundant_taxa <- sort(taxa_sums(juvTpT), TRUE)[1:15] # Identify 15 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),juvTpT) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus") + scale_fill_manual(values=pairBiome)
# AeroFed 
juvAf<-subset_samples(juv1F,Treatment%in%c("AeroFed"))
most_abundant_taxa <- sort(taxa_sums(juvAf), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAf <- prune_taxa(names(most_abundant_taxa),juvAf) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAf,fill="Genus") + scale_fill_manual(values=pairBiome)
# AeroFed + 
juvAplus<-subset_samples(juv,Treatment%in%c("AeroFed_AeroFed","AeroFed_Juv","AeroFed_MucFed"))
most_abundant_taxa <- sort(taxa_sums(juvAplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAplus <- prune_taxa(names(most_abundant_taxa),juvAplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAplus,fill="Genus")
# AeroMuc + 
juvAMplus<-subset_samples(juv,Treatment%in%c("AeroMuc_AeroFed","AeroMuc_Juv","AeroMucFed"))
most_abundant_taxa <- sort(taxa_sums(juvAMplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAMplus <- prune_taxa(names(most_abundant_taxa),juvAMplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAMplus,fill="Genus")
# Juv+ #?!?!
juvPlus<-subset_samples(juv,Treatment%in%c("juv_AeroFed","juv_AeroMuc","juv_MucFed"))
most_abundant_taxa <- sort(taxa_sums(juvPlus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedPlus <- prune_taxa(names(most_abundant_taxa),juvPlus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedPlus,fill="Genus")
# Muc+ 
juvMplus<-subset_samples(juv,Treatment%in%c("MucFed","Muc_AeroFed","MucFed_Juv"))
most_abundant_taxa <- sort(taxa_sums(juvMplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMplus <- prune_taxa(names(most_abundant_taxa),juvMplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMplus,fill="Genus")
# AeroWater 
juvAw<-subset_samples(juv1F,Treatment%in%c("AeroWater"))
most_abundant_taxa <- sort(taxa_sums(juvAw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAw <- prune_taxa(names(most_abundant_taxa),juvAw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAw,fill="Genus") + scale_fill_manual(values=pairBiome)
# AeroWaterDilute # this should have been in 1223 data
juvAwd<-subset_samples(juv,Treatment%in%c("AeroWater2","AeroWater3","AeroWater4","AeroWater5","AeroWater6"))
most_abundant_taxa <- sort(taxa_sums(juvAwd), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAwd <- prune_taxa(names(most_abundant_taxa),juvAwd) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAwd,fill="Genus") + scale_fill_manual(values=pairBiome)
# MucWater # no samples over 10,000 reads
juvMw<-subset_samples(physeq,Treatment=="MucWater")
most_abundant_taxa <- sort(taxa_sums(juvMw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMw <- prune_taxa(names(most_abundant_taxa),juvMw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMw,fill="Genus") + scale_fill_manual(values=pairBiome)
### Pure Culture treatments ###




# Dialysis 
juvDi<-subset_samples(juv,Treatment=="Dialysis")
most_abundant_taxa <- sort(taxa_sums(juvDi), TRUE)[1:10] # Identify 10 most abundant taxa
prunedDi <- prune_taxa(names(most_abundant_taxa),juvDi) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedDi,fill="Genus") + scale_fill_manual(values=pairBiome)

# MacroAdult
juvMa<-subset_samples(juv,Treatment=="MacroAdult")
juvMaI<-subset_samples(juvMa,Sample_Type=="ILF")
juvMaI<-subset_samples(juvMaI,!sample_names(juvMaI)%in%c("MAj4dF110514EMb.iD133","MAj4dF110514EMd.iD105")) # Remove duplicates
most_abundant_taxa <- sort(taxa_sums(juvMaI), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMa <- prune_taxa(names(most_abundant_taxa),juvMaI) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMa,fill="Genus") + 
  scale_fill_manual(values=pairBiome)

# USAwater # carefully check MiSeq run dates
juvUw<-subset_samples(juv,Treatment=="USAwater")
juvUw1F<-subset_samples(juvUw,Da2F=="none")
most_abundant_taxa <- sort(taxa_sums(juvUw1F), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUw <- prune_taxa(names(most_abundant_taxa),juvUw1F) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUw,fill="Genus") + scale_fill_manual(values=pairBiome)
# fWater
juvFw<-subset_samples(juv,Treatment=="fWater")
juvFw1F<-subset_samples(juvFw,Da1F%in%c("0","1","2","4","7","110"))
most_abundant_taxa <- sort(taxa_sums(juvFw1F), TRUE)[1:10] # Identify 10 most abundant taxa
prunedfw <- prune_taxa(names(most_abundant_taxa),juvFw1F) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedfw,fill="Genus") + scale_fill_manual(values=pairBiome)
# Pure culture feed
juvTreatPure<-subset_samples(juvTreat,Subproject=="PureCulture")
juvTpT<-transform_sample_counts(juvTreatPure, function(x) x/sum(x)) # transform raw counts to fraction
most_abundant_taxa <- sort(taxa_sums(juvTpT), TRUE)[1:10] # Identify 10 most abundant taxa
pruned <- prune_taxa(names(most_abundant_taxa),juvTpT) # Create a subset of data including only 10 most abundant taxa
plot_bar(pruned,fill="Genus")

# Dialysis # no samples over 10,000 reads
juvDi<-subset_samples(juv,Treatment=="Dialysis")
most_abundant_taxa <- sort(taxa_sums(juvDi), TRUE)[1:10] # Identify 10 most abundant taxa
prunedDi <- prune_taxa(names(most_abundant_taxa),juvDi) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedDi,fill="Genus")

# fWater
juvfw<-subset_samples(juv,Treatment=="fWater")
most_abundant_taxa <- sort(taxa_sums(juvfw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedfw <- prune_taxa(names(most_abundant_taxa),juvfw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedfw,fill="Genus")

# MacroAdult
juvMa<-subset_samples(juv,Treatment=="MacroAdult")
most_abundant_taxa <- sort(taxa_sums(juvMa), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMa <- prune_taxa(names(most_abundant_taxa),juvMa) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMa,fill="Genus")

# MucWater # no samples over 10,000 reads
juvMw<-subset_samples(juv,Treatment=="MucWater")
most_abundant_taxa <- sort(taxa_sums(juvMw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMw <- prune_taxa(names(most_abundant_taxa),juvMw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMw,fill="Genus")

# StarveAdult # should this experiment be trashed?
juvSa<-subset_samples(juv,Treatment=="StarveAdult")
most_abundant_taxa <- sort(taxa_sums(juvSa), TRUE)[1:10] # Identify 10 most abundant taxa
prunedSa <- prune_taxa(names(most_abundant_taxa),juvSa) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedSa,fill="Genus")

# USAhatch # should this experiment be trashed?
juvUh<-subset_samples(juv,Treatment=="USAhatch")
most_abundant_taxa <- sort(taxa_sums(juvUh), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUh <- prune_taxa(names(most_abundant_taxa),juvUh) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUh,fill="Genus")

# USArock # carefully check MiSeq run dates
juvUr<-subset_samples(juv,Treatment=="USArock")
most_abundant_taxa <- sort(taxa_sums(juvUr), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUr <- prune_taxa(names(most_abundant_taxa),juvUr) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUr,fill="Genus")

# USAwater # carefully check MiSeq run dates
juvUw<-subset_samples(juv,Treatment=="USAwater")
most_abundant_taxa <- sort(taxa_sums(juvUw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUw <- prune_taxa(names(most_abundant_taxa),juvUw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUw,fill="Genus")

# MealType 
juvMt<-subset_samples(juv,Treatment%in%c("Blood","InactiveBlood","Plasma"))
most_abundant_taxa <- sort(taxa_sums(juvMt), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMt <- prune_taxa(names(most_abundant_taxa),juvMt) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMt,fill="Genus") + scale_fill_manual(values=pairBiome)
plot_bar(prunedMt,fill="Genus")

# AeroHatch # meal lot contaminated with Salmonella?
juvAh<-subset_samples(juv,Treatment%in%c("AeroHatch"))
most_abundant_taxa <- sort(taxa_sums(juvAh), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAh <- prune_taxa(names(most_abundant_taxa),juvAh) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAh,fill="Number") + scale_fill_manual(values=pairBiome)

# USArock # carefully check MiSeq run dates
juvUr<-subset_samples(juv,Treatment=="USArock")
most_abundant_taxa <- sort(taxa_sums(juvUr), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUr <- prune_taxa(names(most_abundant_taxa),juvUr) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUr,fill="Genus") + scale_fill_manual(values=pairBiome)

# WlotRock # GREAT! Ready to go
juvWr<-subset_samples(juv,Treatment%in%c("WlotRock"))
juvWr1F<-subset_samples(juvWr,Da1F%in%c("0","1","2","4","7","110"))
most_abundant_taxa <- sort(taxa_sums(juvWr1F), TRUE)[1:10] # Identify 10 most abundant taxa
prunedWr1F <- prune_taxa(names(most_abundant_taxa),juvWr1F) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedWr1F,fill="Genus") + scale_fill_manual(values=pairBiome)




# Adult bladder
bladder<-subset_samples(physeqAn,Sample_Type=="Bladder") #
bladderT<-transform_sample_counts(bladder, function(x) x/sum(x)) # transform raw cou
most_abundant_taxa <- sort(taxa_sums(bladderT), TRUE)[1:15] # Identify 10 most abundant taxa
prunedBl <- prune_taxa(names(most_abundant_taxa),bladderT) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedBl,fill="Genus") + scale_fill_manual(values=pairBiome)


#!#!#!#!#!#!#!# PROBLEMS HERE #!#!#!#!#!#!#!#
# StarveAdult # should this experiment be trashed?
juvSa<-subset_samples(juv,Treatment=="StarveAdult")
most_abundant_taxa <- sort(taxa_sums(juvSa), TRUE)[1:10] # Identify 10 most abundant taxa
prunedSa <- prune_taxa(names(most_abundant_taxa),juvSa) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedSa,fill="Genus") + scale_fill_manual(values=pairBiome)

# USAhatch # should this experiment be trashed?
juvUh<-subset_samples(juv,Treatment=="USAhatch")
most_abundant_taxa <- sort(taxa_sums(juvUh), TRUE)[1:10] # Identify 10 most abundant taxa
prunedUh <- prune_taxa(names(most_abundant_taxa),juvUh) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedUh,fill="Genus") + scale_fill_manual(values=pairBiome)
#!#!#!#!#!#!#!# PROBLEMS HERE #!#!#!#!#!#!#!#
plot_bar(prunedAh,fill="Genus")

# WlotRock # GREAT! Ready to go
juvWr<-subset_samples(juv,Treatment%in%c("WlotRock"))
most_abundant_taxa <- sort(taxa_sums(juvWr), TRUE)[1:10] # Identify 10 most abundant taxa
prunedWr <- prune_taxa(names(most_abundant_taxa),juvWr) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedWr,fill="Genus")

# AeroFed 
juvAf<-subset_samples(juv,Treatment%in%c("AeroFed"))
most_abundant_taxa <- sort(taxa_sums(juvAf), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAf <- prune_taxa(names(most_abundant_taxa),juvAf) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAf,fill="Genus")

# AeroFed + 
juvAplus<-subset_samples(juv,Treatment%in%c("AeroFed_AeroFed","AeroFed_Juv","AeroFed_MucFed"))
most_abundant_taxa <- sort(taxa_sums(juvAplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAplus <- prune_taxa(names(most_abundant_taxa),juvAplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAplus,fill="Genus")

# AeroMuc + 
juvAMplus<-subset_samples(juv,Treatment%in%c("AeroMuc_AeroFed","AeroMuc_Juv","AeroMucFed"))
most_abundant_taxa <- sort(taxa_sums(juvAMplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAMplus <- prune_taxa(names(most_abundant_taxa),juvAMplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAMplus,fill="Genus")

# AeroWater 
juvAw<-subset_samples(juv,Treatment%in%c("AeroWater"))
most_abundant_taxa <- sort(taxa_sums(juvAw), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAw <- prune_taxa(names(most_abundant_taxa),juvAw) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAw,fill="Genus")

# AeroWaterDilute # this should have been in 1223 data
juvAwd<-subset_samples(juv,Treatment%in%c("AeroWater2","AeroWater3","AeroWater4","AeroWater5","AeroWater6"))
most_abundant_taxa <- sort(taxa_sums(juvAwd), TRUE)[1:10] # Identify 10 most abundant taxa
prunedAwd <- prune_taxa(names(most_abundant_taxa),juvAwd) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedAwd,fill="Genus")

# Juv+ #?!?!
juvPlus<-subset_samples(juv,Treatment%in%c("juv_AeroFed","juv_AeroMuc","juv_MucFed"))
most_abundant_taxa <- sort(taxa_sums(juvPlus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedPlus <- prune_taxa(names(most_abundant_taxa),juvPlus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedPlus,fill="Genus")

# Muc+ 
juvMplus<-subset_samples(juv,Treatment%in%c("MucFed","Muc_AeroFed","MucFed_Juv"))
most_abundant_taxa <- sort(taxa_sums(juvMplus), TRUE)[1:10] # Identify 10 most abundant taxa
prunedMplus <- prune_taxa(names(most_abundant_taxa),juvMplus) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedMplus,fill="Genus")

# Adult bladder
bladder<-subset_samples(physeqAn,SampleSite=="Bladder") #
bladderT<-transform_sample_counts(bladder, function(x) x/sum(x)) # transform raw cou
most_abundant_taxa <- sort(taxa_sums(bladderT), TRUE)[1:10] # Identify 10 most abundant taxa
prunedBl <- prune_taxa(names(most_abundant_taxa),bladderT) # Create a subset of data including only 10 most abundant taxa
plot_bar(prunedBl,fill="Genus")
