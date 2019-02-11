# Thank you to jeffkimbrel for the next two lines! #
ls.taxaDC <- setdiff(taxa_names(phyT.Tform), taxa_names(phyR.allc)) # find the difference between taxa in contaminant list and physeqAn (taxa decontam)
phyT.DC <- prune_taxa(ls.taxaDC, phyT.Tform) # keep only taxa not in taxaNC (physeq decontam)

phyT.Adult<-subset_samples(phyT.DC,Age=="A") # keep adult samples (prune Adult)
phyT.Hv<-subset_samples(phyT.DC,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A042317EMe.u","A042317EMe.i","A042317EMf.u","A042317EMh.u","A042117JGbb","A042317EMi","A042317EMj","A042317EMj.b","A042317EMj.i","A042317EMj.u","A042317EMk","A042317EMk.b","A042317EMk.i","A043017EMl","A043017EMm.u","A42hF050317EMb","A050217EMo","A050217EMo.u","A050217EMr.i","A042117JGb.b","A042117JGb","A042317EMi.b","A042317EMi","A043017EMl","A043017EMl.b","A050217EMq"))
phyT.Md<-subset_samples(phyT.Adult,Taxonomic_ID%in%c("Mdecora","Unk")) # subset containing only Macrobdella samples
phyT.mdCT<-subset_samples(phyT.Md,AnimalSource=="Wlot") # subset containing only W-lot samples (CT Macrobdella)
phyT.mdMA<-subset_samples(phyT.Md,AnimalSource=="GrotonMA") # subset containing only MA samples (MA Macrobdella)
phyT.mdNY<-subset_samples(phyT.Md,AnimalSource=="CarogaNY") # subset containing only NY samples (NY Macrobdella)
phyT.mdVT<-subset_samples(phyT.Md,AnimalSource=="MtSnowVT") # subset containing only VT samples (VT Macrobdella)
#phyT.mdSB<-subset_samples(phyT.Md,AnimalSource=="Schoolhouse_Brook") # subset containing only Schoolhouse Brook samples (VT Macrobdella)
phyT.start<-merge_phyloseq(phyT.Md,phyT.Hv) # merge Macrobdella and Hirudo phyloseq objects (physeq control)

cMin<-0.001 # define minimum to count prevalence in a sample

### Hirudo ###
phyT.hvBlad<-subset_samples(phyT.Hv,Sample_Type=="Bladder") # (Hirudo verbana bladder)
fm.hvBlad = fast_melt(phyT.hvBlad) # (fast melt Hirudo verbana bladder)
prevdt.HvBlad = fm.hvBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
phyT.hvInt<-subset_samples(phyT.Hv,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.hvInt = fast_melt(phyT.hvInt) # (fast melt Hirudo verbana intestinum)
prevdt.HvInt = fm.hvInt[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
phyT.hvILF<-subset_samples(phyT.Hv,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.hvILF = fast_melt(phyT.hvILF) # (fast melt Hirudo verbana ILF)
prevdt.HvILF = fm.hvILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###Macrobdella		
phyT.mdBlad<-subset_samples(phyT.Md,Sample_Type=="Bladder") # (Hirudo verbana bladder)
phyT.mdBlad2<-subset_samples(phyT.mdBlad,!Da1F%in%c("1","4","7"))
phyT.mdBlad3<-subset_samples(phyT.mdBlad2,!sample_names(phyT.mdBlad2)%in%c("MdMA0dF040618EMfBD","MdMA0dF040618EMdBD722","MdMA0dF040618EMeBD"))
fm.mdBlad = fast_melt(phyT.mdBlad3) # (fast melt Hirudo verbana bladder)
prevdt.mdBlad = fm.mdBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
phyT.mdInt<-subset_samples(phyT.Md,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.mdInt = fast_melt(phyT.mdInt) # (fast melt Hirudo verbana intestinum)
prevdt.mdInt = fm.mdInt[, list(Prevalence = sum(count >= cMin),
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
phyT.mdILF<-subset_samples(phyT.Md,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.mdILF = fast_melt(phyT.mdILF) # (fast melt Hirudo verbana ILF)
prevdt.mdILF = fm.mdILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###CT
phyT.ctBlad<-subset_samples(phyT.mdBlad3,AnimalSource=="Wlot") # (Hirudo verbana bladder)
fm.ctBlad = fast_melt(phyT.ctBlad) # (fast melt Hirudo verbana bladder)
prevdt.ctBlad = fm.ctBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
phyT.ctInt<-subset_samples(phyT.mdCT,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.ctInt = fast_melt(phyT.ctInt) # (fast melt Hirudo verbana intestinum)
prevdt.ctInt = fm.ctInt[, list(Prevalence = sum(count >= cMin),
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
phyT.ctILF<-subset_samples(phyT.mdCT,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.ctILF = fast_melt(phyT.ctILF) # (fast melt Hirudo verbana ILF)
prevdt.ctILF = fm.ctILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###MA
phyT.maBlad<-subset_samples(phyT.mdBlad3,AnimalSource=="GrotonMA") # (Massachusetts bladder)
fm.maBlad = fast_melt(phyT.maBlad) # (fast melt Massachusetts bladder)
prevdt.maBlad = fm.maBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts bladder)
phyT.maInt<-subset_samples(phyT.mdMA,Sample_Type=="Intestinum") # (Massachusetts intestinum)
fm.maInt = fast_melt(phyT.maInt) # (fast melt Massachusetts intestinum)
prevdt.maInt = fm.maInt[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts intestinum)
phyT.maILF<-subset_samples(phyT.mdMA,Sample_Type=="ILF") # (Massachusetts ILF)
fm.maILF = fast_melt(phyT.maILF) # (fast melt Massachusetts ILF)
prevdt.maILF = fm.maILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts ILF)
###NY
phyT.nyBlad<-subset_samples(phyT.mdNY,Sample_Type=="Bladder") # (New York bladder)
fm.nyBlad = fast_melt(phyT.nyBlad) # (fast melt New York bladder)
prevdt.nyBlad = fm.nyBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York bladder)
#phyT.nyInt<-subset_samples(phyT.mdNY,Sample_Type=="Intestinum") # (New York intestinum)
#fm.nyInt = fast_melt(phyT.nyInt) # (fast melt New York intestinum)
#prevdt.nyInt = fm.nyInt[, list(Prevalence = sum(count >= cMin),TotalPer = sum(count),MinCount = min(count),MaxCount = max(count)),by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York intestinum)
phyT.nyILF<-subset_samples(phyT.mdNY,Sample_Type=="ILF") # (New York ILF)
fm.nyILF = fast_melt(phyT.nyILF) # (fast melt New York ILF)
prevdt.nyILF = fm.nyILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York ILF)
###VT
phyT.vtBlad<-subset_samples(phyT.mdVT,Sample_Type=="Bladder") # (Vermont bladder)
fm.vtBlad = fast_melt(phyT.vtBlad) # (fast melt Vermont bladder)
prevdt.vtBlad = fm.vtBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont bladder)
#phyT.vtInt<-subset_samples(phyT.mdVT,Sample_Type=="Intestinum") # (Vermont intestinum)
#fm.vtInt = fast_melt(phyT.vtInt) # (fast melt Vermont intestinum)
#prevdt.vtInt = fm.vtInt[, list(Prevalence = sum(count >= cMin),TotalPer = sum(count),MinCount = min(count),MaxCount = max(count)),by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont intestinum)
phyT.vtILF<-subset_samples(phyT.mdVT,Sample_Type=="ILF") # (Vermont ILF)
fm.vtILF = fast_melt(phyT.vtILF) # (fast melt Vermont ILF)
prevdt.vtILF = fm.vtILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont ILF)

##########################
##### List core OTUs #####
##########################
cpR<-c(.9) # define the level for determining core OTUs, reported as a percent of the samples tested (core percent)
### Calculated with 85% cut-off rather than 95% cut-off due to low sample size. 85% = 6/7 ###

### Hirudo verbana
ls.coreHvILF = prevdt.HvILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo ILF
ls.coreHvInt = prevdt.HvInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo Intestinum
ls.coreHvBlad = prevdt.HvBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo Bladder
ls.coreHv = unique(c(ls.coreHvBlad,ls.coreHvILF,ls.coreHvInt))
### CT Macrobdella decora
ls.coreCtILF = prevdt.ctILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.ctILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo ILF 
ls.coreCtInt = prevdt.ctInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.ctInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo Intestinum
ls.coreCtBlad = prevdt.ctBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBlad3,AnimalSource=="Wlot")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Hirudo Bladder
ls.coreCT = unique(c(ls.coreCtBlad,ls.coreCtILF,ls.coreCtInt))
### MA Macrobdella decora
ls.coreMaILF = prevdt.maILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.maILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts ILF (
ls.coreMaInt = prevdt.maInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.maInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts Intestinum
ls.coreMaBlad = prevdt.maBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBlad3,AnimalSource=="GrotonMA")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts Bladder
ls.coreMA = unique(c(ls.coreMaBlad,ls.coreMaILF,ls.coreMaInt))
### NY Macrobdella decora
ls.coreNyILF = prevdt.nyILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.nyILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella New York ILF
#ls.coreNyInt = prevdt.nyInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.nyInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella New York Intestinum
ls.coreNyBlad = prevdt.nyBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.nyBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella New York Bladder
ls.coreNY = unique(c(ls.coreNyBlad,ls.coreNyILF))
### VT Macrobdella decora
ls.coreVtILF = prevdt.vtILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.vtILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont ILF (
#ls.coreVtInt = prevdt.vtInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.vtInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont Intestinum
ls.coreVtBlad = prevdt.vtBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.vtBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont Bladder
ls.coreVT = unique(c(ls.coreVtBlad,ls.coreVtILF))
### Compiled Macrobdella decora
ls.coreMdILF = prevdt.mdILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreMdInt = prevdt.mdInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont Intestinum
ls.coreMdBlad = prevdt.mdBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBlad3,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont Bladder
ls.coreMd = unique(c(ls.coreCT,ls.coreMA,ls.coreNY)) # removed ls.coreVT due to this being only 1 sample
ls.coreFeedILF = unique(c(ls.coreCtILF,ls.coreMaILF))
ls.coreFeedInt = unique(c(ls.coreCtInt,ls.coreMaInt))
ls.coreTot = unique(c(ls.coreMd,ls.coreHv,ls.coreCT,ls.coreMA,ls.coreNY,ls.coreMdILF,ls.coreMdInt,ls.coreMdBlad)) # removed ls.coreVT due to this being only 1 sample

############################
##### List common OTUs #####
############################
cpM<-c(.75) # define the level for determining common OTUs, reported as a percent of the samples tested (core percent)

### Hirudo verbana ###
ls.comHvILF = prevdt.HvILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo ILF
ls.comHvInt = prevdt.HvInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo Intestinum
ls.comHvBlad = prevdt.HvBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo Bladder
ls.comHv = unique(c(ls.comHvBlad,ls.comHvILF,ls.comHvInt))
### CT Macrobdella decora ###
ls.comCtILF = prevdt.ctILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.ctILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo ILF 
ls.comCtInt = prevdt.ctInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.ctInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo Intestinum
ls.comCtBlad = prevdt.ctBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.ctBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Hirudo Bladder
ls.comCT = unique(c(ls.comCtBlad,ls.comCtILF,ls.comCtInt))
### MA Macrobdella decora ###
ls.comMaILF = prevdt.maILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.maILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of com OTUs for Macrobdella Massachusetts ILF (
ls.comMaInt = prevdt.maInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.maInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Massachusetts Intestinum
ls.comMaBlad = prevdt.maBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.maBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Massachusetts Bladder
ls.comMA = unique(c(ls.comMaBlad,ls.comMaILF,ls.comMaInt))
### NY Macrobdella decora ###
ls.comNyILF = prevdt.nyILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.nyILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella New York ILF
#ls.comNyInt = prevdt.nyInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.nyInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella New York Intestinum
ls.comNyBlad = prevdt.nyBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.nyBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella New York Bladder
ls.comNY = unique(c(ls.comNyBlad,ls.comNyILF))
### VT Macrobdella decora ###
ls.comVtILF = prevdt.vtILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.vtILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont ILF (
#ls.comVtInt = prevdt.vtInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.vtInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont Intestinum
ls.comVtBlad = prevdt.vtBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.vtBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont Bladder
ls.comVT = unique(c(ls.comVtBlad,ls.comVtILF))
### Compiled Macrobdella decora ###
ls.comMdILF = prevdt.mdILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comMdInt = prevdt.mdInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont Intestinum
ls.comMdBlad = prevdt.mdBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont Bladder
ls.comMd = unique(c(ls.comCT,ls.comMA,ls.comNY)) # removed ls.comVT due to this being only 1 sample
ls.comTot = unique(c(ls.comMd,ls.comHv,ls.comCT,ls.comMA,ls.comNY,ls.comMdILF,ls.comMdInt,ls.comMdBlad)) # removed ls.comVT due to this being only 1 sample

#################
##### Plots #####
#################

##### Core plot #####
phyT.Md2<-subset_samples(phyT.Md,Da1F%in%c("0","1","31"))
phyT.Md1<-subset_samples(phyT.Md2,Sample_Type!="Bladder")
phyT.Con<-merge_phyloseq(phyT.Md1,phyT.Hv) # merge Macrobdella and Hirudo phyloseq objects (physeq control)
phyT.base<-merge_phyloseq(phyT.Con,phyT.mdBlad)

########## Macrobdella Sites plot with 'Other' category ##########
mat.base <- sort(taxa_sums(phyT.base), TRUE)[1:35] # Identify 27 most abundant taxa
phyT.basePru <- prune_taxa(names(mat.base),phyT.base) # Create a subset of data including only 27 most abundant taxa
pBar.matCore <- plot_bar(phyT.basePru,x="Sample",fill="Family") + 
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") + 
  scale_fill_manual(values=pairBiome)
pBar.matCore

phyT.baseGI<-subset_samples(phyT.base,Sample_Type!="Bladder")
dt.core<-data.table(psmelt(phyT.baseGI))
dt.core$Genus3<-with(dt.core,
                    ifelse(Number%in%c(ls.coreTot),as.character(Genus2),
                           as.character("other")))

# Define levels for OTUs
dt.core$Genus3 <- factor(dt.core$Genus3, levels=c("Aeromonas","Aeromonas2",
                                                "Mucinivorans","Bacteroides","Millionella-like","Mucinivorans-like",
                                                "Alkaliphilus-like","Clostridium","Proteocatella-like","unk_Peptostreptococcaceae","Butyricicoccus","Papillibacter-like","unk_Ruminococcaceae",
                                                "Ochrobactrum","Ensifer","Rhizobium","Rhizobium-like",
                                                "Ramlibacter","Methylopumilus-like",
                                                "Bdellovibrio-like","Desulfovibrio","Cystobacter-like",
                                                "Azospirillum","Insolitispirillum-like","Phreatobacter-like",
                                                "other"))
# Define colors for OTUs. Reference = "rainbow"
genus.color<-c(Aeromonas="#267326",Aeromonas2="#39ac39",
               Mucinivorans="#B80000",Bacteroides="#F00000","Millionella-like"="#FF7777","Mucinivorans-like"="#ffcccc",
               "Alkaliphilus-like"="#000080",Clostridium="#0000cd","Proteocatella-like"="#8282ff","unk_Peptostreptococcaceae"="#cfcfff",Butyricicoccus="#0090b4","Papillibacter-like"="#01cdff",unk_Ruminococcaceae="#80e6ff",
               Ochrobactrum="#b34700",Ensifer="#ff6600",Rhizobium="#ff9933","Rhizobium-like"="#ffcc99",
               Ramlibacter="#600080","Methylopumilus-like"="#ac00e6",
               "Bdellovibrio-like"="#cccc00",Desulfovibrio="#ffff00","Cystobacter-like"="#ffffb3",
               Azospirillum="#9a0066","Insolitispirillum-like"="#ff4ec5","Phreatobacter-like"="#ff9ade",
               other="#808080")
##### Plot #####
pBar.Core <- ggplot(dt.core, aes(x=Replicate,y=Abundance, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.text.x=element_blank(),axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pBar.Core # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.Core), size = "last")), filename="Plots/plotBarStack_Core.png", width=12,height=8)
##### Table #####
write.table(tax_table(phyT.base), "tableTax_Core.csv", sep=",")
write.table(dt.core,"table_dtCore.csv",sep=",")

##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.mdBlad), size = "last")), filename="Plots/plotBarStack_CoreBladderOrder.png", width=12,height=8)


########## Evaluation of non-core taxa ##########
phyT.nonCore<-subset_taxa(phyT.baseGI,!Number%in%c(ls.coreTot))
fm.nonCore = fast_melt(phyT.nonCore) # (fast melt non core)
prevdt.nonCore = fm.nonCore[, list(Prevalence = sum(count >= cMin),
                                   TotalPer = sum(count),
                                   MinCount = min(count),
                                   MaxCount = max(count)),
                            by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
ls.hiNonCore = prevdt.nonCore[(MaxCount >= .01), TaxaID] # Make list of OTUs present in more than 2 samples and having at least maxNeg reads in at least one sample (list.high prevalence)
phyT.pruNonCore <- prune_taxa(ls.hiNonCore,phyT.nonCore) # remove identified taxa from phyloseq object (physeq.Pruned 2)
##### Plot #####
pBar.nonCore<-plot_bar(phyT.pruNonCore,x="Replicate",fill="Genus2") + 
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") + 
  theme(text=element_text(size=10), axis.text.x=element_blank(),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.nonCore
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.nonCore), size = "last")), filename="Plots/plotBarStack_nonCore.png", width=12,height=8)
##### Table #####
write.table(tax_table(phyT.pruNonCore), "tableTax_nonCore.csv", sep=",")

########## Print table of core taxa with columns indicating which sample type each is present in ##########
### Add core columns to taxa table ###
phyT.coreTot<-prune_taxa(ls.coreTot,phyT.base) # keep only taxa from the identified 'Core' 
tax.core<-tax_table(phyT.coreTot) # pull taxonomy data from phyloseq object
taxdt.Core<-as.data.table(tax.core) # Change taxa_table to a data.table
taxdt.Core$Number<-as.character(taxdt.Core$Number) # force TAXcore$Number to be a character

# Macrobdella core data
taxdt.mdInt<-as.data.table(ls.coreMdInt)
taxdt.mdInt$Md.Intestinum<-"x"
setnames(taxdt.mdInt,"ls.coreMdInt", "Number") # First column = OTU names. Change column heading
taxdt.mdIntCT<-as.data.table(ls.coreCtInt)
taxdt.mdIntCT$Md.Intestinum.CT<-"x"
setnames(taxdt.mdIntCT,"ls.coreCtInt", "Number") # First column = OTU names. Change column heading
taxdt.mdIntMA<-as.data.table(ls.coreMaInt)
taxdt.mdIntMA$Md.Intestinum.MA<-"x"
setnames(taxdt.mdIntMA,"ls.coreMaInt", "Number") # First column = OTU names. Change column heading

taxdt.mdILF<-as.data.table(ls.coreMdILF)
taxdt.mdILF$Md.ILF<-"x"
setnames(taxdt.mdILF,"ls.coreMdILF", "Number") # First column = OTU names. Change column heading
taxdt.mdILFct<-as.data.table(ls.coreCtILF)
taxdt.mdILFct$Md.ILF.CT<-"x"
setnames(taxdt.mdILFct,"ls.coreCtILF", "Number") # First column = OTU names. Change column heading
taxdt.mdILFma<-as.data.table(ls.coreMaILF)
taxdt.mdILFma$Md.ILF.MA<-"x"
setnames(taxdt.mdILFma,"ls.coreMaILF", "Number") # First column = OTU names. Change column heading

taxdt.mdBlad<-as.data.table(ls.coreMdBlad)
taxdt.mdBlad$Md.Bladder<-"x"
setnames(taxdt.mdBlad,"ls.coreMdBlad", "Number") # First column = OTU names. Change column heading

#Hirudo core data
taxdt.hvInt<-as.data.table(ls.coreHvInt)
taxdt.hvInt$Hv.Intestinum<-"x"
setnames(taxdt.hvInt,"ls.coreHvInt", "Number") # First column = OTU names. Change column heading

taxdt.hvILF<-as.data.table(ls.coreHvILF)
taxdt.hvILF$Hv.ILF<-"x"
setnames(taxdt.hvILF,"ls.coreHvILF", "Number") # First column = OTU names. Change column heading

taxdt.hvBlad<-as.data.table(ls.coreHvBlad)
taxdt.hvBlad$Hv.Bladder<-"x"
setnames(taxdt.hvBlad,"ls.coreHvBlad", "Number") # First column = OTU names. Change 

# Merge tables
ls.taxdt <- list(taxdt.Core,taxdt.mdILF,taxdt.mdILFct,taxdt.mdILFma,taxdt.mdInt,taxdt.mdIntCT,taxdt.mdIntMA,taxdt.mdBlad,taxdt.hvILF,taxdt.hvInt,taxdt.hvBlad) # list of data.tables to merge (list of taxonomy data.tables)
lapply(ls.taxdt, function(i) setkey(i, Number)) # set key for merge function
dt.TAXcore2 <- Reduce(function(...) merge(..., all = T), ls.taxdt) # merge list of data.tables, keeping values even if not present in each table (Taxonomy table core 2)
dt.TAXcore2[is.na(dt.TAXcore2)] <- "" # replace <NA> values with empty values
write.table(dt.TAXcore2, "tableTax_CoreCompile.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded


########## Print table of common taxa with columns indicating which sample type each is present in ##########
### Add core columns to taxa table ###
phyT.comTot<-prune_taxa(ls.comTot,phyT.base) # keep only taxa from the identified 'Core' 
tax.com<-tax_table(phyT.comTot) # pull taxonomy data from phyloseq object
taxdt.Com<-as.data.table(tax.com) # Change taxa_table to a data.table
taxdt.Com$Number<-as.character(taxdt.Com$Number) # force TAXcom$Number to be a character

# Macrobdella common data
taxdt.mdInt<-as.data.table(ls.comMdInt)
taxdt.mdInt$Md.Intestinum<-"x"
setnames(taxdt.mdInt,"ls.comMdInt", "Number") # First column = OTU names. Change column heading

taxdt.mdILF<-as.data.table(ls.comMdILF)
taxdt.mdILF$Md.ILF<-"x"
setnames(taxdt.mdILF,"ls.comMdILF", "Number") # First column = OTU names. Change column heading

taxdt.mdBlad<-as.data.table(ls.comMdBlad)
taxdt.mdBlad$Md.Bladder<-"x"
setnames(taxdt.mdBlad,"ls.comMdBlad", "Number") # First column = OTU names. Change column heading

#Hirudo com data
taxdt.hvInt<-as.data.table(ls.comHvInt)
taxdt.hvInt$Hv.Intestinum<-"x"
setnames(taxdt.hvInt,"ls.comHvInt", "Number") # First column = OTU names. Change column heading

taxdt.hvILF<-as.data.table(ls.comHvILF)
taxdt.hvILF$Hv.ILF<-"x"
setnames(taxdt.hvILF,"ls.comHvILF", "Number") # First column = OTU names. Change column heading

taxdt.hvBlad<-as.data.table(ls.comHvBlad)
taxdt.hvBlad$Hv.Bladder<-"x"
setnames(taxdt.hvBlad,"ls.comHvBlad", "Number") # First column = OTU names. Change 

# Merge tables
ls.taxdt <- list(dt.taxcom,taxdt.mdILF,taxdt.mdInt,taxdt.mdBlad,taxdt.hvILF,taxdt.hvInt,taxdt.hvBlad) # list of data.tables to merge (list of taxonomy data.tables)
lapply(ls.taxdt, function(i) setkey(i, Number)) # set key for merge function
dt.TAXcom2 <- Reduce(function(...) merge(..., all = T), ls.taxdt) # merge list of data.tables, keeping values even if not present in each table (Taxonomy table core 2)
dt.TAXcom2[is.na(dt.TAXcom2)] <- "" # replace <NA> values with empty values
write.table(dt.TAXcom2, "tableTax_comCompile.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded
### Add core columns to taxa table ###


###############################################
################ Practice area ################ 
############################################### 

