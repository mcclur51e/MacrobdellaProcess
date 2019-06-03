# Thank you to jeffkimbrel for the next two lines! #
ls.taxaDC <- setdiff(taxa_names(phyT.Tform), taxa_names(phyR.allc)) # find the difference between taxa in contaminant list and physeqAn (taxa decontam)
phyT.DC <- prune_taxa(ls.taxaDC, phyT.Tform) # keep only taxa not in taxaNC (physeq decontam)

phyT.adult<-subset_samples(phyT.DC,Age=="A") # keep adult samples (prune Adult)
phyT.Nelson<-subset_samples(phyT.adult,SourceCollector=="NelsonM")
phyT.hvNel<-subset_samples(phyT.Nelson,Taxonomic_ID=="Hverbana")
phyT.USA<-subset_samples(phyT.adult,AnimalSource=="USA" & Subproject!="Deposit" & Da2F%in%c("none","0"))
phyT.hv2<-subset_samples(phyT.DC,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042317EMe.u","A042317EMf.u","A042317EMh.u","A042117JGbb","A042317EMi","A042317EMj","A042317EMj.b","A042317EMj.u","A042317EMk","A042317EMk.b","A043017EMm.u","A42hF050317EMb","A050217EMo","A050217EMo.u","A042117JGb.b","A042117JGb","A042317EMi.b","A042317EMi","A043017EMl.b","A050217EMq"))
phyT.hv<-merge_phyloseq(phyT.USA,phyT.hv2,phyT.hvNel)
phyT.md<-subset_samples(phyT.adult,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
phyT.mdCT<-subset_samples(phyT.md,AnimalSource=="Wlot") # subset containing only W-lot samples (CT Macrobdella)
phyT.mdMA<-subset_samples(phyT.md,AnimalSource=="GrotonMA") # subset containing only MA samples (MA Macrobdella)
phyT.mdNY<-subset_samples(phyT.md,AnimalSource=="CarogaNY") # subset containing only NY samples (NY Macrobdella)
phyT.mdVT<-subset_samples(phyT.md,AnimalSource=="MtSnowVT") # subset containing only VT samples (VT Macrobdella)
#phyT.mdSB<-subset_samples(phyT.md,AnimalSource=="Schoolhouse_Brook") # subset containing only Schoolhouse Brook samples (VT Macrobdella)
phyT.leech<-merge_phyloseq(phyT.md,phyT.hv) # merge Macrobdella and Hirudo phyloseq objects (physeq control)

dt.leech = fast_melt(phyT.leech) # make data table from physeq object (data table. physeq Pruned)
prev.leech = dt.leech[, list(Prevalence = sum(count > 0), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                             by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.Pres2 = prev.leech[(Prevalence > 0 & TotalPer >= .001), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)
phyT.start <- prune_taxa(ls.Pres2,phyT.leech) # remove identified taxa from phyloseq object (physeq.Pruned 2)


cMin<-0.001 # define minimum to count prevalence in a sample
### Hirudo ###
phyT.hvBlad<-subset_samples(phyT.hv,Sample_Type=="Bladder") # (Hirudo verbana bladder)
fm.hvBlad = fast_melt(phyT.hvBlad) # (fast melt Hirudo verbana bladder)
prevdt.hvBlad = fm.hvBlad[,list(Prevalence = sum(count >= cMin), 
  TotalCount = sum(count),
  MinCount = min(count),
  MaxCount = max(count),
  MedCount = median(count),
  Med2Count = median(count[count!=0]),
  PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
colnames(prevdt.hvBlad)[colnames(prevdt.hvBlad)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.hvBlad,tax_table(phyT.hvBlad)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.hvBlad <- Reduce(function(...) merge(..., all = T), list(prevdt.hvBlad,tax_table(phyT.hvBlad)[,c("Number","Genus2")])) # m     
prevdt.hvBladhi<-prevdt.hvBlad[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.hvBladhi, "tablePrev_hvBlad.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.hvInt<-subset_samples(phyT.hv,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.hvInt = fast_melt(phyT.hvInt) # (fast melt Hirudo verbana intestinum)
prevdt.hvInt = fm.hvInt[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.hvInt)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
colnames(prevdt.hvInt)[colnames(prevdt.hvInt)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.hvInt,tax_table(phyT.hvInt)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.hvInt <- Reduce(function(...) merge(..., all = T), list(prevdt.hvInt,tax_table(phyT.hvInt)[,c("Number","Genus2")])) # m     
prevdt.hvInthi<-prevdt.hvInt[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.hvInthi, "tablePrev_hvInt.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.hvILF<-subset_samples(phyT.hv,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.hvILF = fast_melt(phyT.hvILF) # (fast melt Hirudo verbana ILF)
prevdt.hvILF = fm.hvILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.hvILF)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.hvILF)[colnames(prevdt.hvILF)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.hvILF,tax_table(phyT.hvILF)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.hvILF <- Reduce(function(...) merge(..., all = T), list(prevdt.hvILF,tax_table(phyT.hvILF)[,c("Number","Genus2")])) # m     
prevdt.hvILFhi<-prevdt.hvILF[which(PrevPer>0),][order(-PrevPer)]
#write.table(prevdt.hvILFhi, "tableTax_hvILF.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.hvILF0<-subset_samples(phyT.hvILF,Da1Fb%in%c("0")) # (Hirudo verbana ILF)
fm.hvILF0 = fast_melt(phyT.hvILF0) # (fast melt Hirudo verbana ILF)
prevdt.hvILF0 = fm.hvILF0[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.hvILF0)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.hvILF0)[colnames(prevdt.hvILF0)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.hvILF0,tax_table(phyT.hvILF0)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.hvILF0 <- Reduce(function(...) merge(..., all = T), list(prevdt.hvILF0,tax_table(phyT.hvILF0)[,c("Number","Genus2")])) # m     
prevdt.hvILF0hi<-prevdt.hvILF0[which(PrevPer>0),][order(-PrevPer)]
#write.table(prevdt.hvILF0hi, "tableTax_hvILF0.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

###Macrobdella		
phyT.mdBlad<-subset_samples(phyT.md,Sample_Type=="Bladder") # (Hirudo verbana bladder)
phyT.mdBlad2<-subset_samples(phyT.mdBlad,!Da1F%in%c("1","4","7"))
phyT.mdBlad3<-subset_samples(phyT.mdBlad2,!sample_names(phyT.mdBlad2)%in%c("MdMA0dF040618EMfBD","MdMA0dF040618EMdBD722","MdMA0dF040618EMeBD"))
fm.mdBlad = fast_melt(phyT.mdBlad3) # (fast melt Hirudo verbana bladder)
prevdt.mdBlad = fm.mdBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.mdBlad3)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana Blad)
colnames(prevdt.mdBlad)[colnames(prevdt.mdBlad)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdBlad,tax_table(phyT.mdBlad)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdBlad <- Reduce(function(...) merge(..., all = T), list(prevdt.mdBlad,tax_table(phyT.mdBlad)[,c("Number","Genus2")])) # m     
prevdt.mdBladhi<-prevdt.mdBlad[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.mdBladhi, "tablePrev_mdBlad.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.mdInt<-subset_samples(phyT.md,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.mdInt = fast_melt(phyT.mdInt) # (fast melt Hirudo verbana intestinum)
prevdt.mdInt = fm.mdInt[, list(Prevalence = sum(count >= cMin),
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.mdInt)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.mdInt)[colnames(prevdt.mdInt)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdInt,tax_table(phyT.mdInt)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdInt <- Reduce(function(...) merge(..., all = T), list(prevdt.mdInt,tax_table(phyT.mdInt)[,c("Number","Genus2")])) # m     
prevdt.mdInthi<-prevdt.mdInt[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.mdInthi, "tablePrev_mdInt.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.mdILF<-subset_samples(phyT.md,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.mdILF = fast_melt(phyT.mdILF) # (fast melt Hirudo verbana ILF)
prevdt.mdILF = fm.mdILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.mdILF)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.mdILF)[colnames(prevdt.mdILF)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdILF,tax_table(phyT.mdILF)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdILF <- Reduce(function(...) merge(..., all = T), list(prevdt.mdILF,tax_table(phyT.mdILF)[,c("Number","Genus2")])) # m     
prevdt.mdILFhi<-prevdt.mdILF[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.mdILFhi, "tablePrev_mdILF.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.mdILF0<-subset_samples(phyT.mdILF,Da1Fb=="0") # (Hirudo verbana ILF)
fm.mdILF0 = fast_melt(phyT.mdILF0) # (fast melt Hirudo verbana ILF)
prevdt.mdILF0 = fm.mdILF0[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   
  Med2Count = median(count[count!=0]),   
  PrevPer = sum(count >= cMin) / nsamples(phyT.mdILF0)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.mdILF0)[colnames(prevdt.mdILF0)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdILF0,tax_table(phyT.mdILF0)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdILF0 <- Reduce(function(...) merge(..., all = T), list(prevdt.mdILF0,tax_table(phyT.mdILF0)[,c("Number","Genus2")])) # m     
prevdt.mdILF0hi<-prevdt.mdILF0[which(PrevPer>0),][order(-PrevPer)]

###CT
phyT.ctBlad<-subset_samples(phyT.mdBlad3,AnimalSource=="Wlot") # (Hirudo verbana bladder)
fm.ctBlad = fast_melt(phyT.ctBlad) # (fast melt Hirudo verbana bladder)
prevdt.ctBlad = fm.ctBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
phyT.ctInt<-subset_samples(phyT.mdCT,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.ctInt = fast_melt(phyT.ctInt) # (fast melt Hirudo verbana intestinum)
prevdt.ctInt = fm.ctInt[, list(Prevalence = sum(count >= cMin),
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
phyT.ctILF<-subset_samples(phyT.mdCT,Sample_Type=="ILF") # (Hirudo verbana ILF)
fm.ctILF = fast_melt(phyT.ctILF) # (fast melt Hirudo verbana ILF)
prevdt.ctILF = fm.ctILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###MA
phyT.maBlad<-subset_samples(phyT.mdBlad3,AnimalSource=="GrotonMA") # (Massachusetts bladder)
fm.maBlad = fast_melt(phyT.maBlad) # (fast melt Massachusetts bladder)
prevdt.maBlad = fm.maBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts bladder)
phyT.maInt<-subset_samples(phyT.mdMA,Sample_Type=="Intestinum") # (Massachusetts intestinum)
fm.maInt = fast_melt(phyT.maInt) # (fast melt Massachusetts intestinum)
prevdt.maInt = fm.maInt[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts intestinum)
phyT.maILF<-subset_samples(phyT.mdMA,Sample_Type=="ILF") # (Massachusetts ILF)
fm.maILF = fast_melt(phyT.maILF) # (fast melt Massachusetts ILF)
prevdt.maILF = fm.maILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts ILF)
###NY
phyT.nyBlad<-subset_samples(phyT.mdNY,Sample_Type=="Bladder") # (New York bladder)
fm.nyBlad = fast_melt(phyT.nyBlad) # (fast melt New York bladder)
prevdt.nyBlad = fm.nyBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York bladder)
#phyT.nyInt<-subset_samples(phyT.mdNY,Sample_Type=="Intestinum") # (New York intestinum)
#fm.nyInt = fast_melt(phyT.nyInt) # (fast melt New York intestinum)
#prevdt.nyInt = fm.nyInt[, list(Prevalence = sum(count >= cMin),TotalPer = sum(count),MinCount = min(count),MaxCount = max(count),   MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York intestinum)
phyT.nyILF<-subset_samples(phyT.mdNY,Sample_Type=="ILF") # (New York ILF)
fm.nyILF = fast_melt(phyT.nyILF) # (fast melt New York ILF)
prevdt.nyILF = fm.nyILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York ILF)
###VT
phyT.vtBlad<-subset_samples(phyT.mdVT,Sample_Type=="Bladder") # (Vermont bladder)
fm.vtBlad = fast_melt(phyT.vtBlad) # (fast melt Vermont bladder)
prevdt.vtBlad = fm.vtBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont bladder)
#phyT.vtInt<-subset_samples(phyT.mdVT,Sample_Type=="Intestinum") # (Vermont intestinum)
#fm.vtInt = fast_melt(phyT.vtInt) # (fast melt Vermont intestinum)
#prevdt.vtInt = fm.vtInt[, list(Prevalence = sum(count >= cMin),TotalPer = sum(count),MinCount = min(count),MaxCount = max(count),   MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont intestinum)
phyT.vtILF<-subset_samples(phyT.mdVT,Sample_Type=="ILF") # (Vermont ILF)
fm.vtILF = fast_melt(phyT.vtILF) # (fast melt Vermont ILF)
prevdt.vtILF = fm.vtILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count),   
  MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont ILF)

##########################
##### List core OTUs #####
##########################
cpR<-c(.9) # define the level for determining core OTUs, reported as a percent of the samples tested (core percent)
### Calculated with 85% cut-off rather than 95% cut-off due to low sample size. 85% = 6/7 ###

### Hirudo verbana
ls.coreHvILF0 = prevdt.hvILF0[(Prevalence >= cpR*nsamples(phyT.hvILF0) & MaxCount >= cMin), Number] # Make list of core OTUs for Hirudo ILF
ls.coreHvILF = prevdt.hvILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of core OTUs for Hirudo ILF
ls.coreHvInt = prevdt.hvInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of core OTUs for Hirudo Intestinum
ls.coreHvBlad = prevdt.hvBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.hvBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of core OTUs for Hirudo Bladder
ls.coreHv = unique(c(ls.coreHvBlad,ls.coreHvILF,ls.coreHvInt,ls.coreHvILF0))
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
ls.coreMdILF = prevdt.mdILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreMdILF0 = prevdt.mdILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdILF,Da1Fb=="0")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreMdInt = prevdt.mdInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont Intestinum
ls.coreMdBlad = prevdt.mdBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBlad3,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont Bladder
ls.coreMd = unique(c(ls.coreCT,ls.coreMA,ls.coreNY)) # removed ls.coreVT due to this being only 1 sample
ls.coreFeedILF = unique(c(ls.coreCtILF,ls.coreMaILF))
ls.coreFeedInt = unique(c(ls.coreCtInt,ls.coreMaInt))
ls.coreTot = unique(c(ls.coreMd,ls.coreHv,ls.coreCT,ls.coreMA,ls.coreNY,ls.coreMdILF,ls.coreMdILF0,ls.coreMdInt,ls.coreMdBlad)) # removed ls.coreVT due to this being only 1 sample

############################
##### List common OTUs #####
############################
cpM<-c(.7) # define the level for determining common OTUs, reported as a percent of the samples tested (core percent)

### Hirudo verbana ###
ls.comHvILF = prevdt.hvILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of common OTUs for Hirudo ILF
ls.comHvILF0 = prevdt.hvILF0[(Prevalence >= cpM*nsamples(phyT.hvILF0) & MaxCount >= cMin), Number] # Make list of core OTUs for Hirudo ILF
ls.comHvInt = prevdt.hvInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of common OTUs for Hirudo Intestinum
ls.comHvBlad = prevdt.hvBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.hvBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of common OTUs for Hirudo Bladder
ls.comHv = unique(c(ls.comHvBlad,ls.comHvILF,ls.comHvInt,ls.comHvILF0))
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
ls.comMdILF = prevdt.mdILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comMdILF0 = prevdt.mdILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdILF,Da1Fb=="0")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comMdInt = prevdt.mdInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdInt,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont Intestinum
ls.comMdBlad = prevdt.mdBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdBlad3,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont Bladder
ls.comMd = unique(c(ls.comCT,ls.comMA,ls.comNY)) # removed ls.comVT due to this being only 1 sample
ls.comTot = unique(c(ls.comMd,ls.comHv,ls.comCT,ls.comMA,ls.comNY,ls.comMdILF,ls.comMdILF0,ls.comMdInt,ls.comMdBlad)) # removed ls.comVT due to this being only 1 sample

#################
##### Plots #####
#################

##### Core plot #####
phyT.md2<-subset_samples(phyT.md,Da1F%in%c("0","1","31"))
phyT.md1<-subset_samples(phyT.md2,Sample_Type!="Bladder")
phyT.hv0<-subset_samples(subset_samples(phyT.hv,Da1F%in%c("0")),Sample_Type=="ILF")
phyT.Con<-merge_phyloseq(phyT.md1,phyT.hv0,phyT.hvInt) # merge Macrobdella and Hirudo phyloseq objects (physeq control)
phyT.base<-merge_phyloseq(phyT.Con,phyT.mdBlad,phyT.hvBlad)

########## Macrobdella Sites plot with 'Other' category ##########
mat.base <- sort(taxa_sums(phyT.base), TRUE)[1:35] # Identify 27 most abundant taxa
phyTmat.base <- prune_taxa(names(mat.base),phyT.base) # Create a subset of data including only 27 most abundant taxa
pBar.matBase <- plot_bar(phyTmat.base,x="Sample",fill="Genus2") + 
  facet_grid(Sample_Type~Taxonomic_ID, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pal.pairBiome)
pBar.matBase

phyT.baseGI<-subset_samples(phyT.base,Sample_Type!="Bladder")
dt.GI<-data.table(psmelt(phyT.baseGI))
dt.GI$Genus3<-with(dt.GI,
                   ifelse(Number%in%c(ls.coreTot),as.character(Genus2),
                   as.character("other")))

# Define levels for OTUs
dt.GI$Genus3 <- factor(dt.GI$Genus3, levels=c("Aeromonas","Aeromonas2",
                                              "Mucinivorans","Bacteroides","Bacteroides-like","Millionella-like",
                                              "Alkaliphilus-like","Clostridium","Clostridium-like","Proteocatella","Proteocatella-like","unk_Peptostreptococcaceae","Butyricicoccus","Butyricicoccus-like","Papillibacter-like","unk_Ruminococcaceae",
                                              "Ochrobactrum","Ensifer","Rhizobium","Rhizobium-like",
                                              "Ramlibacter","Methylopumilus-like",
                                              "Bdellovibrio-like","Desulfovibrio","Cystobacter-like","Pedobacter",
                                              "Azospirillum","Insolitispirillum-like","Phreatobacter-like",
                                              "other"))
# Define colors for OTUs. Reference = "rainbow"
genus.color<-c(Aeromonas="#267326",Aeromonas2="#39ac39",
               Mucinivorans="#B80000",Bacteroides="#F00000","Bacteroides-like"="#FF7777","Millionella-like"="#ffcccc",
               "Alkaliphilus-like"="#000080",Clostridium="#0000cd","Clostridium-like"="#8282ff",Proteocatella="#cfcfff","Proteocatella-like"="#005f6c",unk_Peptostreptococcaceae="#00a4bb",Butyricicoccus="#1ee3ff","Butyricicoccus-like"="#bbf7ff",
               Ochrobactrum="#623800",Ensifer="#c47000",Rhizobium="#ff9914","Rhizobium-like"="#ffddb1",
               Ramlibacter="#430059","Methylopumilus-like"="#7d00a7","Papillibacter-like"="#cc32ff","unk_Ruminococcaceae"="#eebbff",
               "Bdellovibrio-like"="#626200",Desulfovibrio="#c4c400","Cystobacter-like"="#ffff14","Pedobacter"="#ffffb1",
               Azospirillum="#750063","Insolitispirillum-like"="#c400a5","Phreatobacter-like"="#ff13da",
               other="#808080") 
##### Plot #####
pBar.GI <- ggplot(dt.GI, aes(x=Replicate,y=Abundance, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pBar.GI # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.GI), size = "last")), filename="Plots/plotBarStack_Core.eps", device="eps", width=12,height=8)
##### Table #####
write.table(tax_table(phyT.base), "tableTax_Core.csv", sep=",")
write.table(dt.GI,"table_GI.csv",sep=",")

##### Plot by Order #####
phyT.baseCore<-subset_taxa(phyT.base,Number%in%c(ls.comTot))
phyT.baseGI<-subset_samples(phyT.base,Sample_Type!="Bladder")
phyT.baseGI<-phyT.base
fm.baseGI = fast_melt(phyT.baseGI) # make data table from physeq object (data table. physeq Pruned)
prev.baseGI = fm.baseGI[, list(Prevalence = sum(count>=.001), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
                               by = Order] # make simple table l
ls.order<-sort(prev.baseGI[(Prevalence > 2), Order], TRUE)[1:40] # this isn't actually sorting by anything I can detect

dt.baseHi<-data.table(psmelt(phyT.baseGI))
dt.baseHi$Order2<-with(dt.baseHi,
                     ifelse(Order%in%c(ls.order),as.character(Order),
                            as.character("Z_other")))

##### Plot #####
pBar.CoreOrder <- ggplot(dt.baseHi, aes(x=Replicate, y=Abundance, fill=Order2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.CoreOrder # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.CoreOrder), size = "last")), filename="Plots/plotBarStack_CoreOrder.eps", device="eps", width=12,height=8)

###### Bladder #####
ls.coreBlad<-unique(c(ls.coreMdBlad,ls.coreHvBlad)) # 
phyT.baseB<-subset_samples(phyT.base,Sample_Type=="Bladder")
phyT.baseBmd<-subset_samples(phyT.baseB,Taxonomic_ID=="Mdecora")
dt.mdBlad<-data.table(psmelt(phyT.baseBmd))
dt.mdBlad$Genus3<-with(dt.mdBlad,
                     ifelse(Number%in%c(ls.coreBlad),as.character(Genus2),
                            as.character("other")))
dt.mdBlad$Order2<-with(dt.mdBlad,
                     ifelse(Number%in%c(ls.coreBlad),as.character(Order),
                            as.character("other")))
# Define levels for OTUs
dt.mdBlad$Genus3 <- factor(dt.core$Genus3, levels=c("Aeromonas","Proteus",
                                                    "Mucinivorans","unk_Bacteroidaceae","Bacteroides","unk_Rikenellaceae",
                                                    "Proteocatella","unk_Peptostreptococcaceae","Butyricicoccus","unk_Ruminococcaceae","unk_Christensenellaceae","unk_Clostridiaceae","Proteiniclasticum","Clostridium",
                                                    "Ochrobactrum","unk_Rhodospirillaceae","unk_Rhodospirillales","unk_Rhizobiales","Ensifer",
                                                    "Ramlibacter","unk_Methylophilaceae","Pedobacter",
                                                    "unk_Bdellovibrionaceae","Desulfovibrio","unk_Desulfovibrionaceae","unk_Cystobacteraceae",
                                                    "Flavobacterium","unk_Spirochaeta",
                                                    "Fusobacterium","unk_Fusobacteriaceae",
                                                    "other"))
# Define colors for OTUs. Reference = "rainbow"
genus.color<-c(Aeromonas="#267326",Proteus="#39ac39",
               Mucinivorans="#B80000",unk_Bacteroidaceae="#F00000",Bacteroides="#FF7777",unk_Rikenellaceae="#ffcccc",
               Proteocatella="#000080",unk_Peptostreptococcaceae="#0000cd",Butyricicoccus="#8282ff",unk_Ruminococcaceae="#cfcfff",unk_Christensenellaceae="#0090b4",unk_Clostridiaceae="#01cdff",Proteiniclasticum="#80e6ff",Clostridium="#ccf5ff",
               Ochrobactrum="#b34700",unk_Rhodospirillaceae="#ff6600",unk_Rhodospirillales="#ff9933",unk_Rhizobiales="#ffcc99",Ensifer="#ffe6cc",
               unk_Bdellovibrionaceae="#600080",Desulfovibrio="#ac00e6",unk_Desulfovibrionaceae="#e599ff",unk_Cystobacteraceae="#f2ccff",
               Ramlibacter="#cccc00",unk_Methylophilaceae="#ffff00",Pedobacter="#ffffb3",
               Flavobacterium="#9a0066",unk_Spirochaeta="#ff4ec5",
               Fusobacterium="#ff9ade",unk_Fusobacteriaceae="#ffe6f7",
               other="#808080")   
##### Plot #####
pBar.mdBlad <- ggplot(dt.mdBlad, aes(x=Sample,y=Abundance, fill=Order)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(~Header, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.mdBlad # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.mdBlad), size = "last")), filename="Plots/plotBarStack_CoreBladderOrder.eps", device="eps", width=12,height=8)
##### Table #####

########## Evaluation of non-core taxa ##########
phyT.nonCore<-subset_taxa(phyT.baseGI,!Number%in%c(ls.coreTot))
fm.nonCore = fast_melt(phyT.nonCore) # (fast melt non core)
prevdt.nonCore = fm.nonCore[, list(Prevalence = sum(count >= cMin),
                                   TotalPer = sum(count),
                                   MinCount = min(count),
                                   MaxCount = max(count),   MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
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
ggsave(grid.draw(rbind(ggplotGrob(pBar.nonCore), size = "last")), filename="Plots/plotBarStack_nonCore.eps", device="eps", width=12,height=8)
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
taxdt.Com.mdInt<-as.data.table(ls.comMdInt)
taxdt.Com.mdInt$Md.Intestinum<-"x"
setnames(taxdt.Com.mdInt,"ls.comMdInt", "Number") # First column = OTU names. Change column heading
taxdt.Com.mdIntCT<-as.data.table(ls.comCtInt)
taxdt.Com.mdIntCT$Md.Intestinum.CT<-"x"
setnames(taxdt.Com.mdIntCT,"ls.comCtInt", "Number") # First column = OTU names. Change column heading
taxdt.Com.mdIntMA<-as.data.table(ls.comMaInt)
taxdt.Com.mdIntMA$Md.Intestinum.MA<-"x"
setnames(taxdt.Com.mdIntMA,"ls.comMaInt", "Number") # First column = OTU names. Change column heading

taxdt.Com.mdILF<-as.data.table(ls.comMdILF)
taxdt.Com.mdILF$Md.ILF<-"x"
setnames(taxdt.Com.mdILF,"ls.comMdILF", "Number") # First column = OTU names. Change column heading
taxdt.Com.mdILFct<-as.data.table(ls.comCtILF)
taxdt.Com.mdILFct$Md.ILF.CT<-"x"
setnames(taxdt.Com.mdILFct,"ls.comCtILF", "Number") # First column = OTU names. Change column heading
taxdt.Com.mdILFma<-as.data.table(ls.comMaILF)
taxdt.Com.mdILFma$Md.ILF.MA<-"x"
setnames(taxdt.Com.mdILFma,"ls.comMaILF", "Number") # First column = OTU names. Change column heading

taxdt.Com.mdBlad<-as.data.table(ls.comMdBlad)
taxdt.Com.mdBlad$Md.Bladder<-"x"
setnames(taxdt.Com.mdBlad,"ls.comMdBlad", "Number") # First column = OTU names. Change column heading

#Hirudo com data
taxdt.Com.hvInt<-as.data.table(ls.comHvInt)
taxdt.Com.hvInt$Hv.Intestinum<-"x"
setnames(taxdt.Com.hvInt,"ls.comHvInt", "Number") # First column = OTU names. Change column heading

taxdt.Com.hvILF<-as.data.table(ls.comHvILF0)
taxdt.Com.hvILF$Hv.ILF<-"x"
setnames(taxdt.Com.hvILF,"ls.comHvILF0", "Number") # First column = OTU names. Change column heading

taxdt.Com.hvBlad<-as.data.table(ls.comHvBlad)
taxdt.Com.hvBlad$Hv.Bladder<-"x"
setnames(taxdt.Com.hvBlad,"ls.comHvBlad", "Number") # First column = OTU names. Change 

# Merge tables
ls.taxdt.Com <- list(taxdt.Com,taxdt.Com.mdILF,taxdt.Com.mdILFct,taxdt.Com.mdILFma,taxdt.Com.mdInt,taxdt.Com.mdIntCT,taxdt.Com.mdIntMA,taxdt.Com.mdBlad,taxdt.Com.hvILF,taxdt.Com.hvInt,taxdt.Com.hvBlad) # list of data.tables to merge (list of taxonomy data.tables)
lapply(ls.taxdt.Com, function(i) setkey(i, Number)) # set key for merge function
dt.TAXcom2 <- Reduce(function(...) merge(..., all = T), ls.taxdt.Com) # merge list of data.tables, keeping values even if not present in each table (Taxonomy table core 2)
dt.TAXcom2[is.na(dt.TAXcom2)] <- "" # replace <NA> values with empty values
write.table(dt.TAXcom2, "tableTax_comCompile.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded
### Add core columns to taxa table ###

write.table(data.table(psmelt(phyT.comTot)),"tableMelt_commonTot.csv",row.names=FALSE,sep=",")
write.table(data.table(psmelt(phyT.start)),"tableMelt_Tot.csv",row.names=FALSE,sep=",")

###############################################
################ Practice area ################ 
############################################### 
mat.hv <- sort(taxa_sums(phyT.hv), TRUE)[1:35] # Identify 27 most abundant taxa
phyTmat.hv <- prune_taxa(names(mat.hv),phyT.hv) # Create a subset of data including only 27 most abundant taxa
pBar.matHv <- plot_bar(phyTmat.hv,x="Replicate",fill="Genus2") + 
  facet_grid(Sample_Type~Da1Fb, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.matHv


mat.mdILF0 <- sort(taxa_sums(phyT.mdILF0), TRUE)[1:35] # Identify 27 most abundant taxa
phyTmat.mdILF0 <- prune_taxa(names(mat.mdILF0),phyT.mdILF0) # Create a subset of data including only 27 most abundant taxa
phyTmat.mdILF0 <- subset_taxa(phyT.mdILF0,Order=="Aeromonadales")
sample_data(phyTmat.mdILF0)$WildMonth = factor(sample_data(phyTmat.mdILF0)$WildMonth, levels = c("April","May","June","July","August","September","October")) # Reorder Da1Fb


pBar.month <- plot_bar(phyTmat.mdILF0,x="Sample",fill="Order") + 
  facet_grid(.~WildMonth, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.month

phyTmat.mdILF0<-subset_samples(phyTmat.mdILF0,AnimalSource!="MtSnowVT")
sample_data(phyTmat.mdILF0)$tYear <-with(sample_data(phyTmat.mdILF0),
                                         ifelse(WildMonth=="April","April",
                                         ifelse(WildMonth=="October","October",
                                         "Summer")))

phyT.mer<-tax_glom(phyTmat.mdILF0,"Order")
dt.stripMelt <- psmelt(phyT.mer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)
lowA<-1e-1 # set ymin
pStrip.tYear<-ggplot(dt.stripMelt, aes(x=tYear, y=Abundance, color=AnimalSource)) + 
  ggtitle("") + 
  scale_y_log10(limits=c(lowA,1)) + 
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", color="grey50", na.rm=TRUE) +
  stat_summary(fun.data="plot.median", geom="crossbar", width=0.5, color="black", na.rm=TRUE) +
  geom_jitter(position=position_jitter(0.3)) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_color_brewer(palette="Set1")
pStrip.tYear

phyT.adonis<-phyT.mer # assign data for permanova testing (phyloseq transformed . adonis function)
bray.md<- phyloseq::distance(phyT.adonis, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.adonis)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.md<-adonis(bray.md ~ tYear, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md





phyR.hvBlad<-subset_samples(phyR.sin,Sample_Type=="Bladder" & Taxonomic_ID=="Hverbana") # (Hirudo verbana bladder)
phyT.hvBlad<-transform_sample_counts(phyR.hvBlad, function(x) x/sum(x)) 

mat.hvBlad <- sort(taxa_sums(phyT.hvBlad), TRUE)[1:35] # Identify 27 most abundant taxa
phyTmat.hvBlad <- prune_taxa(names(mat.hvBlad),phyT.hvBlad) # Create a subset of data including only 27 most abundant taxa
pBar.matHvBlad <- plot_bar(phyTmat.hvBlad,x="Sample",fill="Genus2") + 
  facet_grid(Sample_Type~., scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.matHvBlad

phyT.adult2<-subset_samples(phyT.adult,Taxonomic_ID=="Hverbana" & Sample_Type=="ILF" & Da2F%in%c("none","0"))
mat.adult <- sort(taxa_sums(phyT.adult2), TRUE)[1:35] # Identify 27 most abundant taxa
phyTmat.adult <- prune_taxa(names(mat.adult),phyT.adult2) # Create a subset of data including only 27 most abundant taxa
pBar.matAdult <- plot_bar(phyTmat.adult,x="Sample",fill="Genus2") + 
  facet_grid(AnimalSource~Da1F, scales="free_x",space="free") + 
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pal.pairBiome)
pBar.matAdult

