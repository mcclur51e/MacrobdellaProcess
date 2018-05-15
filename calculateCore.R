# Thank you to jeffkimbrel for the next two lines! #
taxaDC <- setdiff(taxa_names(physeqAn), taxa_names(ps.allc)) # find the difference between taxa in contaminant list and physeqAn (taxa decontam)
physeqDC <- prune_taxa(taxaDC, physeqAn) # keep only taxa not in taxaNC (physeq decontam)

prAdult<-subset_samples(physeqDC,Age=="A") # keep adult samples (prune Adult)
hiAdult<-prune_taxa(taxa_sums(prAdult)>.00001,prAdult) # keep taxa with at least .001% of 1 sample (number chosen to reduce # taxa as low as possible while bar chart still appears to add to 1. 
hiro<-subset_samples(hiAdult,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A050217EMr.i","A042317EMf.u","A042317EMh.u","A042317EMj.u"))
macro<-subset_samples(hiAdult,Taxonomic_ID%in%c("Mdecora","Unk")) # subset containing only Macrobdella samples
macroN<-subset_samples(macro,!Replicate%in%c("MN2","MN3","MN4"))
macroS<-subset_samples(macroN,!sample_names(macroN)%in%c("MAa0dF110814EMb.iD691","MAa0dF110814EMc.iD55","MAa0dF110814EMb.uD551","Ma5wF082117EMc.uD695","MAa4dF110514EMa.uD571","MAa4dF110514EMc.uD693","MAa7dF110814EMa.uD555","MAa7dF110814EMc.uD553","MAa4dF110514EMa.iD24","MAa4dF110514EMa.iD18","MAa4dF110514EMb.iD25","MAa4dF110514EMb.uD536","MAa4dF110514EMb.iD24","MAa4dF110514EMb.iD488","MAa0dF110814EMb.iD34","MAa0dF110814EMc.iD450","MAa2dF110314EMc.iD54","Ma061817EMc.iD577","MAa4dF110514EMc.iD91","MAa2dF110314EMc.uD581")) # Remove duplicate samples
ctMacro<-subset_samples(macroS,AnimalSource=="Wlot") # subset containing only CT samples (CT Macrobdella)
maMacro<-subset_samples(macroS,AnimalSource=="GrotonMA") # subset containing only MA samples (MA Macrobdella)
nyMacro<-subset_samples(macroS,AnimalSource=="CarogaNY") # subset containing only NY samples (NY Macrobdella)
vtMacro<-subset_samples(macroS,AnimalSource=="MtSnowVT") # subset containing only VT samples (VT Macrobdella)

cMin<-0.01 # define minimum to count prevalence in a sample
cp<-c(.5) # define the level for determining core OTUs, reported as a percent of the samples tested (core percent)
### Calculated with 75% cut-off rather than 95% cut-off due to low sample size ###

###Hirudo
hvBlad<-subset_samples(hiro,Sample_Type=="Bladder") # (Hirudo verbana bladder)
fmhvBlad = fast_melt(hvBlad) # (fast melt Hirudo verbana bladder)
prevdtHvBlad = fmhvBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
hvInt<-subset_samples(hiro,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fmhvInt = fast_melt(hvInt) # (fast melt Hirudo verbana intestinum)
prevdtHvInt = fmhvInt[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
hvILF<-subset_samples(hiro,Sample_Type=="ILF") # (Hirudo verbana ILF)
fmhvILF = fast_melt(hvILF) # (fast melt Hirudo verbana ILF)
prevdtHvILF = fmhvILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###Macrobdella		
macBlad<-subset_samples(macroS,Sample_Type=="Bladder") # (Hirudo verbana bladder)
fmmacBlad = fast_melt(macBlad) # (fast melt Hirudo verbana bladder)
prevdtmacBlad = fmmacBlad[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
macInt<-subset_samples(macroS,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fmmacInt = fast_melt(macInt) # (fast melt Hirudo verbana intestinum)
prevdtmacInt = fmmacInt[, list(Prevalence = sum(count >= cMin),
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
macILF<-subset_samples(macroS,Sample_Type=="ILF") # (Hirudo verbana ILF)
fmmacILF = fast_melt(macILF) # (fast melt Hirudo verbana ILF)
prevdtmacILF = fmmacILF[, list(Prevalence = sum(count >= cMin), 
  TotalPer = sum(count),
  MinCount = min(count),
  MaxCount = max(count)),
  by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###CT
ctBlad<-subset_samples(ctMacro,Sample_Type=="Bladder") # (Hirudo verbana bladder)
fmctBlad = fast_melt(ctBlad) # (fast melt Hirudo verbana bladder)
prevdtctBlad = fmctBlad[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
ctInt<-subset_samples(ctMacro,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fmctInt = fast_melt(ctInt) # (fast melt Hirudo verbana intestinum)
prevdtctInt = fmctInt[, list(Prevalence = sum(count >= cMin),
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana intestinum)
ctILF<-subset_samples(ctMacro,Sample_Type=="ILF") # (Hirudo verbana ILF)
fmctILF = fast_melt(ctILF) # (fast melt Hirudo verbana ILF)
prevdtctILF = fmctILF[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
###MA
maBlad<-subset_samples(maMacro,Sample_Type=="Bladder") # (Massachusetts bladder)
fmmaBlad = fast_melt(maBlad) # (fast melt Massachusetts bladder)
prevdtmaBlad = fmmaBlad[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts bladder)
maInt<-subset_samples(maMacro,Sample_Type=="Intestinum") # (Massachusetts intestinum)
fmmaInt = fast_melt(maInt) # (fast melt Massachusetts intestinum)
prevdtmaInt = fmmaInt[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts intestinum)
maILF<-subset_samples(maMacro,Sample_Type=="ILF") # (Massachusetts ILF)
fmmaILF = fast_melt(maILF) # (fast melt Massachusetts ILF)
prevdtmaILF = fmmaILF[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Massachusetts ILF)
###NY
nyBlad<-subset_samples(nyMacro,Sample_Type=="Bladder") # (New York bladder)
fmnyBlad = fast_melt(nyBlad) # (fast melt New York bladder)
prevdtnyBlad = fmnyBlad[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # nyke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York bladder)
nyInt<-subset_samples(nyMacro,Sample_Type=="Intestinum") # (New York intestinum)
fmnyInt = fast_melt(nyInt) # (fast melt New York intestinum)
prevdtnyInt = fmnyInt[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # nyke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York intestinum)
nyILF<-subset_samples(nyMacro,Sample_Type=="ILF") # (New York ILF)
fmnyILF = fast_melt(nyILF) # (fast melt New York ILF)
prevdtnyILF = fmnyILF[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # nyke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table New York ILF)
###VT
vtBlad<-subset_samples(vtMacro,Sample_Type=="Bladder") # (Vermont bladder)
fmvtBlad = fast_melt(vtBlad) # (fast melt Vermont bladder)
prevdtvtBlad = fmvtBlad[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # vtke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont bladder)
vtInt<-subset_samples(vtMacro,Sample_Type=="Intestinum") # (Vermont intestinum)
fmvtInt = fast_melt(vtInt) # (fast melt Vermont intestinum)
prevdtvtInt = fmvtInt[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # vtke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont intestinum)
vtILF<-subset_samples(vtMacro,Sample_Type=="ILF") # (Vermont ILF)
fmvtILF = fast_melt(vtILF) # (fast melt Vermont ILF)
prevdtvtILF = fmvtILF[, list(Prevalence = sum(count >= cMin), 
                             TotalPer = sum(count),
                             MinCount = min(count),
                             MaxCount = max(count)),
                      by = TaxaID] # vtke simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont ILF)

coreHvILF = prevdtHvILF[(Prevalence >= cp*nsamples(subset_samples(hvILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo ILF
coreHvInt = prevdtHvInt[(Prevalence >= cp*nsamples(subset_samples(hvInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo Intestinum
coreHvBlad = prevdtHvBlad[(Prevalence >= cp*nsamples(subset_samples(hvBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo Bladder

coreCtILF = prevdtctILF[(Prevalence >= cp*nsamples(subset_samples(ctILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo ILF 
coreCtInt = prevdtctInt[(Prevalence >= cp*nsamples(subset_samples(ctInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo Intestinum
coreCtBlad = prevdtctBlad[(Prevalence >= cp*nsamples(subset_samples(ctBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Hirudo Bladder

coreMaILF = prevdtmaILF[(Prevalence >= cp*nsamples(subset_samples(maILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts ILF (
coreMaInt = prevdtmaInt[(Prevalence >= cp*nsamples(subset_samples(maInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts Intestinum
coreMaBlad = prevdtmaBlad[(Prevalence >= cp*nsamples(subset_samples(maBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Massachusetts Bladder

coreNyILF = prevdtnyILF[(Prevalence >= cp*nsamples(subset_samples(nyILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella New York ILF
coreNyInt = prevdtnyInt[(Prevalence >= cp*nsamples(subset_samples(nyInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella New York Intestinum
coreNyBlad = prevdtnyBlad[(Prevalence >= cp*nsamples(subset_samples(nyBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella New York Bladder

coreVtILF = prevdtvtILF[(Prevalence >= cp*nsamples(subset_samples(vtILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont ILF (
coreVtInt = prevdtvtInt[(Prevalence >= cp*nsamples(subset_samples(vtInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont Intestinum
coreVtBlad = prevdtvtBlad[(Prevalence >= cp*nsamples(subset_samples(vtBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont Bladder

coreMacILF = prevdtmacILF[(Prevalence >= cp*nsamples(subset_samples(macILF,Sample_Type=="ILF")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont ILF (
coreMacInt = prevdtmacInt[(Prevalence >= cp*nsamples(subset_samples(macInt,Sample_Type=="Intestinum")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont Intestinum
coreMacBlad = prevdtmacBlad[(Prevalence >= cp*nsamples(subset_samples(macBlad,Sample_Type=="Bladder")) & MaxCount >= .01), TaxaID] # Make list of core OTUs for Macrobdella Vermont Bladder

coreTot = unique(c(coreHvILF,coreHvInt,coreHvBlad,coreCtILF,coreCtInt,coreMaILF,coreMaInt,coreMaBlad,coreNyILF,coreNyBlad,coreVtBlad))
coreMac = unique(c(coreMacILF,coreMacInt,coreMacBlad))
coreMdHv = unique(c(coreMacILF,coreMacInt,coreMacBlad,coreHvILF,coreHvInt,coreHvBlad))

###Core plot
phyCon<-merge_phyloseq(macroS,hiro) # merge Macrobdella and Hirudo phyloseq objects (physeq control)
phyBIU<-subset_samples(phyCon,Sample_Type!="Ovary") # remove ovary samples (physeq bladder, ILF, intestinum) 
phyBase<-subset_samples(phyBIU,Da1F%in%c("0","35","99","101","110")) # remove days after feeding data to leave 0DaF and 5w (physeq base)
phyBasC<-subset_samples(phyBase,!Sample_ID%in%c("Ma5wF082117EMa.i","Ma5wF082117EMa.u","Ma5wF082117EMb.i","Ma5wF082117EMb.u","Ma5wF082117EMc.i","Ma5wF082117EMc.u")) # Remove ILF and intestinum samples from 5wF (physeq base clean)

### Add new column to map data. 'Header' column will be used for pCore figure
MAPcore<-sample_data(phyBasC) # pull map data from phyloseq object
MAPcore$Header<-with(MAPcore,
  ifelse(Taxonomic_ID=="Hverbana","Germany",
  ifelse(AnimalSource=="GrotonMA","MA",
  ifelse(AnimalSource=="CarogaNY","NY",
  ifelse(AnimalSource=="Wlot","CT",
  ifelse(AnimalSource=="MtSnowVT","VT",
  as.character(AnimalSource)))))))
corePhy = merge_phyloseq(phyBasC,MAPcore) # return map data to the phyloseq object (Core phyloseq)

coreTab<-prune_taxa(coreTot,corePhy) # keep only taxa from the identified 'Core' 

########## Macrobdella Sites plot with 'Other' category ##########
#matCore <- sort(taxa_sums(corePhy), TRUE)[1:35] # Identify 27 most abundant taxa
#pruCore <- prune_taxa(names(matCore),corePhy) # Create a subset of data including only 27 most abundant taxa

dtCore<-data.table(psmelt(corePhy))
dtCore$Number<-as.character(dtCore$Number)
dtCore$Genus<-as.character(dtCore$Genus)
dtCore[!dtCore$Number%in%as.character(coreTot),]$Genus<-""

dtCore$Genus<-with(dtCore,
  ifelse(Genus=="","other",
  as.character(Genus)))

dtCore$Genus <- factor(dtCore$Genus, levels=c("Aeromonas",
  "Bacteroides","PW3","AF12",
  "unk_Ruminococcaceae","unk_Clostridiales",
  "Ochrobactrum","Pedobacter","unk_Comamonadaceae","Bdellovibrio",
  "unk_Rhodospirillaceae","unk_Rhodospirillales",
  "unk_Myxococcales","unk_Rhodocyclaceae",
  "unk_Rhizobiales","unk_Rhizobiaceae",
  "Agrobacterium","unk_Desulfovibrionaceae",
  "unk_Alpha","unk_Beta","other"))

genus.color<-c(Aeromonas="#39ac39",
  Bacteroides="#B80000",PW3="#F00000",AF12="#FF7777",
  unk_Ruminococcaceae="#000080",unk_Clostridiales="#3333ff",
  Ochrobactrum="#ff6600",Pedobacter="#ff9933",unk_Comamonadaceae="#ffbf80",Bdellovibrio="#FFCC99",
  unk_Rhodospirillaceae="#600080",unk_Rhodospirillales="#ac00e6",
  unk_Myxococcales="#ffff00",unk_Rhodocyclaceae="#ffffb3",
  unk_Rhizobiales="#006680",unk_Rhizobiaceae="#00ccff",
  Agrobacterium="#ff00aa",unk_Desulfovibrionaceae="#ffb3e6",
  unk_Alpha="#d9d9d9",unk_Beta="#ffffff", other="#808080")


pCore <- ggplot(dtCore, aes(x=Replicate, y=Abundance, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pCore # print plot

##### Figure + Table #####
ggsave(grid.draw(rbind(ggplotGrob(pCore), size = "last")), filename="plotCore.png", width=12,height=8)
write.table(tax_table(pruCore), "taxTable.csv", sep=",")
##### Figure + Table #####

### Clean-up using .999 with NMDS plot (removes 3 ILF and 1 intestinum samples) ###
clean99<-subset_samples(corePhy,!sample_names(coreTab)%in%c("MAa0dF110814EMc.iD506","Ma061817EMc.iD486","MAa0dF110814EMa.iD493","W66B.N","MAa0dF110814EMb.iD449"))

dtCore<-data.table(psmelt(clean99))
dtCore$Number<-as.character(dtCore$Number)
dtCore$Genus<-as.character(dtCore$Genus)
dtCore[!dtCore$Number%in%as.character(coreTot),]$Genus<-""

dtCore$Genus<-with(dtCore,
  ifelse(Genus=="","other",
  as.character(Genus)))

dtCore$Genus <- factor(dtCore$Genus, 
  levels=c("Aeromonas",
           "Bacteroides","PW3","AF12",
           "unk_Ruminococcaceae","unk_Clostridiales",
           "Ochrobactrum","Pedobacter","unk_Comamonadaceae","Bdellovibrio",
           "unk_Rhodospirillaceae","unk_Rhodospirillales",
           "unk_Myxococcales","unk_Rhodocyclaceae",
           "unk_Rhizobiales","unk_Rhizobiaceae",
           "Agrobacterium","unk_Desulfovibrionaceae",
           "unk_Alpha","unk_Beta","other"))

genus.color<-c(Aeromonas="#39ac39",
               Bacteroides="#B80000",PW3="#F00000",AF12="#FF7777",
               unk_Ruminococcaceae="#000080",unk_Clostridiales="#3333ff",
               Ochrobactrum="#ff6600",Pedobacter="#ff9933",unk_Comamonadaceae="#ffbf80",Bdellovibrio="#FFCC99",
               unk_Rhodospirillaceae="#600080",unk_Rhodospirillales="#ac00e6",
               unk_Myxococcales="#ffff00",unk_Rhodocyclaceae="#ffffb3",
               unk_Rhizobiales="#006680",unk_Rhizobiaceae="#00ccff",
               Agrobacterium="#ff00aa",unk_Desulfovibrionaceae="#ffb3e6",
               unk_Alpha="#d9d9d9",unk_Beta="#ffffff", other="#808080")

pCore99 <- ggplot(dtCore, aes(x=Replicate,y=Abundance, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pCore99 # print plot

### Add core columns to taxa table ###
coreMdHv<-prune_taxa(coreMdHv,corePhy) # keep only taxa from the identified 'Core' 
TAXcore<-tax_table(coreMdHv) # pull taxonomy data from phyloseq object
dt.TAXcore<-as.data.table(TAXcore) # Change taxa_table to a data.table
dt.TAXcore$Number<-as.character(dt.TAXcore$Number) # force TAXcore$Number to be a character

# Macrobdella core data
taxdt.MacInt<-as.data.table(coreMacInt)
taxdt.MacInt$Md.Intestinum<-"x"
setnames(taxdt.MacInt,"coreMacInt", "Number") # First column = OTU names. Change column heading

taxdt.MacILF<-as.data.table(coreMacILF)
taxdt.MacILF$Md.ILF<-"x"
setnames(taxdt.MacILF,"coreMacILF", "Number") # First column = OTU names. Change column heading

taxdt.MacBlad<-as.data.table(coreMacBlad)
taxdt.MacBlad$Md.Bladder<-"x"
setnames(taxdt.MacBlad,"coreMacBlad", "Number") # First column = OTU names. Change column heading

#Hirudo core data
taxdt.HvInt<-as.data.table(coreHvInt)
taxdt.HvInt$Hv.Intestinum<-"x"
setnames(taxdt.HvInt,"coreHvInt", "Number") # First column = OTU names. Change column heading

taxdt.HvILF<-as.data.table(coreHvILF)
taxdt.HvILF$Hv.ILF<-"x"
setnames(taxdt.HvILF,"coreHvILF", "Number") # First column = OTU names. Change column heading

taxdt.HvBlad<-as.data.table(coreHvBlad)
taxdt.HvBlad$Hv.Bladder<-"x"
setnames(taxdt.HvBlad,"coreHvBlad", "Number") # First column = OTU names. Change 

# Merge tables
list.taxdt <- list(dt.TAXcore,taxdt.MacILF,taxdt.MacInt,taxdt.MacBlad,taxdt.HvILF,taxdt.HvInt,taxdt.HvBlad) # list of data.tables to merge (list of taxonomy data.tables)
lapply(tbls, function(i) setkey(i, Number)) # set key for merge function
TAXcore2 <- Reduce(function(...) merge(..., all = T), list.taxdt) # merge list of data.tables, keeping values even if not present in each table (Taxonomy table core 2)
TAXcore2[is.na(TAXcore2)] <- "" # replace <NA> values with empty values
write.table(TAXcore2, "taxTableCore.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded
### Add core columns to taxa table ###
