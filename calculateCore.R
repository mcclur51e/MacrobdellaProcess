prAdult<-subset_samples(physeqAn,Age=="A") # keep adult samples (prune Adult)
hiAdult<-prune_taxa(taxa_sums(prAdult)>.01,prAdult) # keep taxa with at least 1% of 1 sample
hiro<-subset_samples(hiAdult,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A050217EMr.i","A042317EMf.u","A042317EMh.u","A042317EMj.u"))
macro<-subset_samples(hiAdult,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
macroS<-subset_samples(macro,!sample_names(macro)%in%c("MAa0dF110814EMb.iD691","MAa0dF110814EMc.iD55","MAa0dF110814EMb.uD551","Ma5wF082117EMc.uD695","MAa4dF110514EMa.uD571","MAa4dF110514EMc.uD693","MAa7dF110814EMa.uD555","MAa7dF110814EMc.uD553","MAa4dF110514EMa.iD24","MAa4dF110514EMa.iD18","MAa4dF110514EMb.iD25","MAa4dF110514EMb.uD536","MAa4dF110514EMb.iD24","MAa0dF110814EMb.iD34","MAa0dF110814EMc.iD450","MAa2dF110314EMc.iD54")) # Remove duplicate samples
ctMacro<-subset_samples(macroS,AnimalSource=="Wlot") # subset containing only CT samples (CT Macrobdella)
maMacro<-subset_samples(macroS,AnimalSource=="GrotonMA") # subset containing only MA samples (MA Macrobdella)
nyMacro<-subset_samples(macroS,AnimalSource=="CarogaNY") # subset containing only NY samples (NY Macrobdella)
vtMacro<-subset_samples(macroS,AnimalSource=="MtSnowVT") # subset containing only VT samples (VT Macrobdella)

cMin<-0.01 # define minimum to count prevalence in a sample
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

cp<-c(.6) # define the level for determining core OTUs, reported as a percent of the samples tested (core percent)
### Calculated with 75% cut-off rather than 95% cut-off due to low sample size ###

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

coreTot = unique(c(coreHvILF,coreHvInt,coreHvBlad,coreCtILF,coreCtInt,coreMaILF,coreMaInt,coreMaBlad,coreNyILF,coreNyBlad))
coreMac = unique(c(coreMacILF,coreMacInt,coreMacBlad))

###Core plot
phyCon<-merge_phyloseq(macroS,hiro) # merge Macrobdella and Hirudo phyloseq objects (physeq control)
phyBIU<-subset_samples(phyCon,Sample_Type!="Ovary") # remove ovary samples (physeq bladder, ILF, intestinum) 
phyBase<-subset_samples(phyBIU,!Da1F%in%c("1","2","4","7")) # remove days after feeding data to leave 0DaF and 5w (physeq base)
phyBasC<-subset_samples(phyBase,!Sample_ID%in%c("Ma5wF082117EMa.i","Ma5wF082117EMa.u","Ma5wF082117EMb.i","Ma5wF082117EMb.u","Ma5wF082117EMc.i","Ma5wF082117EMc.u")) # Remove ILF and intestinum samples from 5wF (physeq base clean)
macro<-subset_samples(phyBasC,Da1F=="0")

### Add new column to map data. 'Header' column will be used for pCore figure
MAPcore<-sample_data(phyBasC) # pull map data from phyloseq object
MAPcore$Header<-with(MAPcore,
                     ifelse(Taxonomic_ID=="Hverbana","Germany",
                            ifelse(AnimalSource=="GrotonMA","MA",
                                   ifelse(AnimalSource=="CarogaNY","NY",
                                          ifelse(AnimalSource=="Wlot","CT",
                                                 ifelse(AnimalSource=="MtSnowVT","VT",
                                                        as.character(AnimalSource)))))))
corePhy = merge_phyloseq(phyBasC,MAPcore) # return map data to the phyloseq object

coreTab<-prune_taxa(coreTot,corePhy)

########## Macrobdella Sites plot with 'Other' category ##########
dtCore<-data.table(psmelt(corePhy))
dtCore$Number<-as.character(dtCore$Number)
dtCore$Genus<-as.character(dtCore$Genus)
dtCore[!dtCore$Number%in%as.character(coreTot),]$Genus<-" "
dtCore$Color<-with(dtCore,
                   ifelse(Genus=="Aeromonas","#39ac39",
                          ifelse(Genus=="Bacteroides","#990000",
                                 ifelse(Genus=="Mucinivorans","#e60000",
                                        ifelse(Genus=="AF12","#ffb3b3",
                                               ifelse(Genus=="Bdellovribio","#000080",
                                                      ifelse(Genus=="unk_Comamonadaceae","#663300",
                                                             ifelse(Genus=="Proteocatella","#600080",
                                                                    ifelse(Genus=="Butyricicoccus","#ac00e6",
                                                                           ifelse(Genus=="unk_Ruminococcaceae","#f2ccff",
                                                                                  ifelse(Genus=="unk_Desulfovibrionaceae","#ffff00",
                                                                                         ifelse(Genus=="Ochrobactrum","#006680",
                                                                                                ifelse(Genus=="unk_Rhizobiales","#00ccff",
                                                                                                       ifelse(Genus=="Azospirillum","#808080",
                                                                                                              ifelse(Genus=="unk_Rhodospirillaceae","#d9d9d9",
                                                                                                                     ifelse(Genus=="Pedobacter","#ccf5ff",
                                                                                                                            ifelse(Genus=="unk_Alphaproteobacteria","#ff00aa",
                                                                                                                                   ifelse(Genus=="unk_Betaproteobacteria","#ffb3e6",
                                                                                                                                          "#1a1a1a"))))))))))))))))))

pCore <- ggplot(dtCore, aes(x=Replicate, y=Abundance, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Taxonomic_ID+Header, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pCore # print plot

##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pCore), size = "last")), filename="plotCore.png", width=8,height=8)
##### Figure #####