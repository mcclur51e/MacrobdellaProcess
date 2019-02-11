########################################
########## Find specific OTUs ##########
########################################
phyT.aaa<-subset_taxa(phyT.ilf04,Number%in%c("denovo187155"))
median(sample_sums(phyT.aaa))
max(sample_sums(phyT.aaa))
min(sample_sums(phyT.aaa))
tax_table(prune_taxa(names(sort(taxa_sums(phyT.aaa), TRUE)[1:50]),phyT.aaa))[,"Genus2"]
sort(sample_sums(phyT.aaa))
phyT.aab<-subset_samples(phyT.aaa,sample_sums(phyT.aaa)>=.001) # keep samples with >=10000 reads (physeq.minimum 10,000)
median(sample_sums(phyT.aab))
nsamples(phyT.aab)
nsamples(phyT.aaa)
nsamples(phyT.aab) / nsamples(phyT.aaa)

######################################################
########## List prevalence of specific OTUs ##########
######################################################
phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="Bladder"),Taxonomic_ID=="Hverbana")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceHvBlad = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                           by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
hvBlad<-prevdt.aaa[PrevalenceHvBlad > 0,list(TaxaID, PrevalenceHvBlad)]

phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="Intestinum"),Taxonomic_ID=="Hverbana")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceHvInt = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                    by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
hvInt<-prevdt.aaa[PrevalenceHvInt > 0,list(TaxaID, PrevalenceHvInt)]

phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="ILF"),Taxonomic_ID=="Hverbana")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceHvILF = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                    by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
hvILF<-prevdt.aaa[PrevalenceHvILF > 0,list(TaxaID, PrevalenceHvILF)]

phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="Bladder"),Taxonomic_ID=="Mdecora")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceMdBlad = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                    by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
mdBlad<-prevdt.aaa[PrevalenceMdBlad > 0,list(TaxaID, PrevalenceMdBlad)]

phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="Intestinum"),Taxonomic_ID=="Mdecora")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceMdInt = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                    by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
mdInt<-prevdt.aaa[PrevalenceMdInt > 0,list(TaxaID, PrevalenceMdInt)]

phyT.aaa<-subset_samples(subset_samples(phyT.base,Sample_Type=="ILF"),Taxonomic_ID=="Mdecora")
fm.aaa = fast_melt(phyT.aaa) # (fast melt Hirudo verbana bladder)
prevdt.aaa = fm.aaa[, list(PrevalenceMdILF = sum(count >= .001), 
                           TotalPer = sum(count),
                           MinCount = min(count),
                           MaxCount = max(count)),
                    by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana bladder)
mdILF<-prevdt.aaa[PrevalenceMdILF > 0,list(TaxaID, PrevalenceMdILF)]

ls.aaa <- list(hvBlad,hvInt,hvILF,mdBlad,mdInt,mdILF) # list of data.tables to merge (list of taxonomy data.tables)
lapply(ls.aaa, function(i) setkey(i, TaxaID)) # set key for merge function
dt.aaa <- Reduce(function(...) merge(..., all = T), ls.aaa) # merge list of data.tables, keeping values even if not present in each table (Taxonomy table core 2)
write.table(dt.aaa,"table_prevalence.csv",sep=",") # Make list of core OTUs for Macrobdella
