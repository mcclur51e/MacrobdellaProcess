hiAdult<-prune_taxa(taxa_sums(prAdult)>.01,prAdult) # keep taxa with at least 1% of 1 sample

ILFadult<-subset_samples(hiAdult,Sample_Type=="ILF")
ILFhv<-subset_samples(ILFadult,Taxonomic_ID=="Hverbana")
ILF1F<-subset_samples(ILFhv,Da2F%in%c("none","NA","0"))
ILFhAdult<-prune_taxa(taxa_sums(ILF1F)>.01,ILF1F) # keep taxa with at least 1% of 1 sample

matAdult <- sort(taxa_sums(ILFhAdult), TRUE)[1:15] # Identify 15 most abundant taxa
pruAdult <- prune_taxa(names(matAdult),ILFhAdult) # Create a subset of data including only 10 most abundant taxa

plot_bar(pruAdult,fill="Genus")

ILFadult<-subset_samples(hiAdult,Sample_Type=="ILF")
hiro<-subset_samples(hiAdult,Sample_ID%in%c("A042117JGa.b","A042117JGd.b","A042317EMg.b","A042117JGa.i","A042317EMg.i","A050217EMr.i","A042317EMf.u","A042317EMh.u","A042317EMj.u"))




macroS<-subset_samples(macroN,!sample_names(macroN)%in%c("MAa0dF110814EMb.iD691","MAa0dF110814EMc.iD55","MAa0dF110814EMb.uD551","Ma5wF082117EMc.uD695","MAa4dF110514EMa.uD571","MAa4dF110514EMc.uD693","MAa7dF110814EMa.uD555","MAa7dF110814EMc.uD553","MAa4dF110514EMa.iD24","MAa4dF110514EMa.iD18","MAa4dF110514EMb.iD25","MAa4dF110514EMb.uD536","MAa4dF110514EMb.iD24","MAa4dF110514EMb.iD488","MAa0dF110814EMb.iD34","MAa0dF110814EMc.iD450","MAa2dF110314EMc.iD54","Ma061817EMc.iD577","MAa4dF110514EMc.iD91","MAa2dF110314EMc.uD581")) # Remove duplicate samples
ILFadult<-subset_samples(macroS,Sample_Type=="Bladder")
ILF1F<-subset_samples(ILFadult,Da2F%in%c("none","NA","0"))
ILFhAdult<-prune_taxa(taxa_sums(ILF1F)>.01,ILF1F) # keep taxa with at least 1% of 1 sample

matAdult <- sort(taxa_sums(ILFhAdult), TRUE)[1:15] # Identify 15 most abundant taxa
pruAdult <- prune_taxa(names(matAdult),ILFhAdult) # Create a subset of data including only 10 most abundant taxa

plot_bar(ILFhAdult,fill="Genus")