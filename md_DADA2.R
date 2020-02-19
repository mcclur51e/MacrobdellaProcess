#####################################################################
########## Processing Set −up #######################################
#####################################################################
########## Call libraries for use ##########
########## Call libraries for use ##########
library("Cairo") # allow use of Times font for graphics # version 1.5-10
library("phyloseq") # version 1.28.0
library("ggplot2") # version 3.2.1
library("ggpubr") # for stat_compare_means command = diversity calculations of plot_richness # version 0.2.3
library("reshape2") # version 1.4.3
library("data.table") # version 1.12.6
library("RColorBrewer") # design new color palette for data # version 1.1-2
library("scales") # for scientific notation in plots # version 1.0.0
library("cowplot") # used to prepare final figures # version 1.0.0
library("decontam") # identify contaminant ASVs # version 1.4.0
library("pairwiseAdonis") # for PERMANOVA calculations # version 0.0.1
library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. # version 1.24.0
library("dada2") # version 1.12.1
library("Biostrings") # 2.52.0
library("tidyr") # version 1.0.0
library("dplyr") # version 0.8.3
library("pheatmap") # version 1.0.12
library("vegan") # version 2.5-6
library("microbiome") # version 1.6.0
#library("ape") #packageVersion("ape")
#library("gridExtra")
#library("structSSI")

########## Call functions for use ##########
source("~/Masters/R_tools/taxa_summary",local=TRUE) # load fast_melt function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
plot.median <- function(x) {
  m <- median(x)
  c(y = m, ymin = m, ymax = m)
}
########## Define color palettes ##########
pal.pairBiome<-c("#1c541c","#2f8f2f","#49c349","#bfeabf",
                 "#B80000","#F00000","#FF7777","#ffcccc",
                 "#000080","#0000cd","#8282ff","#cfcfff",
                 "#623800","#c47000","#ff9914","#ffddb1",
                 "#430059","#7d00a7","#cc32ff","#eebbff",
                 "#626200","#c4c400","#ffff14","#ffffb1",
                 "#005f6c","#00a4bb","#1ee3ff","#bbf7ff",
                 "#750063","#c400a5","#ff13da","#ffb0f3",
                 "#1a1a1a","#808080","#d9d9d9","#ffffff","#808080")
pal.pairMini<-c("#A6CEE3","#1F78B4",
                "#B2DF8A","#33A02C",
                "#FB9A99","#E31A1C",
                "#FDBF6F","#FF7F00",
                "#CAB2D6","#6A3D9A",
                "#FFFF99","#eded09",
                "#9aebff","#01cdff",
                "#e6e6e6","#c0c0c0")
pal.CB<-c("#e69f00","#56b4e9","#009e73","#f0e442","#0072b2","#D55e00","#cc79a7","#999999","#000000") # colorblind-friendly color palette

setwd("~/Desktop/md_DADA2/") # set working directory
#####################################################################
##### ASV assignment with DADA2 #####################################
#####################################################################
path <- "~/Desktop/md_DADA2/" # directory containing fastq files

fnFs <- sort(list.files(path, pattern="R1_001", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2_001", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
plotQualityProfile(fnFs[1:2]) # plot quality profile of forward sequences
plotQualityProfile(fnRs[1:2]) # plot quality profile of reverse sequences

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # 
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 240:263]
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/md_DADA2/silva_nr_v132_train_set.fa", multithread=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

taxa2<-tibble::rownames_to_column(as.data.frame(taxa),"sequence")
taxa3<-taxa
setDT(as.data.frame(taxa3), keep.rownames = "sequence")[]

#####################################################################
########## Transfer data into Phyloseq ##############################
#####################################################################
OTU = otu_table(seqtab.nochim, taxa_are_rows=FALSE) # assigns ASV table from DADA output
TAX = tax_table(taxa) # reads data.frame into matrix without converting to ASVs # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with taxa data separated by columns and headers added for Kingdom, Phylum, Class, etc. First column = ASV IDs
MAP = sample_data(data.frame(read.csv("DataFiles/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads csv file into data.frame with row names in column 1
physeq = phyloseq(OTU, TAX, MAP)
### create and assign ASV numbers to consensus sequences for easier reference ###
ASV <- paste0("ASV", seq(ntaxa(physeq))) 
Sequence <- row.names(tax_table(physeq)) 
bind.asv <- cbind(tax_table(physeq),ASV)
bind.seq <- cbind(bind.asv,Sequence) 
TAX2 = tax_table(as.matrix(bind.seq)) # define new taxonomy table
physeq = phyloseq(OTU, TAX2, MAP) # redefine physeq with new taxonomy table

### These three lines are for use if data has been processed once, then ASV taxonomy modified by hand
setwd("~/Desktop/md_DADA2/") # set working directory for data processing. Must contain folder "DataFiles" containing "table_ASV.csv", "table_tax.csv", and "table_map.csv" and/or "physeq_manuscript.RData"
#write.csv(tax_table(physeq),"DataFiles/table_tax.csv") # output current tax_table to .csv file and allow for editing by hand
TAXm = tax_table(as.matrix.data.frame(read.csv("DataFiles/table_tax_modify.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads data.frame into matrix without converting to numbers # reads csv file into data.frame with row names in column 1 # use to input modified tax_table
physeq = phyloseq(OTU, TAXm, MAP) # recreate physeq object with modified tax_table

save(physeq,file=("DataFiles/physeq_current.RData")) # Save the phyloseq data object in a .RData file 
load("DataFiles/physeq_current.RData") # Use to re-load saved phyloseq data object 
#####################################################################
########## Pre-process data to identify contaminant ASVs  ###########
#####################################################################
phyR.PCRneg <- subset_samples(physeq,Sample_Type%in%c("PCRneg","Reagents")) # keep PCR negative control samples
phyR.taxaNeg <- subset_taxa(phyR.PCRneg,taxa_sums(phyR.PCRneg)>1) # keep only ASVs with >1 read total
var.maxNeg <- max(otu_table(phyR.taxaNeg)) # find the greatest count of any ASV in a single negative control
is.na(otu_table(phyR.taxaNeg)) <- otu_table(phyR.taxaNeg)==0 # change 0 values to NA
var.medNeg <- median(as.vector(otu_table(phyR.taxaNeg)),na.rm=TRUE) # calculate median on values that ar not NA (0)

### Trim unuseable samples ###
phyR.negCount <- subset_samples(physeq,!PostPCRDNA_ng_uL%in%c("Unk") & sample_sums(physeq)>=1) # Remove any samples where i) Post PCR [DNA] was not calculated and that ii) have >=1 reads
phyR.1min <- prune_taxa(taxa_sums(phyR.negCount)>1,phyR.negCount) # keep only taxa with >1 reads

### Map ###
sample_data(phyR.1min)$Conc <- with(sample_data(phyR.1min),
                                    ifelse(PostPCRDNA_ng_uL%in%c("zero","Unk",""), "0",
                                    as.character(PostPCRDNA_ng_uL)))
sample_data(phyR.1min)$c <- as.numeric(as.character(sample_data(phyR.1min)$Conc)) * 10 + 1 # math to make decontam work
sample_data(phyR.1min)$q2 <- as.numeric(as.character(sample_data(phyR.1min)$c))

### Identify contaminants by frequency ###
contamdf.freq <- isContaminant(phyR.1min, method="frequency", conc="q2")
phyR.contamFreq <- prune_taxa(contamdf.freq$contaminant, phyR.1min)
phyR.allc <- prune_taxa(taxa_names(phyR.contamFreq), phyR.1min) # (physeq all contaminants)
phyR.nc <- subset_taxa(physeq,!taxa_names(physeq)%in%c(taxa_names(phyR.contamFreq)))
#####################################################################
########## Histogram of reads counts  ###############################
#####################################################################
sdt = data.table(as(sample_data(phyR.nc), "data.frame"),
                 TotalReads = sample_sums(physeq), 
                 keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth + facet_wrap(~Sample_Type) + scale_x_log10()
#####################################################################
########## Rarefy and describe data #################################
#####################################################################
rarecurve(t(otu_table(phyR.nc)), step = 20, xlab = "Reads", ylab = "ASVs", label = FALSE, xlim=c(0,10000))
phyR.describe <- subset_samples(phyR.nc, !Sample_Type%in%c("PCRneg","Reagents"))
min(sample_sums(phyR.describe)) # minimum reads in a sample -> use for determining rarefaction depth
max(sample_sums(phyR.describe)) # maximum reads in a sample
median(sample_sums(phyR.describe)) # median reads per sample

#####################################################################
########## Set-up directories for processing data  ##################
#####################################################################
dir.create("Process",showWarnings = FALSE) # create folder for process data
setwd("~/Desktop/md_DADA2/Process/") # set working directory for data processing. Must contain folder "DataFiles" containing "table_ASV.csv", "table_tax.csv", and "table_map.csv" and/or "physeq_manuscript.RData"
dir.create("NMDS",showWarnings = FALSE) # create folder for NMDS data output(s)
dir.create("Plots",showWarnings = FALSE) # create folder for plot output(s)

#####################################################################
########## Add columns to map data for use in plotting ##############
#####################################################################
phyR.sin <- subset_taxa(phyR.describe,taxa_sums(phyR.describe)>1) # keep only ASVs with >1 read total
### 'Header' column will be used for many figures ###
sample_data(phyR.sin)$Header <- with(sample_data(phyR.sin),
                                     ifelse(AnimalSource=="USA","",
                                     ifelse(AnimalSource=="BBEZ","",
                                     ifelse(AnimalSource=="GrotonMA","MA",
                                     ifelse(AnimalSource=="CarogaNY","NY",
                                     ifelse(AnimalSource=="Wlot","CT",
                                     ifelse(AnimalSource=="MtSnowVT","VT",
                                     ifelse(AnimalSource=="Schoolhouse_Brook","CT",
                                     as.character(AnimalSource)))))))))
### 'Da1fb' column will be used to merge some time points together. ###
sample_data(phyR.sin)$Da1Fb <- with(sample_data(phyR.sin),
                                    ifelse(Da1F=="8","7",
                                    ifelse(Da1F=="31","30+",
                                    ifelse(Da1F=="35","30+",
                                    ifelse(Da1F=="76","30+",
                                    ifelse(Da1F=="82","30+",
                                    ifelse(Da1F=="90","90+",
                                    ifelse(Da1F=="92","90+",
                                    ifelse(Da1F=="94","90+",
                                    ifelse(Da1F=="96","90+",
                                    ifelse(Da1F=="97","90+",
                                    ifelse(Da1F=="99","90+",
                                    ifelse(Da1F=="101","90+",
                                    ifelse(Da1F=="108","90+",
                                    ifelse(Da1F=="110","90+",
                                    ifelse(Da1F=="112","90+",
                                    ifelse(Da1F=="113","90+",
                                    ifelse(Da1F=="116","90+",
                                    ifelse(Da1F=="118","90+",
                                    ifelse(Da1F=="130","90+",
                                    ifelse(Da1F=="165","90+",
                                    ifelse(Da1F=="166","90+",
                                    ifelse(Da1F=="167","90+",
                                    ifelse(Da1F=="169","90+",
                                    ifelse(Da1F=="172","90+",
                                    ifelse(Da1F=="193","90+",
                                    ifelse(Da1F=="215","90+",
                                    ifelse(Da1F=="265","90+",
                                    as.character(Da1F))))))))))))))))))))))))))))) 
### 'spp' column will be used to identify leech species ###
sample_data(phyR.sin)$spp <- with(sample_data(phyR.sin),
                                  ifelse(Taxonomic_ID=="Hverbana","Hirudo verbana",
                                  ifelse(Taxonomic_ID=="Mdecora","Macrobdella decora",
                                  as.character(Taxonomic_ID))))
### 'Season' column will be used to identify season in which leech was collected ###
sample_data(phyR.sin)$Season <- with(sample_data(phyR.sin),
                                     ifelse(WildMonth=="April","April",
                                     ifelse(WildMonth=="June","Warm",
                                     ifelse(WildMonth=="July","Warm",
                                     ifelse(WildMonth=="August","Warm",
                                     ifelse(WildMonth=="September","Warm",
                                     ifelse(WildMonth=="October","October",
                                     as.character(AnimalSource))))))))
sample_data(phyR.sin)$TaxType<-with(sample_data(phyR.sin),paste(Taxonomic_ID,Sample_Type))
########## Define expressions for plot labels ##########
exp.Hv<-expression(italic("Hirudo verbana"))
exp.HvILF<-expression(paste(italic("H. verbana"), " ILF"))
exp.HvInt<-expression(paste(italic("H. verbana"), " Intestinum"))
exp.Md<-expression(italic("Macrobdella decora"))
exp.MdILF<-expression(paste(italic("M. decora"), " ILF"))
exp.MdInt<-expression(paste(italic("M. decora"), " Intestinum"))

#####################################################################
########## Rarefy Sample Data #######################################
#####################################################################
cRare <- 5000 # value for rarefaction; should be the same in next line (constant Rarefy)
phyR.rare<-rarefy_even_depth(phyR.sin,sample.size=min(5000),rngseed=FALSE,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
cMin <- var.medNeg/cRare # define minimum to count presence in a trasnformed sample (count Minimum)
cpR<-c(.9) # define the level for determining core ASVs, reported as a percent of the samples tested (constant percent coRe)
cpM<-c(.7) # define the level for determining common ASVs, reported as a percent of the samples tested (constant percent coMmon)

#####################################################################
########## Pre-process data to remove low count ASVs, ###############
########## contaminant ASVs, and duplicate samples    ###############
#####################################################################
dt.phyPru = fast_melt(phyR.rare) # make data table from physeq object (data table. physeq Pruned)
prev.phyPru = dt.phyPru[, list(Prevalence = sum(count > 1),
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = ASV] # make simple table listing 'ASV, Prevalence, TotalPer, MinCount, and MaxCount' (prevalence . physeq Pruned)
ls.Pres = prev.phyPru[(Prevalence > 1), ASV] # Make list of ASVs present in dataset at least once (list . Present)
phyR.pru2 <- subset_taxa(phyR.rare,ASV%in%c(ls.Pres)) # remove taxa not in ls.Pres from phyloseq object (physeq Raw . pruned 2)
phyT.leech <- transform_sample_counts(phyR.pru2, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . leech)

#####################################################################
########## Initial sample subsetting ################################
#####################################################################
### Hirudo verbana samples ###
phyT.hv <- subset_samples(phyT.leech, Taxonomic_ID=="Hverbana") # subset containing only Hirudo verbana samples (physeq Transform .  hirudo verbana)
phyT.hvILF <- subset_samples(phyT.leech, Taxonomic_ID=="Hverbana" & Sample_Type=="ILF") # subset containing only Hirudo verbana ILF samples (physeq Transform . Hirudo verbana ILF)
phyT.hvILFbbez <- subset_samples(phyT.hvILF, AnimalSource=="BBEZ") # subset containing only Hirudo verbana ILF samples from animals from supplier 1 (physeq Transform . hirudo verbana ILF bbez)
phyT.hvILFusa <- subset_samples(phyT.hvILF, AnimalSource=="USA") # subset containing only Hirudo verbana ILF samples from animals from supplier 2 (physeq Transform . hirudo verbana ILF usa)
phyT.hvInt <- subset_samples(phyT.leech, Taxonomic_ID=="Hverbana" & Sample_Type=="Intestinum") # subset containing only Hirudo verbana intestinum samples (physeq Transform . Hirudo verbana intestinum)
phyT.hvGut <- subset_samples(phyT.hv, Sample_Type%in%c("ILF","Intestinum")) # subset containing only Hirudo verbana ILF and intestinum samples (physeq Transform . hirudo verbana Gut)
### Macrobdella decora samples ###
phyT.md <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora") # subset containing only Macrobdella decora samples (physeq Transform . macrobdella decora)
phyT.mdCT <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & AnimalSource=="Wlot") # subset containing only Macrobdella decora CT samples (physeq Transform . macrobdella decora CT)
phyT.mdMA <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & AnimalSource=="GrotonMA") # subset containing only Macrobdella decora MA samples (physeq Transform . macrobdella decora MA)
phyT.mdNY <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & AnimalSource=="CarogaNY") # subset containing only Macrobdella decora NY samples (physeq Transform . macrobdella decora NY)
phyT.mdVT <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & AnimalSource=="MtSnowVT") # subset containing only Macrobdella decora VT samples (physeq Transform . macrobdella decora VT)
phyT.mdDaFilf <- subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA),Sample_Type=="ILF" & Da2F%in%c("none","Unk","0")) # merge fed Mdecora ILF samples (phyloseq Transformed . macrobdella decora Days after Feeding ILF)
phyT.mdILF <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples (physeq Transform . macrobdella decora ILF)
phyT.mdInt <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="Intestinum") # subset containing only Macrobdella decora intestinum samples (physeq Transform . macrobdella decora Intestinum)
phyT.mdGut <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type%in%c("ILF","Intestinum")) # subset containing only Macrobdella decora ILF and intestinum samples (physeq Transform . macrobdella decora Gut)

#####################################################################
########## Effect of Extraction Method  #############################
#####################################################################
### Only Mdecora ILF samples were extracted by more than one method ###
### Effecvt of extraction method on Mdecora ILF samples ###
sd.mdILF <- data.frame(sample_data(phyT.mdILF)) # separate sample data from phyloseq object (sample data . macrobdella decora ILF) 
bray.mdILF <- phyloseq::distance(phyT.mdILF, method = "bray") # calculate bray curtis distance matrix for Mdecora ILF samples (bray metric . macrobdella decora ILF) 
perm.md.exILFbray <- adonis(bray.mdILF ~ Da1Fb * Kit_Name, data = sd.mdILF, method = "bray") # Adonis test to determine effect of extraction kit on Mdecora ILF samples (permanova . macrobdella decora . extraction ILF Bray)
perm.md.exILFbray$aov.tab # print results
binary.mdILF <- phyloseq::distance(phyT.mdILF, method = "binary") # calculate binary Ochiai distance matrix for Mdecora ILF samples(binary . macrobdella decora ILF)
perm.md.exILFOchiai <- adonis(binary.mdILF ~ Da1Fb * Kit_Name, data = sd.mdILF, method = "binary") # Adonis test to determine effect of extraction kit on Mdecora ILF samples (permanova . macrobdella decora . extraction ILF Ochiai)
perm.md.exILFOchiai$aov.tab # print results

#####################################################################
########## Effect of Leech Species  #################################
#####################################################################
sd.leech <- data.frame(sample_data(phyT.leech)) # separate sample data from phyloseq object (sample data . leech)
bray.leech <- phyloseq::distance(phyT.leech, method = "bray") # calculate bray curtis distance matrix for leech samples (bray . leech)
perm.taxaBray <- adonis(bray.leech ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test to determine effect of leech species on samples (permanova . taxa Bray)
perm.taxaBray$aov.tab["Taxonomic_ID",] # print results
binary.leech <- phyloseq::distance(phyT.leech, method = "binary") # calculate Binary Ochiai distance matrix for leech samples (binary . leech)
perm.taxaOchiai <- adonis(binary.leech ~ Taxonomic_ID, data = sd.leech, method = "binary") # Adonis test to determine effect of leech species on samples (permanova . taxa Ochiai)
perm.taxaOchiai$aov.tab["Taxonomic_ID",] # print results

#####################################################################
########## TABLE S2 #################################################
#####################################################################
### ASVs grouped at various taxonomic levels
sd.leech <- data.frame(sample_data(phyT.leech)) # separate sample data from phyloseq object (sample data . leech)
perm.taxaPhylum <- adonis(phyloseq::distance(tax_glom(phyT.leech,taxrank="Phylum"), method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped at the Phylum level (permanova . taxa Phylum)
perm.taxaClass <- adonis(phyloseq::distance(tax_glom(phyT.leech,taxrank="Class"), method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped at the Class level (permanova . taxa Class)
perm.taxaOrder <- adonis(phyloseq::distance(tax_glom(phyT.leech,taxrank="Order"), method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped at the Order level(permanova . taxa Order)
perm.taxaFamily <- adonis(phyloseq::distance(tax_glom(phyT.leech,taxrank="Family"), method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped at the Family level(permanova . taxa Family)
perm.taxaGenus <- adonis(phyloseq::distance(tax_glom(phyT.leech,taxrank="Genus"), method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped at the Genus level (permanova . taxa Genus)
perm.taxaASV <- adonis(phyloseq::distance(phyT.leech, method = "bray") ~ Taxonomic_ID, data = sd.leech, method = "bray") # Adonis test of leech sample difference when taxa are grouped by ASV (permanova . taxa ASV)

perm.taxaPhylum$aov.tab["Taxonomic_ID",] # print results
perm.taxaClass$aov.tab["Taxonomic_ID",] # print results
perm.taxaOrder$aov.tab["Taxonomic_ID",] # print results
perm.taxaFamily$aov.tab["Taxonomic_ID",] # print results
perm.taxaGenus$aov.tab["Taxonomic_ID",] # print results
perm.taxaASV$aov.tab["Taxonomic_ID",] # print results
########## end TABLE S2 #############################################

#####################################################################
########## FIGURE 1 #################################################
#####################################################################
#### Host Species and Sample Type NMDS ###
dsNMDS <- phyT.leech  # define data for analysis (data set NMDS)
sample_data(dsNMDS)$Sample_Type = factor(sample_data(dsNMDS)$Sample_Type, levels = c("ILF","Intestinum","Bladder")) # define Sample_Type as a factor...because R 
mNMDS <- "bray" # define metric for analysis (metric NMDS)
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances (distance Ordination)
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination (ordination)
nmds.SpeciesType <- plot_ordination(dsNMDS, ord, shape="Sample_Type",color="spp") +
  labs(x="ASV",y=element_blank()) +
  scale_color_manual(values=pal.CB) +
  scale_shape_manual(values=c(16,17,3)) +
  theme_bw() +
  theme(text=element_text(family="Times New Roman", size=12),legend.position='top',legend.title = element_blank()) +
  guides(color = "none", size="none",fill="none")
nmds.SpeciesType$layers[[1]]<-NULL # remove layer with geom_points
nmds.SpeciesType<-nmds.SpeciesType + 
  geom_point(aes(fill=spp),size=2,alpha=0.5) +  
  stat_ellipse(type = "t", level = 0.95, linetype = 2) # re-write geom_point layer and add elipses
#nmds.SpeciesType # print NMDS plot
### Host Species and Sample Type by Order NMDS ###
dsNMDS.order <- tax_glom(dsNMDS,taxrank="Order") # define data for analysis (data set NMDS . order)
distOrdO = phyloseq::distance(dsNMDS.order, method = c(mNMDS)) # calculate distances (distance Ordination Order)
ordO = ordinate(dsNMDS.order, method = "NMDS", distance = distOrdO) # calculate ordination (ordination Order)
nmds.SpeciesTypeOrder <- plot_ordination(dsNMDS.order, ordO, shape="Sample_Type",color="spp") + 
  labs(x="Order",y=element_blank()) +
  scale_color_manual(values=pal.CB) +
  scale_shape_manual(values=c(16,17,3)) +
  theme_bw() +
  theme(text=element_text(family="Times", size=12),legend.position='top',legend.title = element_blank(),legend.text = element_text(face = "italic")) +
  guides(shape = "none")
nmds.SpeciesTypeOrder$layers[[1]]<-NULL # remove layer with geom_points
nmds.SpeciesTypeOrder<-nmds.SpeciesTypeOrder + 
  geom_point(aes(fill=spp),size=2, alpha=0.5) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) # re-write geom_point layer and add elipses
#nmds.SpeciesTypeOrder # print NMDS plot
ggsave(plot_grid(nmds.SpeciesType, nmds.SpeciesTypeOrder, labels = "AUTO",rel_widths = c(1,1), ncol=2), filename="NMDS/Mdfig1_plotNMDS.pdf", device=cairo_pdf, dpi="retina",width=6.87,height=4,units="in") # Save Figure 1 to pdf file
########## end FIGURE 1 #############################################

#####################################################################
########## TABLE 1 ##################################################
#####################################################################
### H.verbana variation by Organ, Feeding, Supplier, and Shipment lot ###
#phyT.hv <- subset_samples(phyT.leech, Taxonomic_ID=="Hverbana") # subset containing only Hirudo verbana samples (physeq Transform .  hirudo verbana)
sd.hv <- data.frame(sample_data(phyT.hv)) # separate sample data from phyloseq object (sample data . hirudo verbana)
dist.adonisHv <- phyloseq::distance(phyT.hv, method = "bray") # calculate bray curtis distance matrix for Hirudo verbana samples (distance . adonis Hirudo verbana)
perm.adonisHv <- adonis(dist.adonisHv ~ Sample_Type * AnimalSource * Da1Fb * ParentLot * MealLot1, data = sd.hv, method = "bray") # Adonis test of effect of Sample_Type, AnimalSource, ParentLot, Da1Fb, and MealLot1 on Hirudo verbana samples (permanova . adonis Hirudo verbana)
perm.adonisHv$aov.tab # print results
### M.decora variation by Organ, Feeding, Supplier, and Shipment lot
#phyT.md <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora") # subset containing only Macrobdella decora samples (physeq Transform . macrobdella decora)
sd.adonisMd <- data.frame(sample_data(phyT.md)) # separate sample data from phyloseq objec (sample data . macrobdella decora)
dist.adonisMd <- phyloseq::distance(phyT.md, method = "bray") # calculate bray curtis distance matrix for macrobdella decora samplea (distance . adonis Macrobdella decora)
perm.adonisMd <- adonis(dist.adonisMd ~ Sample_Type * Da1Fb * WildMonth * AnimalSource * ParentLot * MealLot1, data = sd.adonisMd, method = "bray") # Adonis test of effect of Sample_Type, WildMonth, Da1Fb, AnimalSource, ParentLot, and MealLot1 on Macrobdella decora samples (permanova . adonis Macrobdella decora)
perm.adonisMd$aov.tab # print results
########## end TABLE 1 ##############################################

#####################################################################
########## Define core and common ASVs  #############################
#####################################################################
#cMin <- var.medNeg/5000 # define minimum to count presence in a trasnformed sample
### Hirudo ###
### Removing Aeromonas and Mucinivorans contaminants from bladder reads (for determining true core/common ASVs) ###
phyR.hv<-subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%c(sample_names(phyT.hv))) # create raw reads phyloseq object of Hirduo verbana samples (physeq Raw . hirudo verbana)
phyR.hvBlad2<-subset_samples(phyR.hv,Sample_Type=="Bladder") # keep only bladder samples (physeq Raw . hirudo verbana Bladder 2)
phyR.hvBlad<-subset_taxa(phyR.hvBlad2,!Genus%in%c("Mucinivorans","Aeromonas","Bacteroides")) # remove Mucinivorans, Bacteroides, and Aeromonas ASVs (physeq Raw . hirudo verbana Bladder)
phyT.hvBlad<-transform_sample_counts(phyR.hvBlad, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . hirudo verbana Bladder)
fm.hvBlad = fast_melt(phyT.hvBlad) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana Bladder)
prevdt.hvBlad = fm.hvBlad[,list(Prevalence = sum(count >= cMin), 
                                TotalCount = sum(count),
                                MinCount = min(count),
                                MaxCount = max(count),
                                MedCount = median(count),
                                Med2Count = median(count[count!=0]),
                                PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana Bladder)

#phyT.hvInt<-subset_samples(phyT.leech, Taxonomic_ID=="Hverbana" & Sample_Type=="Intestinum") # subset containing only Hirudo verbana intestinum samples (physeq Transform . Hirudo verbana intestinum)
fm.hvInt = fast_melt(phyT.hvInt) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana Intestinum)
prevdt.hvInt = fm.hvInt[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.hvInt)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table Hirudo verbana intestinum)

#phyT.hvILF<-subset_samples(phyT.leech, Taxonomic_ID=="Hverbana" & Sample_Type=="ILF") # subset containing only Hirudo verbana ILF samples (physeq Transform . Hirudo verbana ILF)
fm.hvILF = fast_melt(phyT.hvILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF)
prevdt.hvILF = fm.hvILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.hvILF)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table Hirudo verbana ILF)

phyT.hvILF0 <- subset_samples(phyT.hvILF,Da1Fb%in%c("0")) # subset containing only Hirudo verbana ILF samples from unfed animals (physeq Transform . hirudo verbana ILF 0)
fm.hvILF0 = fast_melt(phyT.hvILF0) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF 0)
prevdt.hvILF0 = fm.hvILF0[, list(Prevalence = sum(count >= cMin), 
                                 TotalPer = sum(count),
                                 MinCount = min(count),
                                 MaxCount = max(count),   
                                 MedCount = median(count),   
                                 Med2Count = median(count[count!=0]),   
                                 PrevPer = sum(count >= cMin) / nsamples(phyT.hvILF0)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana ILF 0)

#phyT.hvILFbbez <- subset_samples(phyT.hvILF,AnimalSource=="BBEZ") # subset containing only Hirudo verbana ILF samples from animals from supplier 1 (physeq Transform . hirudo verbana ILF bbez)
fm.hvILFbbez = fast_melt(phyT.hvILFbbez) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF bbez)
prevdt.hvILFbbez = fm.hvILFbbez[, list(Prevalence = sum(count >= cMin), 
                                       TotalPer = sum(count),
                                       MinCount = min(count),
                                       MaxCount = max(count),   
                                       MedCount = median(count),   
                                       Med2Count = median(count[count!=0]),   
                                       PrevPer = sum(count >= cMin) / nsamples(phyT.hvILFbbez)),
                                by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana ILF bbez)

phyT.hvILFbbez0<-subset_samples(phyT.hvILF,AnimalSource=="BBEZ" & Da1Fb%in%c("0","90+")) # subset containing only Hirudo verbana ILF samples from unfed animals from supplier 1 (physeq Transform . hirudo verbana ILF bbez 0)
fm.hvILFbbez0 = fast_melt(phyT.hvILFbbez0) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF bbez 0)
prevdt.hvILFbbez0 = fm.hvILFbbez0[, list(Prevalence = sum(count >= cMin), 
                                         TotalPer = sum(count),
                                         MinCount = min(count),
                                         MaxCount = max(count),   
                                         MedCount = median(count),   
                                         Med2Count = median(count[count!=0]),   
                                         PrevPer = sum(count >= cMin) / nsamples(phyT.hvILFbbez0)),
                                  by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana ILF bbez 0)

#phyT.hvILFusa <- subset_samples(phyT.hvILF,AnimalSource=="USA") # subset containing only Hirudo verbana ILF samples from animals from supplier 2 (physeq Transform . hirudo verbana ILF usa)
fm.hvILFusa = fast_melt(phyT.hvILFusa) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF usa)
prevdt.hvILFusa = fm.hvILFusa[, list(Prevalence = sum(count >= cMin), 
                                     TotalPer = sum(count),
                                     MinCount = min(count),
                                     MaxCount = max(count),   
                                     MedCount = median(count),   
                                     Med2Count = median(count[count!=0]),   
                                     PrevPer = sum(count >= cMin) / nsamples(phyT.hvILFusa)),
                              by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana ILF usa)

phyT.hvILFusa0<-subset_samples(phyT.hvILF,AnimalSource=="USA" & Da1F%in%c("0","90+")) # subset containing only Hirudo verbana ILF samples from unfed animals from supplier 2 (physeq Transform . hirudo verbana ILF usa 0)
fm.hvILFusa0 = fast_melt(phyT.hvILFusa0) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana ILF usa 0)
prevdt.hvILFusa0 = fm.hvILFusa0[, list(Prevalence = sum(count >= cMin), 
                                       TotalPer = sum(count),
                                       MinCount = min(count),
                                       MaxCount = max(count),   
                                       MedCount = median(count),   
                                       Med2Count = median(count[count!=0]),   
                                       PrevPer = sum(count >= cMin) / nsamples(phyT.hvILFusa0)),
                                by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . hirudo verbana ILF usa 0) 

### Macrobdella ###
### Removing Aeromonas and Bacteroides contaminants from bladder reads (for determining true core/common ASVs) ###
phyR.md<-subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%c(sample_names(phyT.md))) # create raw reads phyloseq object of Macrobdella decora samples (physeq Raw . macrobdella decora)
phyR.mdBlad2<-subset_samples(phyR.md,Sample_Type=="Bladder") # keep only bladder samples (physeq Raw . macrobdella decora Bladder 2)
phyR.mdBlad<-subset_taxa(phyR.mdBlad2,!Genus%in%c("Mucinivorans","Aeromonas","Bacteroides")) # remove Mucinivorans, Bacteroides, and Aeromonas ASVs (physeq Raw . macrobdella decora Bladder)
phyT.mdBlad<-transform_sample_counts(phyR.mdBlad, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . macrobdella decora Bladder)
fm.mdBlad = fast_melt(phyT.mdBlad) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . macrobdella decora bladder)
prevdt.mdBlad = fm.mdBlad[,list(Prevalence = sum(count >= .001), 
                                TotalCount = sum(count),
                                MinCount = min(count),
                                MaxCount = max(count),
                                MedCount = median(count),
                                Med2Count = median(count[count!=0]),
                                PrevPer = sum(count >= cMin) / nsamples(phyT.mdBlad)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . macrobdella decora Bladder)

fm.mdInt = fast_melt(phyT.mdInt) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . macrobdella decora Intestinum)
prevdt.mdInt = fm.mdInt[, list(Prevalence = sum(count >= cMin),
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.mdInt)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . macrobdella decora Intestinum)

fm.mdILF = fast_melt(phyT.mdILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . macrobdella decora ILF)
prevdt.mdILF = fm.mdILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.mdILF)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . macrobdella decora ILF)

phyT.mdILF0<-subset_samples(phyT.mdILF,Da1Fb=="0") # subset containing only Macrobdella decora ILF samples from unfed animals (physeq Transform . macrobdella decora ILF 0)
fm.mdILF0 = fast_melt(phyT.mdILF0) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . macrobdella decora ILF 0)
prevdt.mdILF0 = fm.mdILF0[, list(Prevalence = sum(count >= cMin), 
                                 TotalPer = sum(count),
                                 MinCount = min(count),
                                 MaxCount = max(count),   
                                 MedCount = median(count),   
                                 Med2Count = median(count[count!=0]),   
                                 PrevPer = sum(count >= cMin) / nsamples(phyT.mdILF0)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . macrobdella decora ILF 0)

### Connecticut Macrobdella decora ###
phyT.ctBlad<-subset_samples(phyT.mdBlad,AnimalSource=="Wlot") # subset containing only Macrobdella decora bladder samples from Connecticut animals (physeq Transform . connecticut Bladder)
fm.ctBlad = fast_melt(phyT.ctBlad) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . connecticut Bladder)
prevdt.ctBlad = fm.ctBlad[, list(Prevalence = sum(count >= cMin), 
                                 TotalPer = sum(count),
                                 MinCount = min(count),
                                 MaxCount = max(count),   
                                 MedCount = median(count),   
                                 Med2Count = median(count[count!=0]),   
                                 PrevPer = sum(count >= cMin) / nsamples(phyT.ctBlad)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . connecticut Bladder)

phyT.ctInt<-subset_samples(phyT.mdCT,Sample_Type=="Intestinum") # subset containing only Macrobdella decora intestinum samples from Connecticut animals (physeq Transform . connecticut Intestinum)
fm.ctInt = fast_melt(phyT.ctInt) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . connecticut Intestinum)
prevdt.ctInt = fm.ctInt[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.ctInt)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . connecticut Intestinum)

phyT.ctILF<-subset_samples(phyT.mdCT,Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples from Connecticut animals (physeq Transform . connecticut ILF)
fm.ctILF = fast_melt(phyT.ctILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . connecticut ILF)
prevdt.ctILF = fm.ctILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.ctILF)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . connecticut ILF)

### Massachusetts Macrobdella decora ###
phyT.maBlad<-subset_samples(phyT.mdBlad,AnimalSource=="GrotonMA") # subset containing only Macrobdella decora ILF samples from Massachusetts animals (physeq Transform . massachusetts Bladder)
fm.maBlad = fast_melt(phyT.maBlad) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . massachusetts Bladder)
prevdt.maBlad = fm.maBlad[, list(Prevalence = sum(count >= cMin), 
                                 TotalPer = sum(count),
                                 MinCount = min(count),
                                 MaxCount = max(count),   
                                 MedCount = median(count),   
                                 Med2Count = median(count[count!=0]),   
                                 PrevPer = sum(count >= cMin) / nsamples(phyT.maBlad)),
                          by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . massachusetts Bladder)

phyT.maInt<-subset_samples(phyT.mdMA,Sample_Type=="Intestinum") # subset containing only Macrobdella decora ILF samples from Massachusetts animals (physeq Transform . massachusetts Intestinum)
fm.maInt = fast_melt(phyT.maInt) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . massachusetts Intestinum)
prevdt.maInt = fm.maInt[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.maInt)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . massachusetts Intestinum)

phyT.maILF<-subset_samples(phyT.mdMA,Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples from Massachusetts animals (physeq Transform . massachusetts ILF)
fm.maILF = fast_melt(phyT.maILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . massachusetts ILF)
prevdt.maILF = fm.maILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.maILF)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . massachusetts ILF)

### New York Macrobdella decora ###
phyT.nyILF<-subset_samples(phyT.mdNY,Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples from New York animals (physeq Transform . new york ILF)
fm.nyILF = fast_melt(phyT.nyILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . new york ILF)
prevdt.nyILF = fm.nyILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.nyILF)),
                        by = ASV] # make simple table listing ASV, Prevalence, TotalCount, MinCount, MaxCount, MedCount, Med2Count, and PrevPer (prevalence data table . new york ILF)

########## List core ASVs ###########################################
#cpR<-c(.9) # define the level for determining core ASVs, reported as a percent of the samples tested (core percent)
#cMin <- var.medNeg/5000 # define minimum to count presence in a trasnformed sample (count Minimum)
### Hirudo verbana
ls.coreHvILF0 = prevdt.hvILF0[(Prevalence >= cpR*nsamples(phyT.hvILF0) & MaxCount >= cMin), ASV] # make list of core ASVs for unfed Hirudo verbana ILF samples (list . core Hirudo verbana ILF 0)
ls.coreHvILF = prevdt.hvILF[(Prevalence >= cpR*nsamples(phyT.hvILF) & MaxCount >= cMin), ASV] # make list of core ASVs for Hirudo verbana ILF samples (list . core Hirudo verbana ILF)
ls.coreHvILFbbez = prevdt.hvILFbbez[(Prevalence >= cpR*nsamples(phyT.hvILFbbez) & MaxCount >= cMin), ASV] # make list of core ASVs for Hirudo verbana ILF of animals from supplier 1 samples (list . core Hirudo verbana ILF bbez)
ls.coreHvILFusa = prevdt.hvILFusa[(Prevalence >= cpR*nsamples(phyT.hvILFusa) & MaxCount >= cMin), ASV] # make list of core ASVs for Hirudo verbana ILF of animals from supplier 2 samples (list . core Hirudo verbana ILF usa)
ls.coreHvInt = prevdt.hvInt[(Prevalence >= cpR*nsamples(phyT.hvInt) & MaxCount >= cMin), ASV] # make list of core ASVs for Hirudo verbana Intestinum samples (list . core Hirudo verbana Intestinum)
ls.coreHvBlad = prevdt.hvBlad[(Prevalence >= cpR*nsamples(phyT.hvBlad) & MaxCount >= cMin), ASV] # make list of core ASVs for Hirudo verbana Bladder samples (list . core Hirudo verbana Bladder)
ls.coreHv = unique(c(ls.coreHvBlad,ls.coreHvILF,ls.coreHvInt,ls.coreHvILF0)) # compile list of core ASVs for Hirudo verbana samples (list . core Hirudo verbana)
### CT Macrobdella decora
ls.coreCtILF = prevdt.ctILF[(Prevalence >= cpR*nsamples(phyT.ctILF) & MaxCount >= cMin), ASV] # make list of core ASVs for Connecticut Macrobdella decora ILF samples (list . core Connecticut ILF)
ls.coreCtInt = prevdt.ctInt[(Prevalence >= cpR*nsamples(phyT.ctInt) & MaxCount >= cMin), ASV] # make list of core ASVs for Connecticut Macrobdella decora intestinum samples (list . core Connecticut Intestinum)
ls.coreCtBlad = prevdt.ctBlad[(Prevalence >= cpR*nsamples(phyT.ctBlad) & MaxCount >= cMin), ASV] # make list of core ASVs for Connecticut Macrobdella decora bladder samples (list . core Connecticut Bladder)
ls.coreCT = unique(c(ls.coreCtBlad,ls.coreCtILF,ls.coreCtInt)) # compile list of core ASVs for Connecticut Macrobdella decora samples (list . core Connecticut)
### MA Macrobdella decora
ls.coreMaILF = prevdt.maILF[(Prevalence >= cpR*nsamples(phyT.maILF) & MaxCount >= cMin), ASV] # make list of core ASVs for Massachusetts Macrobdella decora ILF samples (list . core Massachusetts ILF)
ls.coreMaInt = prevdt.maInt[(Prevalence >= cpR*nsamples(phyT.maInt) & MaxCount >= cMin), ASV] # make list of core ASVs for Massachusetts Macrobdella decora intestinum samples (list . core Massachusetts Intestinum)
ls.coreMaBlad = prevdt.maBlad[(Prevalence >= cpR*nsamples(phyT.maBlad) & MaxCount >= cMin), ASV] # make list of core ASVs for Massachusetts Macrobdella decora bladder samples (list . core Massachusetts Bladder)
ls.coreMA = unique(c(ls.coreMaBlad,ls.coreMaILF,ls.coreMaInt)) # compile list of core ASVs for Massachusetts Macrobdella decora samples (list . core Massachusetts)
### NY Macrobdella decora
ls.coreNyILF = prevdt.nyILF[(Prevalence >= cpR*nsamples(phyT.nyILF) & MaxCount >= cMin), ASV] # make list of core ASVs for New York Macrobdella decora ILF samples (list . core New york ILF)
ls.coreNY = unique(c(ls.coreNyILF)) # compile list of core ASVs for New York Macrobdella decora samples (list . core New York)
### Compiled Macrobdella decora
ls.coreMdILF = prevdt.mdILF[(Prevalence >= cpR*nsamples(phyT.mdILF) & MaxCount >= cMin), ASV] # make list of core ASVs for Macrobdella decora ILF samples (list . core Macrobdella decora ILF)
ls.coreMdILF0 = prevdt.mdILF0[(Prevalence >= cpR*nsamples(phyT.mdILF0) & MaxCount >= cMin), ASV] # make list of core ASVs for Macrobdella decora ILF samples from unfed animals (list . core Macrobdella decora ILF 0)
ls.coreMdInt = prevdt.mdInt[(Prevalence >= cpR*nsamples(phyT.mdInt) & MaxCount >= cMin), ASV] # make list of core ASVs for Macrobdella decora intestinum samples (list . core Macrobdella decora Intestinum)
ls.coreMdBlad = prevdt.mdBlad[(Prevalence >= cpR*nsamples(phyT.mdBlad) & MaxCount >= cMin), ASV] # make list of core ASVs for Macrobdella decora bladder samples (list . core Macrobdella decora Bladder)
ls.coreMd = unique(c(ls.coreCT,ls.coreMA)) # compile list of core ASVs for Massachusetts and Connecticut Macrobdella decora samples (list . core Macrobdella decora)
ls.coreDaF<-unique(c(ls.coreMdILF,ls.coreMdInt,ls.coreHvILF0,ls.coreHvILF)) # compile list of core ASVs for Hirudo verbana and Macrobdella decora ILF and intestinum samples to track after feeding (list . core Days after Feeding)
ls.coreGut <- unique(c(ls.coreMdILF,ls.coreMdInt,ls.coreHvILF,ls.coreHvInt)) # compile list of common ASVs for Hirudo verbana and Macrobdella decora ILF and intestinum samples (list . common Gut)
ls.coreTot = unique(c(ls.coreMd,ls.coreHv,ls.coreCT,ls.coreMA,ls.coreMdILF,ls.coreMdILF0,ls.coreMdInt,ls.coreMdBlad)) # compile list of core ASVs for Hirudo verbana and Macrobdella decora samples (list . core Total)

########## List common ASVs #########################################
#cMin <- var.medNeg/5000 # define minimum to count presence in a trasnformed sample (count Minimum)
#cpM<-c(.7) # define the level for determining common ASVs, reported as a percent of the samples tested (core percent)
### Hirudo verbana
ls.comHvILF0 = prevdt.hvILF0[(Prevalence >= cpM*nsamples(phyT.hvILF0) & MaxCount >= cMin), ASV] # make list of common ASVs for unfed Hirudo verbana ILF samples (list . common Hirudo verbana ILF 0)
ls.comHvILF = prevdt.hvILF[(Prevalence >= cpM*nsamples(phyT.hvILF) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana ILF samples (list . common Hirudo verbana ILF)
ls.comHvILFbbez = prevdt.hvILFbbez[(Prevalence >= cpM*nsamples(phyT.hvILFbbez) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana ILF of animals from supplier 1 samples (list . common Hirudo verbana ILF bbez)
ls.comHvILFusa = prevdt.hvILFusa[(Prevalence >= cpM*nsamples(phyT.hvILFusa) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana ILF of animals from supplier 2 samples (list . common Hirudo verbana ILF usa)
ls.comHvILFusa0 = prevdt.hvILFusa0[(Prevalence >= cpM*nsamples(phyT.hvILFusa0) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana ILF of unfed animals from supplier 2 samples (list . common Hirudo verbana ILF usa 0)
ls.comHvILFbbez0 = prevdt.hvILFbbez0[(Prevalence >= cpM*nsamples(phyT.hvILFbbez0) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana ILF of unfed animals from supplier 1 samples (list . common Hirudo verbana ILF bbez 0)
ls.comHvInt = prevdt.hvInt[(Prevalence >= cpM*nsamples(phyT.hvInt) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana Intestinum samples (list . common Hirudo verbana Intestinum)
ls.comHvBlad = prevdt.hvBlad[(Prevalence >= cpM*nsamples(phyT.hvBlad) & MaxCount >= cMin), ASV] # make list of common ASVs for Hirudo verbana Bladder samples (list . common Hirudo verbana Bladder)
ls.comHv = unique(c(ls.comHvBlad,ls.comHvILF,ls.comHvInt,ls.comHvILF0)) # compile list of common ASVs for Hirudo verbana samples (list . common Hirudo verbana)
### CT Macrobdella decora
ls.comCtILF = prevdt.ctILF[(Prevalence >= cpM*nsamples(phyT.ctILF) & MaxCount >= cMin), ASV] # make list of common ASVs for Connecticut Macrobdella decora ILF samples (list . common Connecticut ILF)
ls.comCtInt = prevdt.ctILF[(Prevalence >= cpM*nsamples(phyT.ctInt) & MaxCount >= cMin), ASV] # make list of common ASVs for Connecticut Macrobdella decora intestinum samples (list . common Connecticut Intestinum)
ls.comCtBlad = prevdt.ctBlad[(Prevalence >= cpM*nsamples(phyT.ctBlad) & MaxCount >= cMin), ASV] # make list of common ASVs for Connecticut Macrobdella decora bladder samples (list . common Connecticut Bladder)
ls.comCT = unique(c(ls.comCtBlad,ls.comCtILF,ls.comCtInt)) # compile list of common ASVs for Connecticut Macrobdella decora samples (list . common Connecticut)
### MA Macrobdella decora
ls.comMaILF = prevdt.maILF[(Prevalence >= cpM*nsamples(phyT.maILF) & MaxCount >= cMin), ASV] # make list of common ASVs for Massachusetts Macrobdella decora ILF samples (list . common Massachusetts ILF)
ls.comMaInt = prevdt.maInt[(Prevalence >= cpM*nsamples(phyT.maInt) & MaxCount >= cMin), ASV] # make list of common ASVs for Massachusetts Macrobdella decora intestinum samples (list . common Massachusetts Intestinum)
ls.comMaBlad = prevdt.maBlad[(Prevalence >= cpM*nsamples(phyT.mdBlad) & MaxCount >= cMin), ASV] # make list of common ASVs for Massachusetts Macrobdella decora bladder samples (list . common Massachusetts Bladder)
ls.comMA = unique(c(ls.comMaBlad,ls.comMaILF,ls.comMaInt)) # compile list of common ASVs for Massachusetts Macrobdella decora samples (list . common Massachusetts)
### NY Macrobdella decora
ls.comNyILF = prevdt.nyILF[(Prevalence >= cpM*nsamples(phyT.nyILF) & MaxCount >= cMin), ASV] # make list of common ASVs for New York Macrobdella decora ILF samples (list . common New york ILF)
ls.comNY = unique(c(ls.comNyILF)) # compile list of common ASVs for New York Macrobdella decora samples (list . common New York)
### Compiled Macrobdella decora
ls.comMdILF = prevdt.mdILF[(Prevalence >= cpM*nsamples(phyT.mdILF) & MaxCount >= cMin), ASV] # make list of common ASVs for Macrobdella decora ILF samples (list . common Macrobdella decora ILF)
ls.comMdILF0 = prevdt.mdILF0[(Prevalence >= cpM*nsamples(phyT.mdILF0) & MaxCount >= cMin), ASV] # make list of common ASVs for Macrobdella decora ILF samples from unfed animals (list . common Macrobdella decora ILF 0)
ls.comMdInt = prevdt.mdInt[(Prevalence >= cpM*nsamples(phyT.mdInt) & MaxCount >= cMin), ASV] # make list of common ASVs for Macrobdella decora intestinum samples (list . common Macrobdella decora Intestinum)
ls.comMdBlad = prevdt.mdBlad[(Prevalence >= cpM*nsamples(phyT.mdBlad) & MaxCount >= cMin), ASV] # make list of common ASVs for Macrobdella decora bladder samples (list . common Macrobdella decora Bladder)
ls.comMd = unique(c(ls.comCT,ls.comMA)) # compile list of common ASVs for Massachusetts and Connecticut Macrobdella decora samples (list . common Macrobdella decora)
ls.comDaF <- unique(c(ls.comCtILF,ls.comMaILF,ls.comMaInt,ls.comHvILFbbez, ls.comHvILFusa)) # compile list of common ASVs for Hirudo verbana and Macrobdella decora ILF and intestinum samples to track after feeding (list . common Days after Feeding)
ls.comGut <- unique(c(ls.comMdILF,ls.comMdInt,ls.comHvILF,ls.comHvInt)) # compile list of common ASVs for Hirudo verbana and Macrobdella decora ILF and intestinum samples (list . common Gut)
ls.comTot = unique(c(ls.comMd,ls.comHv,ls.comCT,ls.comMA,ls.comMdILF,ls.comMdILF0,ls.comMdInt,ls.comMdBlad)) # compile list of common ASVs for Hirudo verbana and Macrobdella decora samples (list . common Total)

#####################################################################
########## FIGURE 2 #################################################
#####################################################################
phyR.leech <- subset_samples(phyR.pru2, sample_names(phyR.sin)%in%c(sample_names(phyT.leech))) # create raw reads phyloseq object of leech samples (physeq Raw . leech)
phyR.ILF <- subset_samples(phyR.leech, Sample_Type%in%c("ILF")) # subset containing only ILF samples (physeq Raw . ILF)
phyR.Int <- subset_samples(phyR.leech, Sample_Type=="Intestinum") # subset containing only intestinum samples (physeq Raw . Intestinum)
########## Statistics ##########
sd.ILF <- data.frame(sample_data(phyR.ILF)) # separate sample data from phyloseq object (sample data . ILF)
bray.ILF <- phyloseq::distance(phyR.ILF, method = "bray") # calculate bray curtis distance matrix for ILF samples (bray . ILF)
perm.taxILFbray <- adonis(bray.ILF ~ Taxonomic_ID, data = sd.ILF, method = "bray") # Adonis test to determine effect of leech species (permanova . taxonomy ILF bray)
perm.taxILFbray # print results
binary.ILF <- phyloseq::distance(phyR.ILF, method = "binary") # calculate Binary Ochiai distance matrix for ILF samples (binary . ILF)
perm.taxILFochiai <- adonis(binary.ILF ~ Taxonomic_ID, data = sd.ILF, method = "binary") # Adonis test to determine effect of leech species (permanova . taxonomy ILF ochiai)
perm.taxILFochiai # print results

########## Alpha Diversity ILF ##########
### Data must not be transformed for alpha diversity calculations
phyR.adGut <- subset_samples(phyR.leech,!Da2F%in%c("1","2","4","7","8","28") & !Da1F%in%c("1","2","4","7") & Sample_Type%in%c("ILF","Intestinum")) # subset containing only ILF and intestinum samples that are less than 1 or more than 30 days after feeding (Phyloseq raw . alpha diversity Gut)) 
my_comp <- list( c("Hverbana ILF", "Mdecora ILF"), c("Hverbana ILF", "Hverbana Intestinum"),c("Hverbana ILF", "Mdecora Intestinum")) # list of comparisons to make (my_comparisons)
### Shannon diversity sample type comparison ###
pDiv.adGutSh <- plot_richness(phyR.adGut, x = "TaxType", measures=c("Shannon")) +
  geom_boxplot(aes(fill=spp)) +
  labs(y="Shannon Diversity") +
  stat_compare_means(comparisons = my_comp,label="p.signif", p.adjust.method = "bonferroni", method = "wilcox.test") +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,3)) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none",strip.text.x=element_blank()) +
  scale_x_discrete(labels=c("Hverbana ILF" = exp.HvILF, "Hverbana Intestinum" = exp.HvInt, "Mdecora ILF" = exp.MdILF, "Mdecora Intestinum" = exp.MdInt)) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
pDiv.adGutSh$layers<-pDiv.adGutSh$layers[-1] # remove layer with underlying geom_points
#pDiv.adGutSh  # uncomment to print plot 
ggsave(pDiv.adGutSh,filename="Plots/Mdfig2_plotAlpha.pdf", dpi="retina",width=6.87,units="in") # Save FIGURE 2 to .pdf file

### Inverse Simpson diversity sample type comparison ###
pDiv.adGutSi <- plot_richness(phyR.adGut, x = "TaxType", measures=c("InvSimpson")) +
  geom_boxplot(aes(fill=spp)) +
  stat_compare_means(comparisons = my_comp,label="p.signif", p.adjust.method = "bonferroni", method = "wilcox.test") +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,10)) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none",strip.text.x=element_blank()) +
  scale_x_discrete(labels=c("Hverbana ILF" = exp.HvILF, "Hverbana Intestinum" = exp.HvInt, "Mdecora ILF" = exp.MdILF, "Mdecora Intestinum" = exp.MdInt)) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
pDiv.adGutSi$layers<-pDiv.adGutSi$layers[-1] # remove layer with underlying geom_points
#pDiv.adGutSi  # uncomment to print plot 
### Chao diversity sample type comparison ###
pDiv.adGutCh <- plot_richness(phyR.adGut, x = "TaxType", measures=c("Chao1")) +
  geom_boxplot(aes(fill=spp)) +
  stat_compare_means(comparisons = my_comp,label="p.signif", p.adjust.method = "bonferroni", method = "wilcox.test") +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,40)) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none",strip.text.x=element_blank()) +
  scale_x_discrete(labels=c("Hverbana ILF" = exp.HvILF, "Hverbana Intestinum" = exp.HvInt, "Mdecora ILF" = exp.MdILF, "Mdecora Intestinum" = exp.MdInt)) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
pDiv.adGutCh$layers<-pDiv.adGutCh$layers[-1] # remove layer with underlying geom_points
#pDiv.adGutCh  # uncomment to print plot 
### Observed ASVs diversity sample type comparison ###
pDiv.adGutO <- plot_richness(phyR.adGut, x = "TaxType", measures=c("Observed")) +
  geom_boxplot(aes(fill=spp)) +
  stat_compare_means(comparisons = my_comp,label="p.signif", p.adjust.method = "bonferroni", method = "wilcox.test") +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,40)) +
  theme_bw() +
  theme(text=element_text(size=12),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none",strip.text.x=element_blank()) +
  scale_x_discrete(labels=c("Hverbana ILF" = exp.HvILF, "Hverbana Intestinum" = exp.HvInt, "Mdecora ILF" = exp.MdILF, "Mdecora Intestinum" = exp.MdInt)) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
pDiv.adGutO$layers<-pDiv.adGutO$layers[-1] # remove layer with underlying geom_points
#pDiv.adGutO  # uncomment to print plot 
########## end FIGURE 2 #############################################

########## Shifts in abundance v. gut location ###########################
#phyT.hvGut <- subset_samples(phyT.hv, Sample_Type%in%c("ILF","Intestinum")) # subset containing only Hirudo verbana ILF and intestinum samples (physeq Transform . hirudo verbana Gut)
phyR.hvGut <- subset_samples(phyR.pru2, sample_names(phyR.pru2)%in%c(sample_names(phyT.hvGut))) # create raw reads phyloseq object of Hirudo verbana gut samples (physeq Raw . hirudo verbana Gut)
phyR.hvGutTrim <- subset_taxa(phyR.hvGut, taxa_sums(phyR.hvGut)>1) # subset containing only taxa present at > cMin (physeq Raw . hirudo verbana Gut Trim)
fm.hvGut = fast_melt(phyR.hvGutTrim) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana Gut)
prevdt.hvGut = fm.hvGut[, list(Prevalence = sum(count >= var.medNeg)), by = ASV] # make simple table listing ASV and Prevalence (prevalence data table . Hirudo verbana Gut)
ls.keepHvGut = prevdt.hvGut[(Prevalence > 2), ASV] # make list of ASVs present in more than 1 sample (list . keep Hirudo verbana Gut)

dsq.start = phyloseq_to_deseq2(phyR.hvGutTrim, ~ Sample_Type) 
geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.1
sigtab = res[(res$padj < alpha & res$baseMean > var.medNeg & res$padj < 0.01), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.hvGutTrim)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 1, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","ASV")]
ls.02<-levels(posigtab$ASV)
sigtabgen = sigtab

ggplot(posigtab, aes(y=Genus, x=log2FoldChange, color=Order)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=pal.pairBiome) 

write.table(sigtabgen, "table_HvGutChange_tbl3data.csv", sep=",",row.names=FALSE,col.names=TRUE) # export table to csv, column names included

#####################################################################
########## Effect of M.decora Collection Month ######################
#####################################################################
phyT.mdSeason <- subset_samples(phyT.md, AnimalSource%in%c("Wlot","GrotonMA")) # subset containing only Macrobdella decora samples from Connecticut and Massachusetts (physeq Transform . macrobdella decora Season)
phyT.mdSeasILF <- subset_samples(phyT.mdSeason, Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples from Connecticut and Massachusetts (physeq Transform . macrobdella decora Season ILF)
phyT.mdSeasInt <- subset_samples(phyT.mdSeason, Sample_Type=="Intestinum") # subset containing only Macrobdella decora intestinum samples from Connecticut and Massachusetts (physeq Transform . macrobdella decora Season Intestinum)

bray.seasILF <- phyloseq::distance(phyT.mdSeasILF, method = "bray") # calculate bray curtis distance matrix (bray . season ILF)
sd.seasILF <- data.frame(sample_data(phyT.mdSeasILF)) # separate sample data from phyloseq object (sample data . season ILF)
pairwise.adonis(bray.seasILF,sd.seasILF$Season) # pairwise Adonis calculation comparing seasons for ILF samples

perm.md.seasILF <- adonis(bray.seasILF ~ Season, data = sd.seasILF, method = "bray") # Adonis test to determine effect of extraction kit on Mdecora ILF samples (permanova . macrobdella decora . extraction ILF Ochiai)
perm.md.seasILF$aov.tab

bray.seasInt <- phyloseq::distance(phyT.mdSeasInt, method = "bray")  # calculate bray curtis distance matrix (bray . season Intestinum)
sd.seasInt <- data.frame(sample_data(phyT.mdSeasInt)) # separate sample data from phyloseq object (sample data . season Intestinum)
pairwise.adonis(bray.seasInt,sd.seasInt$Season) # pairwise Adonis calculation comparing seasons for intestinum samples

#####################################################################
########## FIGURE 4 #################################################
#####################################################################
### Collection Month NMDS - ILF ###
dsNMDS.mdSeasILF <- phyT.mdSeasILF # define data for analysis
mNMDS.Seas <- "binary" # define metric for analysis
distOrdILF = phyloseq::distance(dsNMDS.mdSeasILF, method = c(mNMDS.Seas)) # calculate distances
ordILF = ordinate(dsNMDS.mdSeasILF, method = "NMDS", distance = distOrdILF) # calculate ordination
nmds.collILF <- plot_ordination(dsNMDS.mdSeasILF, ordILF, color = "Season") + 
  geom_point(aes(fill=Season),size=2) +
  labs(x="ILF",y=element_blank()) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  theme_bw() +
  theme(text=element_text(size=12),legend.position="none")
#nmds.collILF # uncomment to print plot
nmds.legend <- plot_ordination(dsNMDS.mdSeasILF, ordILF, color = "Season") + 
  scale_color_manual(values=pal.CB) +
  theme_bw() +
  theme(text=element_text(size=12),legend.position='top',legend.title = element_blank())
### Collection Month NMDS - Intestinum ###
dsNMDS.mdSeasInt <- phyT.mdSeasInt # define data for analysis
distOrdInt = phyloseq::distance(dsNMDS.mdSeasInt, method = c(mNMDS.Seas)) # calculate distances
ordInt = ordinate(dsNMDS.mdSeasInt, method = "NMDS", distance = distOrdInt) # calculate ordination
nmds.collInt <- plot_ordination(dsNMDS.mdSeasInt, ordInt, color = "Season") + 
  geom_point(aes(fill=Season),size=2) +
  labs(x="Intestinum",y=element_blank()) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  theme_bw() +
  theme(text=element_text(size=12),legend.position="none")
#nmds.collInt # uncomment to print plot
ggsave(plot_grid(plot_grid(nmds.collILF, nmds.collInt, labels = "AUTO"),cowplot::get_legend(nmds.legend),nrow=2,rel_heights = c(5, 1)), filename="NMDS/Mdfig4_plotNMDS.pdf", device=cairo_pdf, dpi="retina",width=6.87,units="in") # Save FIGURE 4 to file
########## end FIGURE 4 #############################################

### Adonis analysis of ILF and intestinum samples
#phyT.mdGut <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type%in%c("ILF","Intestinum")) # subset containing only Macrobdella decora ILF and intestinum samples (physeq Transform . macrobdella decora Gut)
bray.mdGut <- phyloseq::distance(phyT.mdGut, method = "bray") # calculate bray curtis distance matrix (bray . macrobdella decora Gut)
sd.mdGut <- data.frame(sample_data(phyT.mdGut)) # separate sample data from phyloseq object (sample data . macrobdella decora Gut)
pairwise.adonis(bray.mdGut,paste(sd.mdGut$Sample_Type, sd.mdGut$AnimalSource)) # pairwise Adonis calculation comparing combination of Sample_Type and AnimalSoure for Macrobdella decora gut samples
### Adonis analysis of just CT and MA samples
bray.mdGutTrim<- phyloseq::distance(subset_samples(phyT.mdGut, AnimalSource%in%c("Wlot","GrotonMA")), method = "bray") # calculate bray curtis distance matrix for Macrobdella decora samples from CT and MA (bray . macrobdella decora Gut Trim)
sd.mdGutTrim <- data.frame(sample_data(subset_samples(phyT.mdGut, AnimalSource%in%c("Wlot","GrotonMA")))) # separate sample data from phyloseq object (sample data . macrobdella decora Gut Trim)
pairwise.adonis(bray.mdGutTrim,paste(sd.mdGutTrim$Sample_Type, sd.mdGutTrim$AnimalSource)) # pairwise Adonis calculation comparing combination of Sample_Type and AnimalSoure for Macrobdella decora gut samples from CT and MA

#####################################################################
########## FIGURE 5 #################################################
#####################################################################
### ILF Mdecora ###
#phyT.mdILF <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples (physeq Transform . macrobdella decora ILF)
phyT.mdILFt <- subset_samples(phyT.mdILF, Da2F%in%c("none","Unk","0") & AnimalSource%in%c("Wlot","GrotonMA")) # subset containing only Macrobdella decora ILF samples from CT and MA (physeq Transform . macrobdella decora ILF trimmed)
phyT.mdILFp <- subset_taxa(phyT.mdILFt, ASV%in%c(ls.coreGut)) # keep only common ILF and intestinum taxa (phyloseq Transformed . macrobdella decora ILF pruned)
ls.daysMdILF <- unique(sample_data(phyT.mdILFp)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days Macrobdella decora ILF)
### Loop to remove ASVs at low prevalence by days after feeding ###
for (d in c(ls.daysMdILF)){
  phyT.AAa<-subset_samples(phyT.mdILFp,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.mdILFp,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= cMin)), by = TaxaID] # make simple table listing 'ASV and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of ASVs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[,eval(f)] <- .00001 # change all prevalence values to ~0 for ASVs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all ASVs in ls.low
  phyT.mdILFp <- merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . macrobdella decora Box+Whisker ILF)
} # repeat this loop for all days after feeding (d in ls.daysMdILF) 

### ILF Hverbana ###
#phyT.hvILF <- subset_samples(phyT.leech, Taxonomic_ID=="Hverbana" & Sample_Type=="ILF") # subset containing only Hirudo verbana ILF samples (physeq Transform . hirudo verbana ILF)
phyT.hvILFt <- subset_samples(phyT.hvILF, Da2F%in%c("none","Unk","0")) # subset containing only Hirudo verbana ILF samples (physeq Transform . hirudo verbana ILF trimmed)
phyT.hvILFp <- subset_taxa(phyT.hvILFt, ASV%in%c(ls.coreGut)) # keep only common ILF and intestinum taxa (phyloseq Transformed . hirudo verbana ILF pruned)
ls.daysHvILF <- unique(sample_data(phyT.hvILFp)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days Hirudo verbana ILF)
### Loop to remove ASVs at low prevalence by days after feeding ###
for (d in c(ls.daysHvILF)){
  phyT.AAa<-subset_samples(phyT.hvILFp,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.hvILFp,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= cMin)), by = TaxaID] # make simple table listing 'ASV and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of ASVs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[,eval(f)] <- .00001 # change all prevalence values to ~0 for ASVs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all ASVs in ls.low
  phyT.hvILFp <- merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . hirudo verbana Box+Whisker ILF)
} # repeat this loop for all days after feeding (d in ls.daysHvILF) 

### Intestinum M. decora ###
#phyT.mdInt <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="Intestinum") # subset containing only Macrobdella decora intestinum samples (physeq Transform . macrobdella decora Intestinum)
phyT.mdIntt <- subset_samples(phyT.mdInt , Da2F%in%c("none","Unk","0") & AnimalSource%in%c("Wlot","GrotonMA")) # subset containing only Macrobdella decora intestinum samples from CT and MA (physeq Transform . macrobdella decora Intestinum trimmed)
phyT.mdIntp <- subset_taxa(phyT.mdIntt, ASV%in%c(ls.coreGut)) # keep only common ILF and intestinum taxa (phyloseq Transformed . macrobdella decora intestinum pruned)
ls.daysMdInt <- unique(sample_data(phyT.mdIntp)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days Macrobdella decora Intestinum)
### Loop to remove ASVs at low prevalence by days after feeding ###
for (d in c(ls.daysMdInt)){
  phyT.AAa<-subset_samples(phyT.mdIntp,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.mdIntp,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= cMin)), by = TaxaID] # make simple table listing 'ASV and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of ASVs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[,eval(f)] <- .00001 # change all prevalence values to ~0 for ASVs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all ASVs in ls.low
  phyT.mdIntp <- merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . macrobdella decora Box+Whisker intestinum)
} # repeat this loop for all days after feeding (d in ls.daysMdInt) 

phyT.gutMerge <- merge_phyloseq(phyT.mdILFp,phyT.hvILFp,phyT.mdIntp) # merge the modified data from Mdecora ILF, Hverbana ILF, and Mdecora intestinum (phyloseq Transformed . gut Merged)
### 'spp2' column will be used to identify leech species for labeling ###
sample_data(phyT.gutMerge)$spp2 <- with(sample_data(phyT.gutMerge),
                                          ifelse(Taxonomic_ID=="Hverbana","H. verbana",
                                          ifelse(Taxonomic_ID=="Mdecora","M. decora",
                                          as.character(Taxonomic_ID))))
sample_data(phyT.gutMerge)$Da1Fb = factor(sample_data(phyT.gutMerge)$Da1Fb, levels = c(0,1,2,4,7,"30+","90+")) # Reorder Da1Fb
sample_data(phyT.gutMerge)$Taxonomic_ID = factor(sample_data(phyT.gutMerge)$Taxonomic_ID, levels = c("Hverbana","Mdecora")) # Reorder Taxonomic_ID
phyT.gutMergep <- subset_taxa(phyT.gutMerge,taxa_sums(phyT.gutMerge) > 0.1) # (phyloseq Transformed . gut Merged pruned)
dt.gutMerge <- psmelt(phyT.gutMerge) # create data.table from phyloseq object (data table . Box+Whisker merged)

pBW.DaF <- ggplot(dt.gutMerge) +
  stat_boxplot(data=dt.gutMerge,aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  theme_bw() +
  scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(10^-3),expression(10^-2),expression(10^-1),""),limits=c(.001,1)) +
  xlab("Days After Feeding") +
  theme(text=element_text(size=10),strip.text.y = element_text(face = "italic",size=7),axis.text.x=element_text(angle=90,hjust=1),legend.position="none",panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=brewer.pal(8,"Greys")) 
pBW.DaF<-pBW.DaF + facet_grid(Phylum+GenusASV~spp2+Sample_Type) # facet plot by Phylum/Genus (y-axis) and leech species/sample type (x-axis)
pBW.DaF # print plot
ggsave(plot=pBW.DaF, filename="Plots/Mdfig5_plotBW.pdf", device="pdf",dpi="retina",width=5,height=9.06,units="in") # Save FIGURE 5 to .pdf file
########## end FIGURE 5 #############################################

#####################################################################
########## FIGURE 6 #################################################
#####################################################################
##### Bray-Curtis distance by Da1F #####
# Thank you to jeffkimbrel ! #
#phyT.mdILF <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples (physeq Transform . macrobdella decora ILF)
sample_data(phyT.mdILF)$name <- sample_names(phyT.mdILF) # add column 'name' to sample_data
dist.mdILF <- phyloseq::distance(phyT.mdILF, method = "bray") # calculate bray curtis distance matrix (bray . macrobdella decora ILF)
# remove self-comparisons
md.mdILF = reshape2::melt(as.matrix(dist.mdILF)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)
# get sample data (S4 error OK and expected)
sd.mdILF = sample_data(phyT.mdILF) %>%
  select("name","Da1Fb","Header","MealLot1","Date1stFeed","WildMonth","DNA_extractor","Season") %>%
  mutate_if(is.factor,as.character) 
# combined distances with sample data
colnames(sd.mdILF) = c("Var1", "Type1","Source","MealLot","Date1stFeed","Month","Initials","Season")
md.mdILF2 = left_join(md.mdILF, sd.mdILF, by = "Var1")
colnames(sd.mdILF) = c("Var2", "Da1F","Source","MealLot","Date1stFeed","WildMonth","Initials","Season")
md.mdILF3 = left_join(md.mdILF2, sd.mdILF, by = "Var2")

wu.sd2<-md.mdILF3[md.mdILF3$Type1==0,]
wu.sd2$Da1F = factor(wu.sd2$Da1F, levels = c(0,1,2,4,7,30,"30+","90+")) # Reorder Da1Fb
colnames(wu.sd2)[colnames(wu.sd2)=="value"] <- "distance"

pBray.mdILFns <- ggplot(wu.sd2, aes(x = Da1F, y = distance,color=Season.y)) +
  theme_bw() +
  geom_violin(aes(x = Da1F, y = distance)) +
  geom_jitter(width = 0.2) +
  ylab("Bray-Curtis Distance") +
  xlab("Days After Feeding") +
  scale_color_manual(values = pal.CB) +
  theme(text=element_text(family="Times New Roman", size=12),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),legend.title = element_blank())
pBray.mdILFns # uncomment to print plot

ggsave(pBray.mdILFns,filename="Plots/Mdfig6_plotMdBrayFeed.pdf", device=cairo_pdf, dpi="retina",width=5,height=4,units="in") # save figure to file
########## end FIGURE 6 #############################################

#####################################################################
########## Statistics for Figure 6 Comparisons ######################
########## TABLE S4 #################################################
#####################################################################
rel.mdDaFilf <- microbiome::transform(phyT.mdDaFilf, "compositional") 
ASV.mdDaFilf <- abundances(rel.mdDaFilf)
meta.mdDaFilf <- meta(rel.mdDaFilf)
dist.mdDaFilf <- vegdist(t(ASV.mdDaFilf))
perm.mdDaFilf <- permutest(betadisper(dist.mdDaFilf, meta.mdDaFilf$Da1Fb), pairwise = TRUE)

write.table(perm.mdDaFilf$pairwise$permuted, "tableSig_DaFMd_ILF_tblS4data.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included
########## end TABLE S4 ############################################

####################################################################
########## SUPPLEMENTAL FIGURE 1 ###################################
####################################################################
### ILF NMDS - Collection Site ###
#phyT.mdILF <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples (physeq Transform . macrobdella decora ILF)
dsNMDS.mdILF <- phyT.mdILF # define data for analysis
sample_data(dsNMDS.mdILF)$AnimalSource = factor(sample_data(dsNMDS.mdILF)$AnimalSource, levels = c("Wlot","GrotonMA","CarogaNY","MtSnowVT")) # force AnimalSource order (instead of alphabetical default)
mNMDS <-"bray" # define metric for analysis (metric NMDS)
dist.mdILF = phyloseq::distance(dsNMDS.mdILF, method = c(mNMDS)) # calculate distances (distance . macrobdella decora ILF)
ord.mdILF = ordinate(dsNMDS.mdILF, method = "NMDS", distance = dist.mdILF) # calculate ordination (ordination . macrobdella decora ILF)
nmds.mdILFstate<-plot_ordination(dsNMDS.mdILF, ord.mdILF, color = "AnimalSource") + 
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#000000")) +
  theme_bw() +
  theme(text=element_text(size=12), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
nmds.mdILFstate$layers[[1]]<-NULL # remove layer with geom_points
nmds.mdILFstate<-nmds.mdILFstate + 
  geom_point(aes(fill=AnimalSource),size=3,alpha=0.7) +  
  stat_ellipse(type = "t", level = 0.95, linetype = 2) # re-write geom_point layer and add elipses
#nmds.mdILFstate # print plot
### Intestinum/WildTime NMDS ###
#phyT.mdInt <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="Intestinum") # subset containing only Macrobdella decora intestinum samples (physeq Transform . macrobdella decora Intestinum)
dsNMDS.mdInt <- phyT.mdInt # define data for analysis (data set NMDS . macrobdella decora Intestinum)
sample_data(dsNMDS.mdInt)$AnimalSource = factor(sample_data(dsNMDS.mdInt)$AnimalSource, levels = c("Wlot","GrotonMA","CarogaNY","MtSnowVT")) # force AnimalSource order (instead of alphabetical default)
dist.mdInt = phyloseq::distance(dsNMDS.mdInt, method = c(mNMDS)) # calculate distances (distance . macrobdella decora Intestinum)
ord.mdInt = ordinate(dsNMDS.mdInt, method = "NMDS", distance = dist.mdInt) # calculate ordination (ordination . macrobdella decora Intestinum)
nmds.mdIntState<-plot_ordination(dsNMDS.mdInt, ord.mdInt, color = "AnimalSource") + 
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#000000")) +
  theme_bw() +
  theme(text=element_text(size=12), axis.title.x=element_blank(),axis.title.y=element_blank(),legend.position="none")
nmds.mdIntState$layers[[1]]<-NULL # remove layer with geom_points
nmds.mdIntState<-nmds.mdIntState + 
  geom_point(aes(fill=AnimalSource),size=3,alpha=0.7) +  
  stat_ellipse(type = "t", level = 0.95, linetype = 2) # re-write geom_point layer and add elipses
#nmds.mdIntState # print plot
ggsave(plot_grid(nmds.mdILFstate, nmds.mdIntState, labels = "AUTO"), filename="NMDS/MdfigS1_plotNMDS.pdf", device=cairo_pdf, dpi="retina",width=6.87,units="in") # Save supplemental FIGURE 1 to file
########## end SUPPLEMENTAL FIGURE 1 ################################

#####################################################################
########## TABLE 4 ##################################################
#####################################################################
phyT.adultB <- merge_phyloseq(phyT.mdBlad,phyT.hvBlad) # create merged phyloseq object with Macrobdella decora and Hirudo verbana bladder samples (physeq Transformed . adult Bladder)
ls.comBlad <- unique(c(ls.comHvBlad,ls.comMdBlad,ls.comMaBlad,ls.comCtBlad)) # compile list of common ASVs for bladder samples (list . common Bladder)
df.adultB <- psmelt(subset_taxa(phyT.adultB,ASV%in%ls.comBlad)) # create data frame from phyloseq object with only ASVs listed (data frame . adult Bladder)

df.Bcounts<-data.frame(Sequence=character(),
                       HvBladabund=numeric(),
                       HvBladmax=numeric(),
                       HvBladstDev=numeric(),
                       HvBladprev=numeric(),
                       MdBladabund=numeric(),
                       MdBladmax=numeric(),
                       MdBladstDev=numeric(),
                       MdBladprev=numeric(),
                       stringsAsFactors=FALSE) # create empty data frame for counting bladder data (data frame . Bladder counts)
### Loop for calculating values to be entered into df.Bcounts ###
for (ASV in unique(df.adultB$ASV)){
  HvBladabund <- mean(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder abundance)
  HvBladmax <- max(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana"],na.rm=TRUE) # calculate maximum abundance for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder max)
  HvBladstDev <- sd(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana"]) # calculate standard deviation for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder standard Deviation)
  HvBladprev <- nrow(df.adultB[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana" & df.adultB$Abundance > cMin,])/nsamples(phyT.hvBlad) # count prevalence for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder prevlance)
  
  MdBladabund <- mean(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora"],na.rm=TRUE) # calculate mean abundance for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder abundance)
  MdBladmax <- max(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora"],na.rm=TRUE) # calculate maximum abundance for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder max)
  MdBladstDev <- sd(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora"]) # calculate standard deviation for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder standard Deviation)
  MdBladprev <- nrow(df.adultB[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora" & df.adultB$Abundance > cMin,])/nsamples(phyT.mdBlad) # count prevalence for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder prevlance)
  
  df.Bcounts[nrow(df.Bcounts)+1,] = list(ASV,
                                       HvBladabund,HvBladmax,HvBladstDev,HvBladprev,
                                       MdBladabund,MdBladmax,MdBladstDev,MdBladprev) # for each ASV, record calculated values in data frame (data frame . Bladder counts)
}

colnames(df.Bcounts)[1] <- "ASV" # rename first column (Sequence) as 'ASV'
df.Bsum <- merge(x=df.Bcounts,y=tax_table(phyT.adultB),by="ASV") # combine taxonomy table with calculated values for ASVs (data frame . Bladder sum)
df.Bfinal <- df.Bsum[,c("ASV","Phylum","Class","Genus",
                      "HvBladabund","HvBladmax","HvBladstDev","HvBladprev",
                      "MdBladabund","MdBladmax","MdBladstDev","MdBladprev")] # choose limited headings to include in data frame (data frame . Bladder final)
write.csv(df.Bfinal,"table_BladCounts_tbl4data.csv", row.names=FALSE) # write data frame to csv without row names
########## end TABLE 4 ##############################################

########## TABLE 2 ##################################################
#phyT.hvGut <- subset_samples(phyT.hv, Sample_Type%in%c("ILF","Intestinum")) # subset containing only Hirudo verbana ILF and intestinum samples (physeq Transform . hirudo verbana Gut)
#phyT.mdGut <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type%in%c("ILF","Intestinum")) # subset containing only Macrobdella decora ILF and intestinum samples (physeq Transform . macrobdella decora Gut)
phyT.gut <- merge_phyloseq(phyT.hvGut,phyT.mdGut) # combine phyloseq objects (physeq Transform . gut)
phyT.gutp <- subset_taxa(phyT.gut, taxa_sums(phyT.gut) > cMin) # keep only taxa present at minimum value (physeq Transform . gut pruned)
df.gut <- psmelt(phyT.gutp) # create data frame from phyloseq object (data frame . gut)

df.totCounts<-data.frame(Sequence=character(),
                         HvILFabund=numeric(),
                         HvILFabundBBEZ=numeric(),
                         HvILFabundUSA=numeric(),
                         HvILFstDev=numeric(),
                         HvILFprev=numeric(),
                         HvIntabund=numeric(),
                         HvIntstDev=numeric(),
                         HvIntprev=numeric(),
                         HvBladabund=numeric(),
                         HvBladstDev=numeric(),
                         HvBladprev=numeric(),
                         MdILFabund=numeric(),
                         MdILFstDev=numeric(),
                         MdILFprev=numeric(),
                         MdIntabund=numeric(),
                         MdIntstDev=numeric(),
                         MdIntprev=numeric(),
                         MdBladabund=numeric(),
                         MdBladstDev=numeric(),
                         MdBladprev=numeric(),
                         MaxAbund=numeric(),
                         stringsAsFactors=FALSE) # create empty data frame for counting data (data frame . total Counts)
### Loop for calculating values to be entered into df.totCounts ###
for (ASV in unique(df.gut$ASV)){
  HvILFabund <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="ILF"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana ILF samples (Hirudo verbana ILF abundance)
  HvILFabundBBEZ <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="ILF" & df.gut$AnimalSource=="BBEZ"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana ILF samples from supplier 1 (Hirudo verbana ILF abundance BBEZ)
  HvILFabundUSA <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="ILF" & df.gut$AnimalSource=="USA"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana ILF samples from supplier 2 (Hirudo verbana ILF abundance USA)
  HvILFstDev <- sd(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="ILF"]) # calculate standard deviation for ASV in Hirudo verbana ILF samples (Hirudo verbana ILF standard Deviation)
  HvILFprev <- nrow(df.gut[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="ILF" & df.gut$Abundance > cMin,])/nsamples(subset_samples(phyT.gut,Taxonomic_ID=="Hverbana" & Sample_Type=="ILF")) # count prevalence for ASV in Hirudo verbana ILF samples (Hirudo verbana ILF prevlance)
  HvIntabund <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="Intestinum"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana intestinum samples (Hirudo verbana Intestinum abundance)
  HvIntstDev <- sd(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="Intestinum"]) # calculate standard deviation for ASV in Hirudo verbana intestinum samples (Hirudo verbana Intestinum standard Deviation)
  HvIntprev <- nrow(df.gut[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Hverbana" & df.gut$Sample_Type=="Intestinum" & df.gut$Abundance > cMin,])/nsamples(subset_samples(phyT.gut,Taxonomic_ID=="Hverbana" & Sample_Type=="Intestinum")) # count prevalence for ASV in Hirudo verbana intestinum samples (Hirudo verbana Intestinum prevlance)
  HvBladabund <- mean(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana"],na.rm=TRUE) # calculate mean abundance for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder abundance)
  HvBladstDev <- sd(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana"]) # calculate standard deviation for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder standard Deviation)
  HvBladprev <- nrow(df.adultB[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Hverbana" & df.adultB$Abundance > cMin,])/nsamples(phyT.hvBlad) # count prevalence for ASV in Hirudo verbana bladder samples (Hirudo verbana Bladder prevlance)
  
  MdILFabund <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="ILF"],na.rm=TRUE) # calculate mean abundance for ASV in Macrobdella decora ILF samples (Macrobdella decora ILF abundance)
  MdILFstDev <- sd(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="ILF"]) # calculate standard deviation for ASV in Macrobdella decora ILF samples (Macrobdella decora ILF standard Deviation)
  MdILFprev <- nrow(df.gut[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="ILF" & df.gut$Abundance > cMin,])/nsamples(subset_samples(phyT.gut,Taxonomic_ID=="Mdecora" & Sample_Type=="ILF")) # count prevalence for ASV in Macrobdella decora ILF samples (Macrobdella decora ILF prevlance)
  MdIntabund <- mean(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="Intestinum"],na.rm=TRUE) # calculate mean abundance for ASV in Macrobdella decora intestinum samples (Macrobdella decora Intestinum abundance)
  MdIntstDev <- sd(df.gut$Abundance[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="Intestinum"]) # calculate standard deviation for ASV in Macrobdella decora intestinum samples (Macrobdella decora Intestinum standard Deviation)
  MdIntprev <- nrow(df.gut[df.gut$ASV==eval(ASV) & df.gut$Taxonomic_ID=="Mdecora" & df.gut$Sample_Type=="Intestinum" & df.gut$Abundance > cMin,])/nsamples(subset_samples(phyT.gut,Taxonomic_ID=="Mdecora" & Sample_Type=="Intestinum")) # count prevalence for ASV in Macrobdella decora intestinum samples (Macrobdella decora Intestinum prevlance)
  MdBladabund <- mean(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora"],na.rm=TRUE) # calculate mean abundance for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder abundance)
  MdBladstDev <- sd(df.adultB$Abundance[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora"]) # calculate standard deviation for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder standard Deviation)
  MdBladprev <- nrow(df.adultB[df.adultB$ASV==eval(ASV) & df.adultB$Taxonomic_ID=="Mdecora" & df.adultB$Abundance > cMin,])/nsamples(phyT.mdBlad) # count prevalence for ASV in Macrobdella decora bladder samples (Macrobdella decora Bladder prevlance)
  
  MaxAbund <- max(df.gut$Abundance[df.gut$ASV==eval(ASV)],na.rm=TRUE) # calculate maximum abundance for ASV in all samples (Maxumum Abundance)
  
  df.totCounts[nrow(df.totCounts)+1,] = list(ASV,
                                           HvILFabund,HvILFabundBBEZ,HvILFabundUSA,HvILFstDev,HvILFprev,
                                           HvIntabund,HvIntstDev,HvIntprev,
                                           HvBladabund,HvBladstDev,HvBladprev,
                                           MdILFabund,MdILFstDev,MdILFprev,
                                           MdIntabund,MdIntstDev,MdIntprev,
                                           MdBladabund,MdBladstDev,MdBladprev,
                                           MaxAbund) # for each ASV, record calculated values in data frame (data frame . total Counts)
}

colnames(df.totCounts)[1] <- "ASV" # rename first column (Sequence) as 'ASV'
df.totSum <- merge(x=df.totCounts,y=tax_table(phyT.gut),by="ASV") # combine taxonomy table with calculated values for ASVs (data frame . total Sum)
df.totFinal <- df.totSum[,c("Phylum","Class","Order","Family","Genus","GenusASV",
                          "HvILFabund","HvILFprev",
                          "HvIntabund","HvIntprev",
                          "HvBladabund","HvBladprev",
                          "MdILFabund","MdILFprev",
                          "MdIntabund","MdIntprev",
                          "MdBladabund","MdBladprev",
                          "HvILFstDev","HvIntstDev","HvBladstDev",
                          "MdILFstDev","MdIntstDev","MdBladstDev",
                          "ASV","Sequence","MaxAbund")] # choose limited headings to include in data frame (data frame . total Final)
write.csv(df.totFinal,"table_ILFIntCounts_tblS3data.csv", row.names=FALSE) # write data frame to csv without row names

df.totCommon <- df.totSum[,c("Phylum","Class","Order","Family","Genus","GenusASV",
                           "HvILFabundBBEZ","HvILFabundUSA","HvIntabund",
                           "MdILFabund","MdIntabund",
                           "ASV","Sequence")] # choose limited headings to include in data frame (data frame . total Common)
df.totComASV <- df.totCommon[which(df.totCommon$ASV %in%c(ls.comTot)), ] # choose limited ASVs to include in data frame (data frame . total Common ASV)
write.csv(df.totComASV,"table_ILFIntCounts_tbl2data.csv", row.names=FALSE) # write data frame to csv without row names
########## end TABLE 2 ##############################################

#####################################################################
########## FIGURE S2 ################################################
#####################################################################
########## Heat map of significant ASVS in Da1Fb groups #############
#phyT.mdILF <- subset_samples(phyT.leech, Taxonomic_ID=="Mdecora" & Sample_Type=="ILF") # subset containing only Macrobdella decora ILF samples (physeq Transform . macrobdella decora ILF)
phyR.mdILF <- subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%sample_names(phyT.mdILF) ) # create data set with raw reads for Macrobdella decora ILF samples
phyG.mdILF <- subset_taxa(phyR.mdILF, taxa_sums(phyR.mdILF) > var.medNeg) # remove any taxa that are not present to reduce calculation time (physeq G . macrobdella decora ILF)
fmG.mdILF = fast_melt(phyG.mdILF) # make data table from physeq object (fast melt . macrobdella decora ILF)
prevG.mdILF = fmG.mdILF[, list(Prevalence = sum(count > 1)),by = ASV] # make simple table listing 'ASV and Prevalencet' (prevalence G . macrobdella decora ILF)
phyG.mdILFp <- subset_taxa(phyG.mdILF,ASV%in%c(prevG.mdILF[(Prevalence > 1), ASV])) # remove taxa not present in more than 1 sample (physeq G . macrobdella decora ILF pruned)
taxa_names(phyG.mdILFp) <- tax_table(phyG.mdILFp)[,"GenusASV"] # identify ASVs by genus instead of sequence

sample_data(phyG.mdILFp) <- as.factor(sample_data(phyG.mdILFp)) # define sample_data as factors...because R 

dsq.mdILF = phyloseq_to_deseq2(phyG.mdILFp, ~ Da1Fb) # change phyloseq data to DESeq2 format contrasting Da1Fb (DeSeq . macrobdella decora ILF)
#dsq.mdILF$Da1Fb <- factor(dsq.mdILF$Da1Fb, levels = c("0","1","4","7","30+","90+")) # force days after feeeding order
geo.mdILF = apply(counts(dsq.mdILF), 1, gm_mean)
esf.mdILF = estimateSizeFactors(dsq.mdILF, geoMeans = geo.mdILF) # estimate size factors for DeSeq object (estimated size factor . macrobdella decora ILF)
dsq.mdILF = DESeq(esf.mdILF) 

res.mdILF = results(dsq.mdILF, contrast=c("Da1Fb","0","4")) # input the contrast variable you will be using, followed by the two values to contrast (results . macrobdella decora ILF)
res.mdILF = res.mdILF[order(res.mdILF$baseMean, na.last=NA), ] 
st.mdILF = res.mdILF[which(res.mdILF$pvalue < .5), ] # keep results with low p values (significance table . macrobdella decora ILF)
ls.keep <- row.names(st.mdILF) # list ASVs from significance table (list . keep)

nt <- normTransform(dsq.mdILF) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[ls.keep, ]
df <- as.data.frame(colData(dsq.mdILF)[,c("Da1Fb","AnimalSource")])

ann_colors = list(
  Da1Fb = c("0"="#ffffff","1"="#f2ceb2","2"="#e18e4c","4"="#D55e00","7"="#954100","30+"="#999999","90+"="#4c4c4c"),
  AnimalSource = c("Wlot"="#56b4e9","GrotonMA"="#009e73","MtSnowVT"="#e69f00","CarogaNY"="#f0e442")
)

pHeat.Da1Fb <- pheatmap::pheatmap(log2.norm.counts, annotation_col=df, annotation_colors = ann_colors,show_colnames = F,fontface_row="italic",fontsize = 8)
ggsave(pHeat.Da1Fb, filename="Plots/MdFigS2_plotHeat.pdf", dpi="retina",width=7.5,height=5,units="in") # Save FIGURE S2 to file
########## end FIGURE S2 ##################################################

#####################################################################
########## Table S5 #################################################
#####################################################################
phyR.mdILF <- subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%sample_names(phyT.mdILF) ) # create data set with raw reads for Macrobdella decora ILF samples
phyG.mdILF <- subset_taxa(phyR.mdILF, taxa_sums(phyR.mdILF) > var.medNeg) # remove any taxa that are not present to reduce calculation time (physeq G . macrobdella decora ILF)
taxa_names(phyG.mdILF) <- tax_table(phyG.mdILF)[,"ASV"] # identify ASVs by ASV
ls.Da1F <- c("1","4","7","30+","90+") # (list . Days after Feeding)

fm.mdASV = fast_melt(phyG.mdILF) # convert ASV table into a long format with three main columns: SampleID, TaxaID, and count (fast melt . hirudo verbana Gut)
prevdt.mdASV = fm.mdASV[, list(Prevalence = sum(count >= var.medNeg)), by = ASV] # make simple table listing ASV and Prevalence (prevalence data table . Hirudo verbana Gut)
ls.mdASV = prevdt.mdASV[(Prevalence > 2), ASV] # make list of ASVs present in more than 1 sample (list . keep Hirudo verbana Gut)
phyG.mdILFt <- subset_taxa(phyG.mdILF,ASV%in%ls.mdASV)
phyG.mdILFt <- subset_samples(phyG.mdILF, Da1Fb!="2")

dsq.mdASV = phyloseq_to_deseq2(phyG.mdILFt, ~ Da1Fb) # change phyloseq data to DESeq2 format. AFter ~ indicates variables to contrast
dsq.mdASV$Da1Fb <- factor(dsq.mdASV$Da1Fb, levels = c("0","1","4","7","30+","90+"))
geo.mdASV = apply(counts(dsq.mdASV), 1, gm_mean)
esf.mdASV = estimateSizeFactors(dsq.mdASV, geoMeans = geo.mdASV)
dsq.mdASV = DESeq(esf.mdASV, fitType="local")

mtx.mdASV <- matrix(ncol=length(ls.Da1F), nrow=length(ls.mdASV)) # create empty matrix (matrix . significant Day after Feeding by asv)
dimnames(mtx.mdASV) = list(c(ls.mdASV),c(ls.Da1F)) # label rows as ASV names and columns as days after feeding

for(day in c(ls.Da1F)){
  loop.res = results(dsq.mdASV, contrast=c("Da1Fb","0",day)) # input the contrast variable you will be using, followed by the two values to contrast
  for(asv in c(ls.mdASV)){
    mtx.mdASV[asv,day] <- loop.res[asv,"padj"] # copy padj value to proper [ASV,day] in matrix
    }
}

write.csv(mtx.mdASV,"table_mdDaFasv_tblS5data.csv") # output matrix of padj values to csv
########## end Table S5 #############################################

dsq.start = phyloseq_to_deseq2(phyR.hvGutTrim, ~ Sample_Type) 
geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.1
sigtab = res[(res$padj < alpha & res$baseMean > var.medNeg & res$padj < 0.01), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.hvGutTrim)[rownames(sigtab), ], "matrix"))
head(sigtab)

posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 1, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","ASV")]
ls.02<-levels(posigtab$ASV)
sigtabgen = sigtab

#####################################################################
########## Figure S3 ################################################
#####################################################################
phyR.mdILF <- subset_samples(phyR.sin,Taxonomic_ID=="Mdecora" & Sample_Type=="ILF" & AnimalSource%in%c("Wlot","GrotonMA"))
sample_data(phyR.mdILF)$Da1Fb = factor(sample_data(phyR.mdILF)$Da1Fb, levels = c(0,1,2,4,7,"30+","90+")) # Reorder Da1Fb

pDiv.adMdILF <- plot_richness(phyR.mdILF, x = "Da1Fb", color="Header",measures=c("Shannon")) +
  geom_boxplot(aes()) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,2.2)) +
  scale_color_manual(values=pal.CB) +
  ylab("Shannon Index") +
  xlab("Days After Feeding") +
  theme_bw() +
  theme(text=element_text(family="Times", size=12),axis.text.x=element_text(angle = 0),strip.text.x=element_blank(),legend.title = element_blank()) 
pDiv.adMdILF  # uncomment to print plot 

ggsave(pDiv.adMdILF,filename="Plots/MdfigS3_plotMdAlphaFeed.pdf",device=cairo_pdf, dpi="retina",width=6,height=5,units="in")
########## End figure S3 ############################################




#####################################################################
# END # END # END # END # END # END # END # END # END # END # END #

sample_data(phyT.mdInt)$name <- sample_names(phyT.mdInt) # add column 'name' to sample_data
dist.mdInt <- phyloseq::distance(phyT.mdInt, method = "bray") # calculate bray curtis distance matrix (bray . macrobdella decora Int)
# remove self-comparisons
md.mdInt = reshape2::melt(as.matrix(dist.mdInt)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)
# get sample data (S4 error OK and expected)
sd.mdInt = sample_data(phyT.mdInt) %>%
  select("name","Da1Fb","Header","MealLot1","Date1stFeed","WildMonth","DNA_extractor","Season") %>%
  mutate_if(is.factor,as.character) 
# combined distances with sample data
colnames(sd.mdInt) = c("Var1", "Type1","Source","MealLot","Date1stFeed","Month","Initials","Season")
md.mdInt2 = left_join(md.mdInt, sd.mdInt, by = "Var1")
colnames(sd.mdInt) = c("Var2", "Da1F","Source","MealLot","Date1stFeed","WildMonth","Initials","Season")
md.mdInt3 = left_join(md.mdInt2, sd.mdInt, by = "Var2")

wu.sd2<-md.mdInt3[md.mdInt3$Type1==0,]
wu.sd2$Da1F = factor(wu.sd2$Da1F, levels = c(0,1,2,4,7,30,"30+","90+")) # Reorder Da1Fb
colnames(wu.sd2)[colnames(wu.sd2)=="value"] <- "distance"

pBray.mdIntns <- ggplot(wu.sd2, aes(x = Da1F, y = distance,color=Season.y)) +
  theme_bw() +
  geom_violin(aes(x = Da1F, y = distance)) +
  geom_jitter(width = 0.2) +
  ylab("Bray-Curtis Distance") +
  xlab("Days After Feeding") +
  scale_color_manual(values = pal.CB) +
  theme(text=element_text(family="Times New Roman", size=12),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),legend.title = element_blank())
pBray.mdIntns # uncomment to print plot

