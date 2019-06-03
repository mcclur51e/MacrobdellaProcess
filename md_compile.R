#######################################
########## Processing Set-up ##########
#######################################
########## Call libraries for use ##########
library("phyloseq") #packageVersion("phyloseq")
library("ggplot2") #packageVersion("ggplot2")
library("ape") #packageVersion("ape")
library("data.table") #packageVersion("data.table")
library("RColorBrewer") # design new color palette for data
library("grid")
library("cowplot") # used to prepare final figures
library("vegan")
library("decontam") # identify contaminant OTUs 
library("dplyr")
library("pairwiseAdonis") # for PERMANOVA calculations
library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. #
library("structSSI")

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
pal.rainbow<-c("#ff00aa","#ac00e6","#3333ff","#0ba29a","#39ac39","#ffff00","#ff9933","#F00000","#800055")
pal.CB<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000") # colorblind-friendly color palette
########## Set working directory where processing will be carried out. Make subdirectories for saving into ##########
setwd("~/Dropbox/Manuscript_Macrobdella/Md_processData/") # set working directory for data processing. Must contain folder "DataFiles" containing "table_otu.csv", "table_tax.csv", and "table_map.csv" and/or "physeq_manuscript.RData"
dir.create("NegControls",showWarnings = FALSE) # create folder for negative control data output(s)
dir.create("NMDS",showWarnings = FALSE) # create folder for NMDS data output(s)
dir.create("Plots",showWarnings = FALSE) # create folder for plot output(s)
########## Input basic files to create phyloseq object ##########
OTU = otu_table(as.matrix.data.frame(read.csv("DataFiles/table_otu.csv", header = TRUE, row.names = 1, check.names=FALSE) ), taxa_are_rows = TRUE) # reads data.frame into matrix without converting to numbers # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with OTU data for each sample
TAX = tax_table(as.matrix.data.frame(read.csv("DataFiles/table_tax.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads data.frame into matrix without converting to numbers # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with taxa data separated by columns and headers added for Kingdom, Phylum, Class, etc. First column = otu IDs
MAP = sample_data(data.frame(read.csv("DataFiles/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads csv file into data.frame with row names in column 1
TREE = read.tree("DataFiles/Rep_Set_tree.tree")
physeq = phyloseq(OTU, TAX, MAP,TREE)
save(physeq,file=("DataFiles/physeq_current.RData")) # Save the phyloseq data object in a .RData file 
########## Keep only samples for publishing ##########
ls.pub<-as.character(read.csv("DataFiles/list_samplePublish.csv",header=FALSE,sep=",")$V1) # Import list of samples 
phyR.McC<-subset_samples(physeq, sample_names(physeq)%in%c(ls.pub))# Keep only samples for publishing (physeq.McClure)
save(phyR.McC,file=("DataFiles/physeq_publish.RData")) # Save phyR.McC in a .RData file 
phyR.trim1<-subset_taxa(phyR.McC,taxa_sums(phyR.McC)>1) # keep only OTUs with >1 read total
########## Load most recent data ##########
### Use this to skip lines 35--48 if the basic phyloseq data sets have previously been produced and saved ###
# load("DataFiles/physeq_current.RData") # Use to re-load saved phyloseq data object 
# load("DataFiles/physeq_publish.RData") # Use to re-load saved phyloseq data object

##################################################################
########## Pre-process data to remove contaminant OTUs  ##########
##################################################################
phyR.PCRneg<-subset_samples(phyR.trim1,Sample_Type%in%c("PCRneg","Reagents")) # keep PCR negative control samples
phyR.taxaNeg<-subset_taxa(phyR.PCRneg,taxa_sums(phyR.PCRneg)>1) # keep only OTUs with >1 read total
var.maxNeg <- max(otu_table(phyR.taxaNeg)) # find the greatest count of any OTU in a single sample

### Trim unuseable samples ###
phyR.negConc<-subset_samples(phyR.trim1,!PostPCRDNA_ng_uL%in%c("Unk")) # Remove any samples where Post PCR [DNA] was not calculated
phyR.negCount<-subset_samples(phyR.negConc,sample_sums(phyR.negConc)>=1) # keep samples with >=1 reads
phyR.1min<-prune_taxa(taxa_sums(phyR.negCount)>1,phyR.negCount) # keep only taxa with >1 reads

### Map ###
map.neg<-sample_data(phyR.1min) # pull map data from phyloseq object
map.neg$Conc<-with(map.neg,
                     ifelse(PostPCRDNA_ng_uL%in%c("zero","Unk",""), "0", 
                            as.character(PostPCRDNA_ng_uL)))
map.neg$c <- as.numeric(as.character(map.neg$Conc)) * 10 + 1 # ridiculous math to make decontam work
map.neg$q2 <- as.numeric(as.character(map.neg$c))
phyR.1min = merge_phyloseq(phyR.1min,map.neg) # return map data to the phyloseq object

### FREQUENCY ###
contamdf.freq <- isContaminant(phyR.1min, method="frequency", conc="q2")
phyR.contamFreq <- prune_taxa(contamdf.freq$contaminant, phyR.1min)

phyR.allc <- prune_taxa(taxa_names(phyR.contamFreq), phyR.1min) # (physeq all contaminants)
phyR.nc <- subset_taxa(phyR.1min,!taxa_names(phyR.1min)%in%c(taxa_names(phyR.contamFreq)))


#########################################################################################################
########## Pre-process data to remove low count OTUs, contaminant OTUs, and duplicate samples  ##########
#########################################################################################################
phyR.sin<-subset_samples(phyR.trim1,Age!="M")
phyR.10min<-subset_samples(phyR.sin,sample_sums(phyR.sin)>=10000) # keep samples with >=10000 reads (physeq.minimum 10,000)
dt.phyPru = fast_melt(phyR.10min) # make data table from physeq object (data table. physeq Pruned)
prev.phyPru = dt.phyPru[, list(Prevalence = sum(count >= var.maxNeg),
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.Pres = prev.phyPru[(Prevalence > 0 & MaxCount >= var.maxNeg), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)
phyR.pru2<-prune_taxa(ls.Pres,phyR.10min) # remove identified taxa from phyloseq object (physeq.Pruned 2)
phyT.Tform<-transform_sample_counts(phyR.pru2, function(x) x/sum(x)) # transform raw counts to fraction (physeq.Transform)

### Add new column to map data. 'Header' column will be used for many figures
sample_data(phyT.Tform)$Header<-with(sample_data(phyT.Tform),
                                     ifelse(AnimalSource=="USA","",
                                            ifelse(AnimalSource=="BBEZ","",
                                                   ifelse(AnimalSource=="GrotonMA","MA",
                                                          ifelse(AnimalSource=="CarogaNY","NY",
                                                                 ifelse(AnimalSource=="Wlot","CT",
                                                                        ifelse(AnimalSource=="MtSnowVT","VT",
                                                                               ifelse(AnimalSource=="Schoolhouse_Brook","CT",     
                                                                                      as.character(AnimalSource)))))))))
### Add new column to map data. 'Da1fb' column will be used to merge some time points together
sample_data(phyT.Tform)$Da1Fb<-with(sample_data(phyT.Tform),
                                    ifelse(Da1F=="8","7",
                                    ifelse(Da1F=="Unk","4",
                                    ifelse(Da1F=="31","30",
                                    ifelse(Da1F=="35","30",
                                    ifelse(Da1F=="90","90+",
                                    ifelse(Da1F=="94","90+",
                                    ifelse(Da1F=="99","90+",
                                    ifelse(Da1F=="100","90+",
                                    ifelse(Da1F=="101","90+",
                                    ifelse(Da1F=="108","90+",
                                    ifelse(Da1F=="110","90+",
                                    ifelse(Da1F=="113","90+",
                                    ifelse(Da1F=="116","90+",
                                    ifelse(Da1F=="130","90+",
                                    ifelse(Da1F=="165","90+",
                                    ifelse(Da1F=="193","90+",
                                    ifelse(Da1F=="215","90+",
                                    ifelse(Da1F=="265","90+",
                                    as.character(Da1F)))))))))))))))))))) 
### Add new column to map data. 'Da1fb' column will be used to merge some time points together
sample_data(phyR.sin)$Da1Fb<-with(sample_data(phyR.sin),
                                  ifelse(Da1F=="8","7",
                                         ifelse(Da1F=="31","30",
                                                ifelse(Da1F=="35","30",
                                                       ifelse(Da1F=="82","90+",
                                                       ifelse(Da1F=="90","90+",
                                                              ifelse(Da1F=="92","90+",
                                                                     ifelse(Da1F=="94","90+",
                                                                            ifelse(Da1F=="96","90+",
                                                                                   ifelse(Da1F=="97","90+",
                                                              ifelse(Da1F=="99","90+",
                                                                     ifelse(Da1F=="100","90+",
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
###############################################
########## Initial sample subsetting ##########
###############################################
# Thank you to jeffkimbrel for the next two lines! #
ls.taxaDC <- setdiff(taxa_names(phyT.Tform), taxa_names(phyR.allc)) # find the difference between taxa in contaminant list and physeqAn (taxa decontam)
phyT.leech <- prune_taxa(ls.taxaDC, phyT.Tform) # keep only taxa not in taxaNC (physeq decontam)
write.table(psmelt(phyT.leech), "tableMelt_Tot.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included

phyT.USA<-subset_samples(phyT.leech,AnimalSource=="USA" & Subproject!="Deposit" & Da2F%in%c("none","0"))
phyT.hv<-subset_samples(phyT.leech,Taxonomic_ID=="Hverbana")
phyT.md<-subset_samples(phyT.leech,Taxonomic_ID=="Mdecora") # subset containing only Macrobdella samples
phyT.mdCT<-subset_samples(phyT.md,AnimalSource=="Wlot") # subset containing only W-lot samples (CT Macrobdella)
phyT.mdMA<-subset_samples(phyT.md,AnimalSource=="GrotonMA") # subset containing only MA samples (MA Macrobdella)
phyT.mdNY<-subset_samples(phyT.md,AnimalSource=="CarogaNY") # subset containing only NY samples (NY Macrobdella)
phyT.mdVT<-subset_samples(phyT.md,AnimalSource=="MtSnowVT") # subset containing only VT samples (VT Macrobdella)

##################################################
########## Effect of extraction method  ##########
##################################################
bray.md<- phyloseq::distance(phyT.md, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(phyT.md)) # make
perm.md<-adonis(bray.md ~ Kit_Name, data = sd.md, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.md

########################################
########## Initial PERMANOVA  ##########
########################################
phyT.adonis1<-phyT.leech # assign data for permanova testing (phyloseq transformed . adonis function)
dist.adonis1<- phyloseq::distance(phyT.adonis1, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.adonis1 <- data.frame(sample_data(phyT.adonis1)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.adonis1<-adonis(dist.adonis1 ~ Taxonomic_ID, data = sd.adonis1, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.adonis1

#######################################################################
########## PERMANOVA for taxonomic groups merged at "Order"  ##########
#######################################################################
phyT.adonisOrder<-tax_glom(phyT.leech,taxrank="Order") # assign data for permanova testing (phyloseq transformed . adonis function)
dist.adonisOrder<- phyloseq::distance(phyT.adonisOrder, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.adonisOrder <- data.frame(sample_data(phyT.adonisOrder)) # make a data frame from the sample_data (sample data . macrobdella decora)
perm.adonisOrder<-adonis(dist.adonisOrder ~ Taxonomic_ID, data = sd.adonisOrder, method = "bray") # Adonis test (permanova . macrobdella decora)
perm.adonisOrder

#####################################################################################
########## FIGURE 1 #################################################################
#####################################################################################
#### Host Species and Sample Type NMDS ###
dsNMDS<-phyT.leech # define data for analysis
sample_data(dsNMDS)$Sample_Type = factor(sample_data(dsNMDS)$Sample_Type, levels = c("ILF","Intestinum","Bladder"))
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination
nmds.SpeciesType<-plot_ordination(dsNMDS, ord, shape="Sample_Type",color="Taxonomic_ID") + 
  geom_point(aes(fill=Taxonomic_ID),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  scale_shape_manual(values=c(16,17,3)) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position='top',legend.title = element_blank()) +
  guides(color = "none", size="none",fill="none")
#nmds.SpeciesType # uncomment to print plot
### Host Species and Sample Type by Order NMDS ###
dsNMDSorder<-tax_glom(dsNMDS,taxrank="Order")
distOrdO = phyloseq::distance(dsNMDSorder, method = c(mNMDS)) # calculate distances
ordO = ordinate(dsNMDSorder, method = "NMDS", distance = distOrdO) # calculate ordination
nmds.SpeciesTypeOrder<-plot_ordination(dsNMDSorder, ordO, shape="Sample_Type",color="Taxonomic_ID") + 
  geom_point(aes(fill=Taxonomic_ID),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=pal.CB) +
  scale_shape_manual(values=c(16,17,3)) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position='top',legend.title = element_blank()) +
  guides(shape = "none")
#nmds.SpeciesTypeOrder # uncomment to print plot
ggsave(plot_grid(nmds.SpeciesType, nmds.SpeciesTypeOrder, labels = "AUTO",ncol=2), filename="NMDS/plotNMDS_fig1.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 1 to .eps file
########## end FIGURE 1 ##############################################################

##########################################
########## H.verbana PERMANOVA  ##########
##########################################
phyT.adonisHv<-subset_samples(phyT.leech, Taxonomic_ID=="Hverbana") # assign data for permanova testing (phyloseq transformed . adonis Hirudo verbana)
dist.adonisHv<- phyloseq::distance(phyT.adonisHv, method = "bray") # Calculate bray curtis distance matrix (distance . adonis Hirudo verbana)
sd.adonisHv <- data.frame(sample_data(phyT.adonisHv)) # make a data frame from the sample_data (sample data . adonis Hirudo verbana)
perm.adonisHv<-adonis(dist.adonisHv ~ Sample_Type * Da1Fb * AnimalSource * ParentLot, data = sd.adonisHv, method = "bray") # Adonis test (permanova . adonis Hirudo verbana)
perm.adonisHv

#########################################
########## M.decora PERMANOVA  ##########
#########################################
phyT.adonisMd<-subset_samples(phyT.leech, Taxonomic_ID=="Mdecora") # assign data for permanova testing (phyloseq transformed . adonis Macrobdella decora)
dist.adonisMd<- phyloseq::distance(phyT.adonisMd, method = "bray") # Calculate bray curtis distance matrix (distance . adonis Macrobdella decora)
sd.adonisMd <- data.frame(sample_data(phyT.adonisMd)) # make a data frame from the sample_data (sample data . adonis Macrobdella decora)
perm.adonisMd<-adonis(dist.adonisMd ~ Sample_Type * WildMonth * Da1Fb * AnimalSource, data = sd.adonisMd, method = "bray") # Adonis test (permanova . adonis Macrobdella decora)
perm.adonisMd

bray.mdInt<- phyloseq::distance(subset_samples(phyT.adonisMd,Sample_Type=="Intestinum"), method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.mdInt <- data.frame(sample_data(subset_samples(phyT.adonisMd,Sample_Type=="Intestinum"))) # make a data frame from the sample_data (sample data . macrobdella decora)
pairwise.adonis(bray.mdInt,sd.mdInt$Da1Fb)

#####################################################################################
########## FIGURE 2 #################################################################
#####################################################################################
phyR.adult<-subset_samples(phyR.pru2,sample_names(phyR.sin)%in%c(sample_names(phyT.leech)))
phyR.ILF<-subset_samples(phyR.adult,Sample_Type%in%c("ILF")) # subset phyR.md to include only ILF samples
phyR.Int<-subset_samples(phyR.adult,Sample_Type=="Intestinum")
########## Alpha Diversity ILF ##########
phyR.adILF<-subset_samples(phyR.ILF,!Da2F%in%c("1","2","4","7","8","28") & !Da1F%in%c("1","2","4","7")) # keep samples that are less than 1 or more than 30 days after feeding (Phyloseq raw . alpha diversity ILF)) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed
pDiv.ILF<-plot_richness(phyR.adILF, x = "Taxonomic_ID",measures=c("Shannon")) +
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,3)) +
  theme_bw() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none") +
  scale_x_discrete(labels=c("Hverbana" = "Hirudo verbana", "Mdecora" = "Macrobdella decora")) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
#pDiv.ILF # uncomment to print plot
########## Alpha Diversity Intestinum ##########
phyR.adInt<-subset_samples(phyR.Int,!Da2F%in%c("1","2","4","7","8","28") & !Da1F%in%c("1","2","4","7")) # keep samples that are less than 1 or more than 30 days after feeding (Phyloseq raw . alpha diversity Intestinum)
pDiv.Int<-plot_richness(phyR.adInt, x="Taxonomic_ID", measures=c("Shannon")) + 
  geom_boxplot(aes(fill=Taxonomic_ID)) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,3)) +
  theme_bw() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none") +
  scale_x_discrete(labels=c("Hverbana" = "Hirudo verbana", "Mdecora" = "Macrobdella decora")) +
  scale_fill_manual(values=c("#808080","#f9f9f9"))
#pDiv.Int # uncomment to print plot
ggsave(plot_grid(pDiv.ILF, pDiv.Int, labels = "AUTO"), filename="Plots/plotAlpha_fig2.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 2 to .eps file
########## end FIGURE 2 ##############################################################


##### Shifts in OTU abundance between H.verbana ILF and M.decora ILF #####
#phyR.adILF<-subset_samples(phyR.ILF,!Da2F%in%c("1","2","4","7","8","28") & !Da1F%in%c("1","2","4","7")) # keep samples that are less than 1 or more than 30 days after feeding (Phyloseq raw . alpha diversity ILF)) # assign data for alpha diversity analysis (alpha diversitydata). This data must be non-transformed
dsq.start = phyloseq_to_deseq2(phyR.adILF, ~ Taxonomic_ID)
geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.adILF)[rownames(sigtab), ], "matrix"))
ls.mdChange<-levels(sigtab$Number)

phyT.mdChange<-prune_taxa(ls.mdChange,subset_samples(phyT.Tform,sample_names(phyT.Tform)%in%c(sample_names(phyR.adILF))))
dt.mdChange = fast_melt(phyT.mdChange) # make data table from physeq object (data table. physeq Pruned)
prev.mdChange = dt.mdChange[, list(Prevalence = sum(count >= .001),
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)
ls.sigChange = prev.mdChange[(Prevalence > 3 & MaxCount >= .01), TaxaID] # Make list of OTUs present in dataset and having at least maxNeg reads in at least one sample (list.high prevalence)

phyR.dsq<-prune_taxa(ls.sigChange,phyR.adILF)
dsq.start = phyloseq_to_deseq2(phyR.dsq, ~ Taxonomic_ID)
geoMeans = apply(counts(dsq.start), 1, gm_mean)
esf.start = estimateSizeFactors(dsq.start, geoMeans = geoMeans)
dsq.start = DESeq(esf.start, fitType="local")

res = results(dsq.start)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.001
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyR.dsq)[rownames(sigtab), ], "matrix"))
posigtab = sigtab[abs(sigtab[, "log2FoldChange"]) > 2, ]
posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Order","Family", "Genus","Genus2","Number","RDP")]
ls.02<-levels(posigtab$Number)

pL2fc.ILF<-ggplot(posigtab, aes(y=Genus2, x=log2FoldChange, color=Order)) +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=6) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_manual(values=pal.pairMini) 
#pL2fc.ILF # uncomment to print plot

##################################################
########## M.decora Adonis #######################
##################################################
### Adonis analysis of ILF and intestinum samples
bray.md<- phyloseq::distance(subset_samples(phyT.md,Sample_Type!="Bladder"), method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(subset_samples(phyT.md,Sample_Type!="Bladder"))) # make a data frame from the sample_data (sample data . macrobdella decora)
pairwise.adonis(bray.md,paste(sd.md$Sample_Type,sd.md$AnimalSource))
### Adonis analysis of just CT and MA samples
bray.md<- phyloseq::distance(subset_samples(phyT.md,Sample_Type!="Bladder" & AnimalSource%in%c("Wlot","GrotonMA")), method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)
sd.md <- data.frame(sample_data(subset_samples(phyT.md,Sample_Type!="Bladder" & AnimalSource%in%c("Wlot","GrotonMA")))) # make a data frame from the sample_data (sample data . macrobdella decora)
pairwise.adonis(bray.md,paste(sd.md$Sample_Type,sd.md$AnimalSource))

########## DaF in a loop ##########
### Macrobdella decora ILF by Da1F ###
phyT.mdDaFilf<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA),Sample_Type=="ILF" & Da2F%in%c("none","Unk","0")) # merge fed Mdecora ILF samples (phyloseq Transformed . macrobdella decora Days after Feeding ILF)
ls.Da1F<-levels(factor(sample_data(phyT.mdDaFilf)$Da1Fb))
mtx.sigDaFILF<-matrix(ncol=length(ls.Da1F),nrow=length(ls.Da1F))
dimnames(mtx.sigDaFILF) = list(c(ls.Da1F),c(ls.Da1F))
for(day1 in c(ls.Da1F)){
  for(day in c(ls.Da1F)){ 
    if(day==day1){next}
    loop.mdILF<-subset_samples(phyT.mdDaFilf,Da1Fb%in%c(day1,day)) # Examine only ILF data
    bray.mdILF<- phyloseq::distance(loop.mdILF, method = "bray") # Calculate bray curtis distance matrix
    df.mdILF <- data.frame(sample_data(loop.mdILF)) # make a data frame from the sample_data
    perm.mdILF<-adonis(bray.mdILF ~ Da1Fb, data = df.mdILF, method = "bray") # Adonis test
    sig.mdILF<-as.data.frame(perm.mdILF$aov.tab)["Da1F", "Pr(>F)"]
    print(perm.mdILF$aov.tab)
    mtx.sigDaFILF[day1,day] <- as.numeric(sig.mdILF)
  }
}
write.table(mtx.sigDaFILF, "tableSig_DaFMd_ILF.csv", sep=",",row.names=TRUE,col.names=TRUE) # export table to csv, row and column names included


##################################################
########## Define core and common OTUs  ##########
##################################################
cMin<-0.03 # define minimum to count prevalence in a sample
### Hirudo ###
phyR.hv<-subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%c(sample_names(phyT.hv))) # create raw reads phyloseq object of Hirduo verbana samples (physeq Raw . hirudo verbana)
phyR.hvBlad2<-subset_samples(phyR.hv,Sample_Type=="Bladder") # keep only bladder samples (physeq Raw . hirudo verbana Bladder 2)
phyR.hvBlad<-subset_taxa(phyR.hvBlad2,!Genus2%in%c("Mucinivorans","Aeromonas","Bacteroides")) # remove Mucinivorans, Bacteroides, and Aeromonas OTUs (physeq Raw . hirudo verbana Bladder)
phyT.hvBlad<-transform_sample_counts(phyR.hvBlad, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . hirudo verbana Bladder)
fm.hvBlad = fast_melt(phyT.hvBlad) # (fast melt Hirudo verbana bladder)
prevdt.hvBlad = fm.hvBlad[,list(Prevalence = sum(count >= .001), 
                                TotalCount = sum(count),
                                MinCount = min(count),
                                MaxCount = max(count),
                                MedCount = median(count),
                                Med2Count = median(count[count!=0]),
                                PrevPer = sum(count >= .001) / nsamples(phyT.hvBlad)),
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

###Macrobdella		
phyR.md<-subset_samples(phyR.pru2,sample_names(phyR.pru2)%in%c(sample_names(phyT.md))) # create raw reads phyloseq object of Macrobdella decora samples (physeq Raw . macrobdella decora)
phyR.mdBlad2<-subset_samples(phyR.md,Sample_Type=="Bladder") # keep only bladder samples (physeq Raw . macrobdella decora Bladder 2)
phyR.mdBlad<-subset_taxa(phyR.mdBlad2,!Genus2%in%c("Mucinivorans","Aeromonas","Bacteroides")) # remove Mucinivorans, Bacteroides, and Aeromonas OTUs (physeq Raw . macrobdella decora Bladder)
phyT.mdBlad<-transform_sample_counts(phyR.mdBlad, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . macrobdella decora Bladder)
phyT.mdBlad2<-subset_samples(phyT.mdBlad,!Da1F%in%c("1","4","7"))
phyT.mdBlad3<-subset_samples(phyT.mdBlad2,!sample_names(phyT.mdBlad2)%in%c("MdMA0dF040618EMfBD","MdMA0dF040618EMdBD722","MdMA0dF040618EMeBD"))
fm.mdBlad = fast_melt(phyT.mdBlad3) # (fast melt Hirudo verbana bladder)
prevdt.mdBlad = fm.mdBlad[,list(Prevalence = sum(count >= .001), 
                                TotalCount = sum(count),
                                MinCount = min(count),
                                MaxCount = max(count),
                                MedCount = median(count),
                                Med2Count = median(count[count!=0]),
                                PrevPer = sum(count >= .001) / nsamples(phyT.hvBlad)),
                          by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana Blad)
colnames(prevdt.mdBlad)[colnames(prevdt.mdBlad)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdBlad,tax_table(phyT.mdBlad)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdBlad <- Reduce(function(...) merge(..., all = T), list(prevdt.mdBlad,tax_table(phyT.mdBlad)[,c("Number","Genus2")])) # m     
prevdt.mdBladhi<-prevdt.mdBlad[which(PrevPer>0),][order(-PrevPer)]
write.table(prevdt.mdBladhi, "tablePrev_mdBlad.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded

phyT.mdBWint<-subset_samples(phyT.md,Sample_Type=="Intestinum") # (Hirudo verbana intestinum)
fm.mdInt = fast_melt(phyT.mdBWint) # (fast melt Hirudo verbana intestinum)
prevdt.mdInt = fm.mdInt[, list(Prevalence = sum(count >= cMin),
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   
                               Med2Count = median(count[count!=0]),   
                               PrevPer = sum(count >= cMin) / nsamples(phyT.mdBWint)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Hirudo verbana ILF)
colnames(prevdt.mdInt)[colnames(prevdt.mdInt)=="TaxaID"]<-"Number" # rename 'TaxaID' column
lapply(list(prevdt.mdInt,tax_table(phyT.mdBWint)[,c("Number","Genus2")]), function(i) setkey(i, Number)) # set key for merge function
prevdt.mdInt <- Reduce(function(...) merge(..., all = T), list(prevdt.mdInt,tax_table(phyT.mdBWint)[,c("Number","Genus2")])) # m     
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

phyT.vtILF<-subset_samples(phyT.mdVT,Sample_Type=="ILF") # (Vermont ILF)
fm.vtILF = fast_melt(phyT.vtILF) # (fast melt Vermont ILF)
prevdt.vtILF = fm.vtILF[, list(Prevalence = sum(count >= cMin), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count),   
                               MedCount = median(count),   Med2Count = median(count[count!=0]),   PrevPer = sum(count >= cMin) / nsamples(phyT.hvBlad)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence data table Vermont ILF)

########## List core OTUs ##########
cpR<-c(.9) # define the level for determining core OTUs, reported as a percent of the samples tested (core percent)

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
ls.coreNyBlad = prevdt.nyBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.nyBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella New York Bladder
ls.coreNY = unique(c(ls.coreNyBlad,ls.coreNyILF))
### VT Macrobdella decora
ls.coreVtILF = prevdt.vtILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.vtILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreVtBlad = prevdt.vtBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.vtBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of core OTUs for Macrobdella Vermont Bladder
ls.coreVT = unique(c(ls.coreVtBlad,ls.coreVtILF))
### Compiled Macrobdella decora
ls.coreMdILF = prevdt.mdILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreMdILF0 = prevdt.mdILF[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdILF,Da1Fb=="0")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont ILF (
ls.coreMdInt = prevdt.mdInt[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBWint,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont Intestinum
ls.coreMdBlad = prevdt.mdBlad[(Prevalence >= cpR*nsamples(subset_samples(phyT.mdBlad3,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of core OTUs for Macrobdella Vermont Bladder
ls.coreMd = unique(c(ls.coreCT,ls.coreMA,ls.coreNY)) # removed ls.coreVT due to this being only 1 sample
ls.coreFeedILF = unique(c(ls.coreCtILF,ls.coreMaILF))
ls.coreFeedInt = unique(c(ls.coreCtInt,ls.coreMaInt))
ls.coreDaF<-unique(c(ls.coreMdILF,ls.coreMdILF0,ls.coreMdInt,ls.coreHvILF0,ls.coreHvILF,ls.coreHvInt))
ls.coreTot = unique(c(ls.coreMd,ls.coreHv,ls.coreCT,ls.coreMA,ls.coreNY,ls.coreMdILF,ls.coreMdILF0,ls.coreMdInt,ls.coreMdBlad)) # removed ls.coreVT due to this being only 1 sample

########## List common OTUs ##########
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
ls.comNyBlad = prevdt.nyBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.nyBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella New York Bladder
ls.comNY = unique(c(ls.comNyBlad,ls.comNyILF))
### VT Macrobdella decora ###
ls.comVtILF = prevdt.vtILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.vtILF,Sample_Type=="ILF")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comVtBlad = prevdt.vtBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.vtBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), TaxaID] # Make list of common OTUs for Macrobdella Vermont Bladder
ls.comVT = unique(c(ls.comVtBlad,ls.comVtILF))
### Compiled Macrobdella decora ###
ls.comMdILF = prevdt.mdILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdILF,Sample_Type=="ILF")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comMdILF0 = prevdt.mdILF[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdILF,Da1Fb=="0")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont ILF (
ls.comMdInt = prevdt.mdInt[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdBWint,Sample_Type=="Intestinum")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont Intestinum
ls.comMdBlad = prevdt.mdBlad[(Prevalence >= cpM*nsamples(subset_samples(phyT.mdBlad,Sample_Type=="Bladder")) & MaxCount >= cMin), Number] # Make list of common OTUs for Macrobdella Vermont Bladder
ls.comMd = unique(c(ls.comCT,ls.comMA,ls.comNY)) # removed ls.comVT due to this being only 1 sample
ls.comFeedILF = unique(c(ls.comCtILF,ls.comMaILF))
ls.comFeedInt = unique(c(ls.comCtInt,ls.comMaInt))
ls.comTot = unique(c(ls.comMd,ls.comHv,ls.comCT,ls.comMA,ls.comNY,ls.comMdILF,ls.comMdILF0,ls.comMdInt,ls.comMdBlad)) # removed ls.comVT due to this being only 1 sample

########## Print table of core taxa with columns indicating which sample type each is present in ##########
### Add core columns to taxa table ###
phyT.coreTot<-prune_taxa(ls.coreTot,phyT.leech) # keep only taxa from the identified 'Core' 
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
write.table(dt.TAXcore2, "tableTax_coreCompile.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded


########## Print table of common taxa with columns indicating which sample type each is present in ##########
### Add core columns to taxa table ###
phyT.comTot<-prune_taxa(ls.comTot,phyT.leech) # keep only taxa from the identified 'Core' 
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
write.table(dt.TAXcom2, "tableTax_commonCompile.csv", row.names=FALSE,sep=",") # export table to csv, row names excluded



#####################################################################################
########## FIGURE 4 #################################################################
#####################################################################################
### Collection Month NMDS - ILF ###
dsNMDS<-subset_samples(phyT.md,Sample_Type=="ILF") # define data for analysis
sample_data(dsNMDS)$WildMonth = factor(sample_data(dsNMDS)$WildMonth, levels = c("April","June","July","August","September","October")) # force Months to be in chronological order (instead of alphabetical default)
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS, method = "NMDS", distance = distOrd) # calculate ordination
nmds.collILF<-plot_ordination(dsNMDS, ord, color = "WildMonth") + 
  geom_point(aes(fill=WildMonth),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#9bcdff","#ffb4b4","#ff0101","#810000","#810000","#005ab4")) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position="none")
#nmds.collILF # uncomment to print plot
nmds.legend<-plot_ordination(dsNMDS, ord, color = "WildMonth") + 
  scale_color_manual(values=c("#9bcdff","#ffb4b4","#ff0101","#810000","#810000","#005ab4")) +
  theme(text=element_text(size=10),legend.position='top',legend.title = element_blank())

### Collection Month NMDS - Intestinum ###
dsNMDSint<-subset_samples(phyT.md,Sample_Type=="Intestinum") # define data for analysis
sample_data(dsNMDSint)$WildMonth = factor(sample_data(dsNMDSint)$WildMonth, levels = c("April","June","July","August","September","October")) # force Months to be in chronological order (instead of alphabetical default)
distOrdInt = phyloseq::distance(dsNMDSint, method = c(mNMDS)) # calculate distances
ordInt = ordinate(dsNMDSint, method = "NMDS", distance = distOrdInt) # calculate ordination
nmds.collInt<-plot_ordination(dsNMDSint, ordInt, color = "WildMonth") + 
  geom_point(aes(fill=WildMonth),size=2) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#9bcdff","#ffb4b4","#ff0101","#810000","#810000","#005ab4")) +
  theme_bw() +
  theme(text=element_text(size=10),legend.position="none")
#nmds.collInt # uncomment to print plot
ggsave(plot_grid(plot_grid(nmds.collILF, nmds.collInt, labels = "AUTO"),cowplot::get_legend(nmds.legend),nrow=2,rel_heights = c(5, 1)), filename="NMDS/plotNMDS_fig4.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 4 to .eps file
########## end FIGURE 4 ##############################################################


#####################################################################################
########## FIGURE 5 #################################################################
#####################################################################################
### ILF Mdecora ###
phyT.DaFilf<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA,phyT.hv),Sample_Type=="ILF" & Da2F%in%c("none","Unk","0")) # merge fed Hverbana and Mdecora ILF samples (phyloseq Transformed . Days after Feeding ILF)
phyT.DaFilfp<-prune_taxa(c(ls.coreDaF),phyT.DaFilf) # keep only core ILF and intestinum taxa (phyloseq Transformed . Days after Feeding ILF pruned)
phyT.mdBWilf<-subset_samples(phyT.DaFilfp,Taxonomic_ID=="Mdecora") # keep only Mdecora samples (phyloseq Transformed . macrobdella decora Box+Whisker ILF)
ls.days<-unique(sample_data(phyT.mdBWilf)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days)

for (d in c(ls.days)){
  phyT.AAa<-subset_samples(phyT.mdBWilf,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.mdBWilf,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= .001)), by = TaxaID] # make simple table listing 'TaxaID and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of OTUs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[eval(f),]<-0 # change all prevalence values to 0 for OTUs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all OTUs in ls.low
  phyT.mdBWilf<-merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . macrobdella decora Box+Whisker ILF)
} # repeat this loop for all days after feeding (d in ls.days) 
### ILF Hverbana ###
phyT.hvBWilf<-subset_samples(phyT.DaFilfp,Taxonomic_ID=="Hverbana") # keep only Hverbana samples (phyloseq Transformed . hirudo verbana Box+Whisker ILF)
ls.days<-unique(sample_data(phyT.hvBWilf)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days)

for (d in c(ls.days)){
  phyT.AAa<-subset_samples(phyT.hvBWilf,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.hvBWilf,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= .001)), by = TaxaID] # make simple table listing 'TaxaID and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of OTUs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[eval(f),]<-0 # change all prevalence values to 0 for OTUs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all OTUs in ls.low
  phyT.hvBWilf<-merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . hirudo verbana Box+Whisker ILF)
} # repeat this loop for all days after feeding (d in ls.days) 
### Intestinum M. decora ###
phyT.mdDaFint<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA),Sample_Type=="Intestinum" & Da2F%in%c("none","Unk","0")) # merge fed Mdecora Intestinum samples (phyloseq Transformed . macrobdella decora Days after Feeding intestinum)
phyT.mdBWint<-prune_taxa(c(ls.coreDaF),phyT.mdDaFint) # keep only core ILF and intestinum taxa (phyloseq Transformed . macrobdella decora Box+Whisker intestinum)
ls.days<-unique(sample_data(phyT.mdBWint)$Da1Fb) # list the days after feeding at which ILF was sampled (list . days)

for (d in c(ls.days)){
  phyT.AAa<-subset_samples(phyT.mdBWint,Da1Fb==eval(d)) # subset samples to those from timepoint d
  phyT.AAb<-subset_samples(phyT.mdBWint,Da1Fb!=eval(d)) # subset samples to those NOT from timepoint d
  dt.AAa = fast_melt(phyT.AAa) # make data table from the samples at timepoint d
  prev.AAa = dt.AAa[, list(Prevalence = sum(count >= .001)), by = TaxaID] # make simple table listing 'TaxaID and Prevalence'
  ls.low = prev.AAa[(Prevalence <= 1), TaxaID] # make list of OTUs present in ≤ 1 sample from timepoint d (list . low prevalence)
  for (f in c(ls.low)){
    otu_table(phyT.AAa)[eval(f),]<-0 # change all prevalence values to 0 for OTUs present in ≤ 1 sample from timepoint d 
  } # repeat this loop for all OTUs in ls.low
  phyT.mdBWint<-merge_phyloseq(phyT.AAa,phyT.AAb) # merge the modified data from timepoint d to the unmodified data from NOT timepoint d (phyloseq Transformed . macrobdella decora Box+Whisker intestinum)
} # repeat this loop for all days after feeding (d in ls.days) 

phyT.BWmerge<-merge_phyloseq(phyT.mdBWilf,phyT.hvBWilf,phyT.mdBWint) # merge the modified data from Mdecora ILF, Hverbana ILF, and Mdecora intestinum (phyloseq Transformed . Box + Whisker merged)
sample_data(phyT.BWmerge)$Da1Fb = factor(sample_data(phyT.BWmerge)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+")) # Reorder Da1Fb
sample_data(phyT.BWmerge)$Taxonomic_ID = factor(sample_data(phyT.BWmerge)$Taxonomic_ID, levels = c("Hverbana","Mdecora")) # Reorder Taxonomic_ID

dt.BWmerge <- psmelt(phyT.BWmerge) # create data.table from phyloseq object (data table . Box+Whisker merged)
pBW.DaF <- ggplot(dt.BWmerge) +
  stat_boxplot(data=dt.BWmerge, aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  facet_grid(Taxonomic_ID+Sample_Type~Order+Genus2, scales="free_x",space="free") +
  theme_bw() +
  geom_text(x=2,y=1,label="*") +
  scale_y_log10(limits=c(.001,1)) + 
  xlab("Days after feeding") +
  theme(text=element_text(size=10),strip.text=element_text(size=6), axis.text.x=element_text(angle=90,hjust=1), legend.position="none") +
  scale_fill_manual(values=brewer.pal(8,"Greys")) 
pBW.DaF # uncomment to print plot
ggsave(plot_grid(pBW.DaF, labels = "AUTO"), filename="Plots/plotBW_Fig5.eps", device="eps", dpi="retina",width=6.87,height=6.87,units="in") # Save FIGURE 5 to .eps file
########## end FIGURE 5 ##############################################################

############################################################################################
########## FIGURE 6 ########################################################################
############################################################################################
########## Alpha Diversity ILF ##########
#phyT.DaFilf<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA,phyT.hv),Sample_Type=="ILF" & Da2F%in%c("none","Unk","0")) # merge fed Hverbana and Mdecora ILF samples (phyloseq Transformed . Days after Feeding ILF)
phyR.mdDaFilf<-subset_samples(phyR.sin,sample_names(phyR.sin)%in%c(sample_names(phyT.DaFilf)) ) # collect raw data values from only Mdecora samples (phyloseq Raw . macrobdella decora Days after Feeding ILF)
sample_data(phyR.mdDaFilf)$Da1Fb = factor(sample_data(phyR.mdDaFilf)$Da1Fb, levels = c(0,1,2,4,7,30,"90+")) # Reorder Da1Fb

pDiv.mdDaFilf<-plot_richness(phyR.mdDaFilf, x = "Da1Fb",measures=c("Shannon")) +
  geom_boxplot(aes(fill=Da1Fb)) +
  geom_jitter(width=0.15) +
  facet_grid(~Taxonomic_ID, scales="free_x",space="free") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none") +
  scale_fill_manual(values=brewer.pal(7,"Greys"))
#pDiv.mdDaFilf # uncomment to print plot

##### weighted Unifrac distance by Da1F #####
# Thank you to jeffkimbrel ! #
phyT.mdILFdaf<-subset_samples(phyT.DaFilf, Taxonomic_ID=="Mdecora")
sample_data(phyT.mdILFdaf)$name<-sample_names(phyT.mdILFdaf)
dist.mdILFdaf<- phyloseq::distance(phyT.mdILFdaf, method = "bray") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)

# remove self-comparisons
md.mdILFdaf = melt(as.matrix(dist.mdILFdaf)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd.mdILFdaf = sample_data(phyT.mdILFdaf) %>%
  select("name","Da1Fb","Header") %>%
  mutate_if(is.factor,as.character) 

# combined distances with sample data
colnames(sd.mdILFdaf) = c("Var1", "Type1","Source")
md.mdILFdaf2 = left_join(md.mdILFdaf, sd.mdILFdaf, by = "Var1")

colnames(sd.mdILFdaf) = c("Var2", "Da1F","Source")
md.mdILFdaf3 = left_join(md.mdILFdaf2, sd.mdILFdaf, by = "Var2")

wu.sd2<-md.mdILFdaf3[md.mdILFdaf3$Type1==0,]
wu.sd2$Da1F = factor(wu.sd2$Da1F, levels = c(0,1,2,4,7,30,"90+")) # Reorder Da1Fb
colnames(wu.sd2)[colnames(wu.sd2)=="value"] <- "distance"

pBray.mdILFdaf<-ggplot(wu.sd2, aes(x = Da1F, y = distance)) +
  theme_bw() +
  geom_violin(aes(x = Da1F, y = distance)) +
  geom_jitter(width = 0.2) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Bray-curtis distance") +
  xlab("Days after feeding") +
  theme_bw()
pBray.mdILFdaf # uncomment to print plot

ggsave(plot_grid(pDiv.mdDaFilf, pBray.mdILFdaf, labels = "AUTO", nrow=2), filename="Plots/plotAlpha_fig6.eps", device="eps", dpi="retina",width=6.87,height=10,units="in") # Save FIGURE 6 to .eps file
########## end FIGURE 6 ##############################################################


#####################################################################################
########## SUPPLEMENTAL FIGURE 1 ####################################################
#####################################################################################
### ILF NMDS - Collection Site ###
dsNMDS.mdILF<-subset_samples(phyT.md,Sample_Type=="ILF") # define data for analysis
sample_data(dsNMDS.mdILF)$AnimalSource = factor(sample_data(dsNMDS.mdILF)$AnimalSource, levels = c("Wlot","GrotonMA","CarogaNY","MtSnowVT")) # force AnimalSource order (instead of alphabetical default)
mNMDS<-"unifrac" # define metric for analysis
distOrd = phyloseq::distance(dsNMDS.mdILF, method = c(mNMDS)) # calculate distances
ord = ordinate(dsNMDS.mdILF, method = "NMDS", distance = distOrd) # calculate ordination
nmds.mdILFstate<-plot_ordination(dsNMDS.mdILF, ord, color = "AnimalSource") + 
  geom_point(aes(fill=AnimalSource),size=3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#000000")) +
  theme_bw() +
  theme(text=element_text(size=10), legend.position="none")
#nmds.mdILFstate # uncomment to print plot
### Intestinum/WildTime NMDS ###
dsNMDS.mdInt<-subset_samples(phyT.md,Sample_Type=="Intestinum") # define data for analysis
sample_data(dsNMDS.mdInt)$AnimalSource = factor(sample_data(dsNMDS.mdInt)$AnimalSource, levels = c("Wlot","GrotonMA","CarogaNY","MtSnowVT")) # force AnimalSource order (instead of alphabetical default)
distOrdInt = phyloseq::distance(dsNMDS.mdInt, method = c(mNMDS)) # calculate distances
ordInt = ordinate(dsNMDS.mdInt, method = "NMDS", distance = distOrdInt) # calculate ordination
nmds.mdIntState<-plot_ordination(dsNMDS.mdInt, ordInt, color = "AnimalSource") + 
  geom_point(aes(fill=AnimalSource),size=3) +
  stat_ellipse(type = "t", level = 0.95, linetype = 2) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73","#000000")) +
  theme_bw() +
  theme(text=element_text(size=10), legend.position="none")
#nmds.mdIntState # uncomment to print plot
ggsave(plot_grid(nmds.mdILFstate, nmds.mdIntState, labels = "AUTO"), filename="NMDS/plotNMDS_supFig1.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 4 to .eps file
########## end SUPPLEMENTAL FIGURE 1 ################################################


### Alpha Diversity Intestinum ###
#phyT.mdDaFint<-subset_samples(merge_phyloseq(phyT.mdCT,phyT.mdMA),Sample_Type=="Intestinum" & Da2F%in%c("none","Unk","0")) # merge fed Mdecora Intestinum samples (phyloseq Transformed . macrobdella decora Days after Feeding intestinum)

phyR.mdDaFint<-subset_samples(phyR.sin,sample_names(phyR.sin)%in%c(sample_names(phyT.mdDaFint)))
sample_data(phyR.mdDaFint)$Da1Fb = factor(sample_data(phyR.mdDaFint)$Da1Fb, levels = c(0,1,2,4,7,30)) # Reorder Da1Fb

pDiv.mdDaFint<-plot_richness(phyR.mdDaFint, x = "Da1Fb",measures=c("Shannon")) +
  geom_boxplot(aes(fill=Da1Fb)) +
  geom_jitter(width=0.15) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  theme(text=element_text(size=10),axis.text.x=element_text(angle = 0),axis.title.x=element_blank(),legend.position = "none") +
  scale_fill_manual(values=brewer.pal(6,"Greys"))
pDiv.mdDaFint # uncomment to print plot

ggsave(plot_grid(pDiv.ILF, pDiv.Int, labels = "AUTO"), filename="Plots/plotAlpha_Int.eps", device="eps", dpi="retina",width=6.87,units="in") # Save FIGURE 2 to .eps file






##############################
# END # END # END # END # END # END 
