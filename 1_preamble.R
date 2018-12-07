#######################################
########## Processing Set-up ##########
#######################################
### Call libraries for use ###
library("phyloseq") #packageVersion("phyloseq")
library("ggplot2") #packageVersion("ggplot2")
library("ape") #packageVersion("ape")
library("data.table") #packageVersion("data.table")
library("RColorBrewer") # design new color palette for data
library("grid")
library("vegan")
#library("plyr") # may need to turn back 'on' and dplyr 'off'
library("dplyr")
### Call functions for use ###
source("~/Masters/R_tools/taxa_summary",local=TRUE) # load fast_melt function
########## Processing Set-up ##########
#Green,Red,Blue,Orange,Purple,Yellow,Teal,Magenta,Grey
pairBiome<-c("#004200","#267326","#39ac39","#b3e6b3",
             "#B80000","#F00000","#FF7777","#ffcccc",
             "#000080","#0000cd","#8282ff","#cfcfff",
             "#b34700","#ff6600","#ff9933","#ffcc99",
             "#600080","#ac00e6","#e599ff","#f2ccff",
             "#cccc00","#ffff00","#ffffb3",
             "#0090b4","#01cdff","#80e6ff","#ccf5ff",
             "#9a0066","#ff4ec5","#ff9ade","#ffe6f7",
             "#1a1a1a","#808080","#d9d9d9","#ffffff","#808080")
rainbow<-c("#ff00aa","#ac00e6","#3333ff","#0ba29a","#39ac39","#ffff00","#ff9933","#F00000","#800055")

########## Set working directory where processing will be carried out. MAke subdirectories for saving into ##########
setwd("~/Desktop/Rprocess_2018_1130/") # set working directory for data processing. Must contain folder "DataFiles" containing "table_otu.csv", "table_tax.csv", and "table_map.csv" and/or "physeq_manuscript.RData"
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
########## Remove duplicate samples from physeq ##########
ls.dups<-as.character(read.csv("DataFiles/list_dups.csv",header=FALSE,sep=",")$V1) # Import list of samples determined to be duplicates. Last evaluated 2018-11-22. (list.duplicates)
phy.sin<-subset_samples(physeq, !sample_names(physeq)%in%c(ls.dups))# Remove duplicates from list (physeq.single)
save(phy.sin,file=("DataFiles/physeq_noDups.RData")) # Save phy.sin in a .RData file 
########## Remove outlier samples identified in list_ouliers.csv ##########
ls.outliers<-as.character(read.csv("DataFiles/list_outliers.csv",header=FALSE,sep=",")$V1) # Import list of samples determined to be duplicates. Last evaluated 2018-11-22. (list.duplicates)
phy.sin<-subset_samples(phy.sin,!sample_names(physeq)%in%c(ls.outliers)) # Remove outliers
save(phy.sin,file=("DataFiles/physeq_noOutliers.RData")) # Save phy.sin in a .RData file 

########## Load most recent data ##########
### Use this to skip lines 35--48 if the basic phyloseq data sets have previously been produced and saved ###
load("DataFiles/physeq_current.RData") # Use to re-load saved phyloseq data object 
load("DataFiles/physeq_noDups.RData") # Use to re-load saved phyloseq data object
load("DataFiles/physeq_noOutliers.RData") # Use to re-load saved phyloseq data object
