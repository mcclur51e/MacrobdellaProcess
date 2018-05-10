# may want to change 'setwd' line 15

########## Processing Set-up ##########
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
#library(“gridExtra”)
### Call functions for use ###
source("~/Masters/R_tools/taxa_summary",local=TRUE) # load fast_melt function
### Set working directory where processing will be carried out ###
setwd("~/Desktop/Process2018_0327/") # set working directory for data processing. Must contain folder "DataFiles" containing "table_otu.csv", "table_tax.csv", and "table_map.csv" and/or "physeq_manuscript.RData"
########## Processing Set-up ##########
#Black,Green,Red,Blue,Orange,Purple,Yellow,Teal,Magenta,Grey
pairBiome<-c("#1a1a1a",
             "#267326","#39ac39","#b3e6b3",
             "#B80000","#F00000","#FF7777","#ffcccc",
             "#000080","#3333ff","#b3b3ff",
             "#ff6600","#ff9933","#ffcc99",
             "#600080","#ac00e6","#f2ccff",
             "#ffff00","#ffffb3",
             "#006680","#00ccff","#ccf5ff",
             "#800055","#ff00aa","#ffb3e6",
             "#808080","#d9d9d9","#ffffff")

rainbow<-c("#ff00aa","#ac00e6","#3333ff","#39ac39","#ffff00","#ff9933","#F00000")

OTU = otu_table(as.matrix.data.frame(read.csv("DataFiles/table_otu.csv", header = TRUE, row.names = 1, check.names=FALSE) ), taxa_are_rows = TRUE) # reads data.frame into matrix without converting to numbers # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with OTU data for each sample
TAX = tax_table(as.matrix.data.frame(read.csv("DataFiles/table_tax.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads data.frame into matrix without converting to numbers # reads csv file into data.frame with row names in column 1 ### .csv file prepped from "otu_table.txt" with taxa data separated by columns and headers added for Kingdom, Phylum, Class, etc. First column = otu IDs
MAP = sample_data(data.frame(read.csv("DataFiles/table_map.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads csv file into data.frame with row names in column 1
TREE = read.tree("DataFiles/Rep_Set_tree.tree")
physeq = phyloseq(OTU, TAX, MAP, TREE)
save(physeq,file=("DataFiles/physeq_current.RData")) # Save the phyloseq data object in a .RData file 

########## Load most recent data ##########
load("DataFiles/physeq_current.RData") # Use to re-load saved phyloseq data object if "DataFiles/physeq_manuscript.RData" existed previously
########## Load most recent data ##########
