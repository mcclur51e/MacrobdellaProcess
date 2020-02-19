# Most current files
table_map_simple.csv # required for phyloseq processing using code file

md_DADA2.R # complete updated and commented code for processing data as reported in "Macrobdella decora: Old World Leech Gut Microbial Community Structure Conserved in a New World Leech"


# MacrobdellaProcess
Started 2018_0309

Preamble = preamble of packages and commands for processing data. Defines pairBiome set of colors 

preprocess_data = preprocessing steps to clean-up data before analysis. Includes: removes OTUs with < 100 reds in entire data set, removes specified outliers, removes samples with < 10000 reads, removes samples with >=1% Halomonas or Anaplasma reads,  transforms all read data to percentage, and merges OTUs by Genus

calculateCore = calculation of core OTUs. Change cMin to change minimum count cut-off (1%). Change cp to change minimum prevalence cut-off (60%).
