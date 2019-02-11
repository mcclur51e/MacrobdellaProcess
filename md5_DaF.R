########## Macrobdella Da1F plot with 'Other' category ##########
phyT.mdDaF<-merge_phyloseq(phyT.mdCT,phyT.mdMA) # merge fed samples from CT and MA
phyT.mdDaFGI<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF","Intestinum")) # Keep only ILF and intestinum samples
phyT.core.mdDaF<-prune_taxa(ls.coreTot,phyT.mdDaFGI) # keep only core taxa
sample_data(phyT.mdDaFGI)$Da1Fb = factor(sample_data(phyT.mdDaFGI)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1Fb


### Make plot
dt.DaF<-data.table(psmelt(phyT.mdDaFGI))
dt.DaF$Number<-as.character(dt.DaF$Number)
dt.DaF$Genus3<-with(dt.DaF,
                     ifelse(Number%in%c(ls.coreTot),as.character(Genus2),
                     as.character("other")))
# Define levels for OTUs
dt.DaF$Genus3 <- factor(dt.DaF$Genus3, levels=c("Aeromonas","Aeromonas-like","Proteus","Morganella","unk_Enterobacteriaceae",
                                                "Mucinivorans","Mucinivorans-like","Bacteroides","Bacteroides-like","Millionella-like","Pedobacter",
                                                "Christenella-like","Alkaliphilus-like","Proteiniclasticum","Clostridium","Clostridium-like",
                                                "Proteocatella","Proteocatella-like","Butyricicoccus","Butyricicoccus-like","Papillibacter-like","Sporobacter-like",
                                                "Ochrobactrum","Ensifer","Rhizobium","Rhizobium-like","Azospirillum",
                                                "Insolitispirillum-like","Phaeospirillum-like","Phreatobacter-like",
                                                "Bacteriovorax-like","Bdellovibrio-like","Desulfovibrio","Desulfovibrio-like","Cystobacter-like",
                                                "Ramlibacter","Methylopumilus-like",
                                                "Fusobacterium",
                                                "other"))
# Define colors for OTUs. Reference = "rainbow"
genus.color<-c(Aeromonas="#1c541c","Aeromonas-like"="#1c541c",Proteus="#2f8f2f",Morganella="#49c349",unk_Enterobacteriaceae="#bfeabf",
               Mucinivorans="#B80000","Mucinivorans-like"="#B80000",Bacteroides="#F00000","Bacteroides-like"="#F00000","Millionella-like"="#FF7777","Pedobacter"="#ffcccc",
               "Christenella-like"="#000080","Alkaliphilus-like"="#0000cd",Proteiniclasticum="#8282ff",Clostridium="#8282ff","Clostridium-like"="#8282ff",
               Proteocatella="#005f6c","Proteocatella-like"="#005f6c",Butyricicoccus="#00a4bb","Butyricicoccus-like"="#00a4bb","Papillibacter-like"="#1ee3ff","Sporobacter-like"="#bbf7ff",
               Ochrobactrum="#623800",Ensifer="#c47000",Rhizobium="#ff9914","Rhizobium-like"="#ff9914",Azospirillum="#ffddb1",
               "Insolitispirillum-like"="#430059","Phaeospirillum-like"="#7d00a7","Phreatobacter-like"="#cc32ff",
               "Bacteriovorax-like"="#626200","Bdellovibrio-like"="#c4c400",Desulfovibrio="#ffff14","Desulfovibrio-like"="#ffff14","Cystobacter-like"="#ffffb1",
               "Ramlibacter"="#750063","Methylopumilus-like"="#c400a5",
               "Fusobacterium"="#ffb0f3",
               other="#808080") 

pBar.DaF <- ggplot(dt.DaF, aes(x=Replicate,y=Abundance, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~AnimalSource+Da1Fb, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pBar.DaF # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaF), size = "last")), filename="Plots/plotBarStack_Da1F.eps", device="eps", width=12,height=8)


### Macrobdella DaF ILF only
phyT.mdDaFILF<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF")) # Keep only ILF samples
phyT.mdDaFgr<-subset_samples(phyT.mdDaFILF,AnimalSource%in%c("GrotonMA","Wlot"))
phyT.mdILFpr<-prune_taxa(taxa_sums(phyT.mdDaFgr)>.1,phyT.mdDaFgr) # keep taxa with at least 1% of 1 sample
#phyT.mdILFfam<-subset_taxa(phyT.mdILFpr, taxa_sums(phyT.mdILFpr)>.05)

phyT.mdILFfam<-subset_taxa(phyT.mdILFpr, taxa_names(phyT.mdILFpr)%in%c(ls.coreFeedILF))

# Convert phyloseq to data table
dt.ILF<-data.table(psmelt(phyT.mdILFfam))
dt.ILF$Number<-as.character(dt.ILF$Number)
dt.ILF$Genus2<-as.character(dt.ILF$Genus2)
dt.ILF[!dt.ILF$Number%in%as.character(ls.coreMd),]$Genus2<-" "
dt.ILF$Da1Fb<-factor(dt.ILF$Da1Fb, levels = c(0,1,2,4,7,30,"90+",100)) # Reorder Da1Fb

# Make plot
pBar.DaFilf <- ggplot(dt.ILF, aes(x=Sample, y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Family~Da1Fb, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(9,"Set1"))
pBar.DaFilf # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaFilf), size = "last")), filename="Plots/plotBar_DaFILF.eps", device="eps", width=12,height=8)



phyT.mdILFblad<-subset_taxa(phyT.mdILFpr, taxa_names(phyT.mdILFpr)%in%c(setdiff(ls.coreMdBlad,ls.coreMdILF)))
sample_data(phyT.mdILFblad)$Da1Fb = factor(sample_data(phyT.mdILFblad)$Da1Fb, levels = c(0,1,4,7,30,31,35,113,215)) # Reorder Da1Fb

# Convert phyloseq to data table
dt.mdILFblad<-data.table(psmelt(phyT.mdILFblad))
dt.mdILFblad$Number0o<-as.character(dt.mdILFblad$Number)
dt.mdILFblad$Genus2<-as.character(dt.mdILFblad$Genus2)
dt.mdILFblad[!dt.mdILFblad$Number%in%as.character(ls.coreMd),]$Genus2<-" "
# Make plot
pBar.DaFilfBlad <- ggplot(dt.mdILFblad, aes(x=Sample, y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Family~Da1Fb, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.DaFilfBlad # print plot










### Macrobdella DaF Intestinum only
phyT.mdDaFInt<-subset_samples(phyT.mdDaF,Sample_Type%in%c("Intestinum")) # Keep only Intestinum samples
phyT.mdDaFpr<-prune_taxa(taxa_sums(phyT.mdDaFInt)>.1,phyT.mdDaFInt) # keep taxa with at least 1% of 1 sample
phyT.mdDaFcore<-subset_taxa(phyT.mdDaFpr, taxa_names(phyT.mdDaFpr)%in%c(ls.coreFeedInt))
sample_data(phyT.mdDaFcore)$Da1Fb = factor(sample_data(phyT.mdDaFcore)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,113,215)) # Reorder Da1Fb
# Convert phyloseq to data table
dt.Int<-data.table(psmelt(phyT.mdDaFcore))
dt.Int$Number<-as.character(dt.Int$Number)
dt.Int$Genus<-as.character(dt.Int$Genus2)
dt.Int[!dt.Int$Number%in%as.character(ls.coreTot),]$Genus2<-" "
# Make plot
pBar.DaFint <- ggplot(dt.Int, aes(x=Sample,y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Order~Da1Fb, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_blank(),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(9,"Set1"))
pBar.DaFint # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaFint), size = "last")), filename="Plots/plotBar_DaFInt.eps", device="eps", width=12,height=8)
##### Table #####
write.table(sample_data(phyT.core.mdDaF), "tableTax_Da1F.csv", sep=",")


#### ILF #####
### Box + whisker ###
library(grid)
library(gtable)

phyT.bwdat<-phyT.mdDaFILF 
sample_data(phyT.bwdat)$Da1Fb = factor(sample_data(phyT.bwdat)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1Fb
phyT.bwdatP<-prune_taxa(c(ls.coreMdILF),phyT.bwdat) # keep only core ILF taxa (box+whisker data pruned)
phyT.bwdatPmer<-tax_glom(phyT.bwdatP,"Number")
ls.bwGen<-as.character(get_taxa_unique(phyT.bwdatPmer, "Number")) # compile list of genera to examine (box + whisker genera)

lowA<-1e-3 # set ymin
dt.bwMelt <- psmelt(phyT.bwdatPmer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)

phyT.bwdatP0<-subset_samples(phyT.bwdatP,Da1Fb=="0")
bwMelt0<-psmelt(phyT.bwdatP0)
# add empty columns to data table
dt.bwMelt$Q25<-0
dt.bwMelt$Q50<-0
dt.bwMelt$Q75<-0
bwMelt0$Q25<-0
bwMelt0$Q50<-0
bwMelt0$Q75<-0
bwMelt0$box2<-0
bwMelt0$box3<-0
bwMelt0$box4<-0
bwMelt0$conf1<-0
bwMelt0$conf2<-0

for (f in c(ls.bwGen)){
  bwMelt0[,f]<-with(bwMelt0,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Number) %>%
  summarise(avg_Q25 = mean(Q25),
            avg_Q50 = mean(Q50),
            avg_Q75 = mean(Q75),
            hinge1 = mean(box2),
            avg_box3 = mean(box3),
            hinge2 = mean(box4),
            conf1 = mean(conf1),
            conf2 = mean(conf2))

# Enter Abundance in new taxa-specific columns in data.table (dt.bwMelt). 
for (f in c(ls.bwGen)){
  dt.bwMelt[,f]<-with(dt.bwMelt,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(dt.bwMelt,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (dt.bwMelt). #   
  dt.bwMelt$Q25<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q25,dt.bwMelt$Q25))
  dt.bwMelt$Q50<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q50,dt.bwMelt$Q50))
  dt.bwMelt$Q75<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q75,dt.bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pbox.ProtILF <- ggplot(dt.bwMelt) +
  ggtitle("Macrobdella ILF Taxa") + 
  stat_boxplot(data=dt.bwMelt, aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  scale_y_log10(limits=c(.0001,1)) + 
  facet_grid(~Number+Genus2, scales="free_x",space="free") +
  geom_line(data=dt.bwMelt,aes(x=as.numeric(Da1Fb),y=Q50),linetype=2) + 
  geom_ribbon(data=dt.bwMelt, aes(x=as.numeric(Da1Fb), ymin=Q25, ymax=Q75),fill="black", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=brewer.pal(8,"Greys"))   
pbox.ProtILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbox.ProtILF), size = "last")), filename="Plots/plotBox_Da1FILF.eps", device="eps", width=12,height=8)

#### INTESTINUM #####
### Box + whisker ###
phyT.bwdat<-phyT.mdDaFInt
sample_data(phyT.bwdat)$Da1Fb = factor(sample_data(phyT.bwdat)$Da1Fb, levels = c(0,1,2,4,7,30,"90+",100)) # Reorder Da1Fb
phyT.bwdatP<-prune_taxa(c(ls.coreMdInt),phyT.bwdat) # keep only core ILF taxa (box+whisker data pruned)
phyT.bwdatPmer<-tax_glom(phyT.bwdatP,"Number")

lowA<-1e-3 # set ymin
dt.bwMelt <- psmelt(phyT.bwdatPmer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)
ls.bwGen<-as.character(get_taxa_unique(phyT.bwdatPmer, "Number")) # compile list of genera to examine (box + whisker genera)

phyT.bwdatP0<-subset_samples(phyT.bwdatP,Da1Fb=="0")
bwMelt0<-psmelt(phyT.bwdatP0)

dt.bwMelt$Q25<-0
dt.bwMelt$Q50<-0
dt.bwMelt$Q75<-0
bwMelt0$Q25<-0
bwMelt0$Q50<-0
bwMelt0$Q75<-0
bwMelt0$box2<-0
bwMelt0$box3<-0
bwMelt0$box4<-0
bwMelt0$conf1<-0
bwMelt0$conf2<-0

for (f in c(ls.bwGen)){
  bwMelt0[,f]<-with(bwMelt0,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Number) %>%
  summarise(avg_Q25  = mean(Q25),
            avg_Q50 = mean(Q50),
            avg_Q75 = mean(Q75),
            hinge1 = mean(box2),
            avg_box3 = mean(box3),
            hinge2 = mean(box4),
            conf1 = mean(conf1),
            conf2 = mean(conf2))

# Enter Abundance in new taxa-specific columns in data.table (dt.bwMelt). 
for (f in c(ls.bwGen)){
  dt.bwMelt[,f]<-with(dt.bwMelt,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(dt.bwMelt,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (dt.bwMelt). #   
  dt.bwMelt$Q25<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q25,dt.bwMelt$Q25))
  dt.bwMelt$Q50<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q50,dt.bwMelt$Q50))
  dt.bwMelt$Q75<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q75,dt.bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pbox.ProtInt <- ggplot(dt.bwMelt) +
  ggtitle("Macrobdella Intestinum Taxa") + 
  stat_boxplot(data=dt.bwMelt, aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  scale_y_log10(limits=c(lowA,1)) + 
  facet_grid(~Number+Genus2, scales="free_x",space="free") +
  geom_line(data=dt.bwMelt,aes(x=as.numeric(Da1Fb),y=Q50),linetype=2) + 
  geom_ribbon(data=dt.bwMelt, aes(x=as.numeric(Da1Fb), ymin=Q25, ymax=Q75),fill="black", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=brewer.pal(6,"Greys"))   
pbox.ProtInt
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbox.ProtInt), size = "last")), filename="Plots/plotBox_Da1FInt.eps", device="eps", width=12,height=8)

###############################################
################ Practice area ################ 
############################################### 
#phyT.hvDaFILF<-subset_samples(phyT.hv,Sample_Type%in%c("ILF")) # Keep only ILF samples

#### ILF #####
### Box + whisker ###
library(grid)
library(gtable)

phyT.bwdat<-merge_phyloseq(phyT.mdDaFILF,phyT.hvILF) 
sample_data(phyT.bwdat)$Da1Fb = factor(sample_data(phyT.bwdat)$Da1Fb, levels = c(0,1,2,4,7,30,31,35,"90+",100,113,215)) # Reorder Da1Fb
phyT.bwdatP<-prune_taxa(union(ls.coreMdILF,ls.comHvILF0),phyT.bwdat) # keep only core ILF taxa (box+whisker data pruned)
phyT.bwdatPmer<-tax_glom(phyT.bwdatP,"Genus2")

ls.bwGen<-as.character(get_taxa_unique(phyT.bwdatPmer, "Number")) # compile list of genera to examine (box + whisker genera)

lowA<-1e-3 # set ymin
dt.bwMelt <- psmelt(phyT.bwdatPmer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)

phyT.bwdatP0<-subset_samples(phyT.bwdatP,Da1Fb=="0")
bwMelt0<-psmelt(phyT.bwdatP0)
# add empty columns to data table
dt.bwMelt$Q25<-0
dt.bwMelt$Q50<-0
dt.bwMelt$Q75<-0
bwMelt0$Q25<-0
bwMelt0$Q50<-0
bwMelt0$Q75<-0
bwMelt0$box2<-0
bwMelt0$box3<-0
bwMelt0$box4<-0
bwMelt0$conf1<-0
bwMelt0$conf2<-0

for (f in c(ls.bwGen)){
  bwMelt0[,f]<-with(bwMelt0,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Number==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Number==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Number) %>%
  summarise(avg_Q25 = mean(Q25),
            avg_Q50 = mean(Q50),
            avg_Q75 = mean(Q75),
            hinge1 = mean(box2),
            avg_box3 = mean(box3),
            hinge2 = mean(box4),
            conf1 = mean(conf1),
            conf2 = mean(conf2))

# Enter Abundance in new taxa-specific columns in data.table (dt.bwMelt). 
for (f in c(ls.bwGen)){
  dt.bwMelt[,f]<-with(dt.bwMelt,ifelse(Number==eval(f),Abundance,0)) 
  abundPro<-subset(dt.bwMelt,Number==eval(f) | Da1Fb==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (dt.bwMelt). #   
  dt.bwMelt$Q25<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q25,dt.bwMelt$Q25))
  dt.bwMelt$Q50<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q50,dt.bwMelt$Q50))
  dt.bwMelt$Q75<-with(dt.bwMelt,ifelse(Number==eval(f),filter(sumA,Number==eval(f))$avg_Q75,dt.bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pbox.ProtILF <- ggplot(dt.bwMelt) +
  ggtitle("Macrobdella ILF Taxa") + 
  stat_boxplot(data=dt.bwMelt, aes(x=Da1Fb, y=Abundance, fill=Da1Fb)) +    
  scale_y_log10(limits=c(.0001,1)) + 
  facet_grid(Taxonomic_ID~Genus2, scales="free_x",space="free") +
  #geom_line(data=dt.bwMelt,aes(x=as.numeric(Da1Fb),y=Q50),linetype=2) + 
  #geom_ribbon(data=dt.bwMelt, aes(x=as.numeric(Da1Fb), ymin=Q25, ymax=Q75),fill="black", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=brewer.pal(8,"Greys"))   
pbox.ProtILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbox.ProtILF), size = "last")), filename="Plots/plotBox_Da1FILF.eps", device="eps", width=12,height=8)

phyT.hvILF
phyTmat.hvILF <- prune_taxa(names(sort(taxa_sums(phyT.hvILF), TRUE)[1:35]),phyT.hvILF) # Create a subset of data including only 27 most abundant taxa
pBar.hvILF <- plot_bar(phyTmat.hvILF,x="Sample",fill="BLASTn") + 
  facet_grid(Sample_Type~Da1F, scales="free_x",space="free") + 
  scale_fill_manual(values=pairBiome)
pBar.hvILF

phyT.aaa<-subset_samples(phyT.hvILF,Da1Fb%in%c("0","90+"))
phyT.hvAer<-prune_taxa(c("denovo211217","denovo81947","denovo102614","denovo174765","denovo188726"),phyT.hvILF)
min(sample_sums(phyT.hvAer))
median(sample_sums(phyT.hvAer))
max(sample_sums(phyT.hvAer))
phyT.hvMuc<-prune_taxa(c("denovo235669","denovo931"),phyT.hvILF)
min(sample_sums(phyT.hvMuc))
median(sample_sums(phyT.hvMuc))
max(sample_sums(phyT.hvMuc))
phyT.hvDom<-prune_taxa(c("denovo211217","denovo81947","denovo102614","denovo174765","denovo188726","denovo235669","denovo931"),phyT.hvILF)
min(sample_sums(phyT.hvDom))
median(sample_sums(phyT.hvDom))
max(sample_sums(phyT.hvDom))
phyT.hvCl<-prune_taxa(c("denovo136656","denovo233988","denovo44870","denovo247165"),phyT.hvILF)
min(sample_sums(phyT.hvCl))
median(sample_sums(phyT.hvCl))
max(sample_sums(phyT.hvCl))

max(sample_sums(prune_taxa(c("denovo172596","denovo80607"),phyT.hvILF)))

dt.hvILF = fast_melt(phyT.hvILF) # make data table from physeq object (data table. physeq Pruned)
prev.hvILF = dt.hvILF[, list(Prevalence = sum(count >= .001), 
                               TotalPer = sum(count),
                               MinCount = min(count),
                               MaxCount = max(count)),
                        by = TaxaID] # make simple table listing 'TaxaID, Prevalence, and TotalPer' (prevalence . physeq Pruned)







###############################################
################ Practice area ################ 
############################################### 
# Thank you to jeffkimbrel ! #
##### weighted Unifrac distance by Da1F #####
##### Unifrac distance by Da1F #####
phyT.aaa<-subset_samples(phyT.start,Sample_Type=="ILF" & Taxonomic_ID=="Mdecora")
sample_data(phyT.aaa)$fuck<-sample_names(phyT.aaa)

dist.md<- phyloseq::distance(phyT.aaa, method = "wunifrac") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)

# remove self-comparisons
md = melt(as.matrix(dist.md)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(phyT.aaa) %>%
  select("fuck","Da1Fb") %>%
  mutate_if(is.factor,as.character) 

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(md, sd, by = "Var1")

colnames(sd) = c("Var2", "Da1F")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd2<-wu.sd[wu.sd$Type1==0,]
wu.sd2$Da1F = factor(wu.sd2$Da1F, levels = c(0,1,2,4,7,30,"90+")) # Reorder Da1Fb
colnames(wu.sd2)[colnames(wu.sd2)=="value"] <- "distance"

# plot
ggplot(wu.sd2, aes(x = Da1F, y = distance)) +
  theme_bw() +
  geom_point() +
  geom_jitter(width = 0.2) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Macrobdella decora Weighted Unifrac")


##### H verbana #####
##### Unifrac distance by Da1F #####
phyT.aaa<-subset_samples(phyT.start,Sample_Type=="ILF" & Taxonomic_ID=="Hverbana")
sample_data(phyT.aaa)$fuck<-sample_names(phyT.aaa)

dist.hv<- phyloseq::distance(phyT.aaa, method = "wunifrac") # Calculate bray curtis distance matrix (bray metric . macrobdella decora)

# remove self-comparisons
md = melt(as.matrix(dist.hv)) %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = sample_data(phyT.aaa) %>%
  select("fuck","Da1Fb") %>%
  mutate_if(is.factor,as.character) %>%
  mutate_if(is.integer,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(md, sd, by = "Var1")

colnames(sd) = c("Var2", "Da1F")
wu.sd = left_join(wu.sd, sd, by = "Var2")

wu.sd2<-wu.sd[wu.sd$Type1==0,]
wu.sd2$Da1F = factor(wu.sd2$Da1F, levels = c(0,1,2,4,7,30,"90+")) # Reorder Da1Fb
colnames(wu.sd2)[colnames(wu.sd2)=="value"] <- "distance"

# plot
ggplot(wu.sd2, aes(x = Da1F, y = distance)) +
  theme_bw() +
  geom_point() +
  geom_jitter(width = 0.2) +
  geom_boxplot() +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Hirudo verbana Weighted Unifrac")

