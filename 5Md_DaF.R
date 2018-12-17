########## Macrobdella Da1F plot with 'Other' category ##########
phyT.mdDaF<-merge_phyloseq(phyT.mdCT,phyT.mdMA) # merge fed samples from CT and MA
# combine 31 DaF and 35 DaF into one 30 DaF group
sdat.mdDaF<-sample_data(phyT.mdDaF)
sdat.mdDaF$Da1F<-with(sdat.mdDaF,
  ifelse(Da1F=="31","30",
  ifelse(Da1F=="35","30",
  ifelse(Da1F=="100","100",
  ifelse(Da1F=="113","100",
  ifelse(Da1F=="215","100",
  as.character(Da1F)))))))  
phyT.mdDaF<-phyloseq(otu_table(phyT.mdDaF),tax_table(phyT.mdDaF),sdat.mdDaF,phy_tree(phyT.mdDaF))
phyT.mdDaFGI<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF","Intestinum")) # Keep only ILF and intestinum samples
phyT.core.mdDaF<-prune_taxa(ls.coreTot,phyT.mdDaFGI) # keep only core taxa
sample_data(phyT.mdDaFGI)$Da1F = factor(sample_data(phyT.mdDaFGI)$Da1F, levels = c(0,1,2,4,7,30,31,35,100,113,215)) # Reorder Da1F

### Make plot
dt.DaF<-data.table(psmelt(phyT.mdDaFGI))
dt.DaF$Number<-as.character(dt.DaF$Number)
dt.DaF$Genus3<-with(dt.DaF,
                     ifelse(Number%in%c(ls.coreTot),as.character(Genus2),
                     as.character("other")))
# Define levels for OTUs
dt.DaF$Genus3 <- factor(dt.DaF$Genus3, levels=c("Aeromonas","Aeromonas2",
                                                  "Mucinivorans","Bacteroides","unk_Rikenellaceae","unk_Bacteroides",
                                                  "unk_Peptostreptococcaceae","Proteocatella","unk_Proteocatella","unk_Butyricicoccus","unk_Ruminococcaceae",
                                                  "Ochrobactrum","unk_Rhodospirillaceae","Rhizobium","unk_Rhizobiales",
                                                  "unk_Comamonadaceae","unk_Methylophilaceae",
                                                  "Bdellovibrio","Desulfovibrio","unk_Myxococcales",
                                                  "Nubsella","Cetobacterium","Fusobacterium",
                                                  "other"))
# Define colors for OTUs. Reference = "rainbow"
genus.color<-c(Aeromonas="#267326",Aeromonas2="#39ac39",
               Mucinivorans="#B80000",Bacteroides="#F00000",unk_Rikenellaceae="#FF7777",unk_Bacteroides="#ffcccc",
               unk_Peptostreptococcaceae="#000080",Proteocatella="#0000cd",unk_Proteocatella="#8282ff",unk_Butyricicoccus="#cfcfff",unk_Ruminococcaceae="#0090b4",
               Ochrobactrum="#b34700",unk_Rhodospirillaceae="#ff6600",Rhizobium="#ff9933",unk_Rhizobiales="#ffcc99",
               unk_Comamonadaceae="#600080",unk_Methylophilaceae="#ac00e6",
               Bdellovibrio="#cccc00",Desulfovibrio="#ffff00",unk_Myxococcales="#ffffb3",
               Nubsella="#9a0066",Cetobacterium="#ffe6f7",Fusobacterium="#ff4ec5",
               other="#808080") 

pBar.DaF <- ggplot(dt.DaF, aes(x=Sample,y=Abundance, fill=Genus3)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~AnimalSource+Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=genus.color)
pBar.DaF # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaF), size = "last")), filename="Plots/plotBarStack_Da1F.png", width=12,height=8)


### Macrobdella DaF ILF only
phyT.mdDaFILF<-subset_samples(phyT.mdDaF,Sample_Type%in%c("ILF")) # Keep only ILF samples
phyT.mdDaFgr<-subset_samples(phyT.mdDaFILF,AnimalSource=="GrotonMA")
phyT.mdILFpr<-prune_taxa(taxa_sums(phyT.mdDaFgr)>.001,phyT.mdDaFgr) # keep taxa with at least 1% of 1 sample
phyT.mdILFfam<-subset_taxa(phyT.mdILFpr, taxa_names(phyT.mdILFpr)%in%c(ls.coreMd))
sample_data(phyT.mdILFfam)$Da1F = factor(sample_data(phyT.mdILFfam)$Da1F, levels = c(0,1,4,7,30,31,35,113,215)) # Reorder Da1F
# Convert phyloseq to data table
dt.ILF<-data.table(psmelt(phyT.mdILFfam))
dt.ILF$Number0o<-as.character(dt.ILF$Number)
dt.ILF$Genus2<-as.character(dt.ILF$Genus2)
dt.ILF[!dt.ILF$Number%in%as.character(ls.coreMd),]$Genus2<-" "
# Make plot
pBar.DaFilf <- ggplot(dt.ILF, aes(x=Sample, y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Family~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.DaFilf # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaFilf), size = "last")), filename="Plots/plotBar_DaFILF.png", width=12,height=8)



phyT.mdILFblad<-subset_taxa(phyT.mdILFpr, taxa_names(phyT.mdILFpr)%in%c(setdiff(ls.coreMdBlad,ls.coreMdILF)))

sample_data(phyT.mdILFblad)$Da1F = factor(sample_data(phyT.mdILFblad)$Da1F, levels = c(0,1,4,7,30,31,35,113,215)) # Reorder Da1F

# Convert phyloseq to data table
dt.mdILFblad<-data.table(psmelt(phyT.mdILFblad))
dt.mdILFblad$Number0o<-as.character(dt.mdILFblad$Number)
dt.mdILFblad$Genus2<-as.character(dt.mdILFblad$Genus2)
dt.mdILFblad[!dt.mdILFblad$Number%in%as.character(ls.coreMd),]$Genus2<-" "
# Make plot
pBar.DaFilfBlad <- ggplot(dt.mdILFblad, aes(x=Sample, y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Family~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pBar.DaFilfBlad # print plot










### Macrobdella DaF Intestinum only
phyT.mdDaFInt<-subset_samples(phyT.mdDaF,Sample_Type%in%c("Intestinum")) # Keep only Intestinum samples
phyT.mdDaFpr<-prune_taxa(taxa_sums(phyT.mdDaFInt)>.01,phyT.mdDaFInt) # keep taxa with at least 1% of 1 sample
phyT.mdDaFcore<-subset_taxa(phyT.mdDaFpr, taxa_names(phyT.mdDaFpr)%in%c(ls.coreFeedInt))
sample_data(phyT.mdDaFcore)$Da1F = factor(sample_data(phyT.mdDaFcore)$Da1F, levels = c(0,1,2,4,7,30,31,35,113,215)) # Reorder Da1F
# Convert phyloseq to data table
dt.Int<-data.table(psmelt(phyT.mdDaFcore))
dt.Int$Number<-as.character(dt.Int$Number)
dt.Int$Genus<-as.character(dt.Int$Genus2)
dt.Int[!dt.Int$Number%in%as.character(ls.coreTot),]$Genus2<-" "
# Make plot
pBar.DaFint <- ggplot(dt.Int, aes(x=Sample,y=Abundance, fill=Genus2)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Family~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10),axis.text.x=element_text(angle=90),axis.title.x=element_blank()) +
  scale_fill_manual(values=brewer.pal(9,"Set1"))
pBar.DaFint # print plot
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pBar.DaFint), size = "last")), filename="Plots/plotBar_DaFInt.png", width=12,height=8)
##### Table #####
write.table(sample_data(phyT.core.mdDaF), "tableTax_Da1F.csv", sep=",")


#### ILF #####
### Box + whisker ###
library(grid)
library(gtable)

phyT.bwdat<-phyT.mdDaFILF 
sample_data(phyT.bwdat)$Da1F = factor(sample_data(phyT.bwdat)$Da1F, levels = c(0,1,4,7,30,100)) # Reorder Da1F
phyT.bwdatP<-prune_taxa(c(ls.coreMdILF),phyT.bwdat) # keep only core ILF taxa (box+whisker data pruned)
phyT.bwdatPmer<-tax_glom(phyT.bwdatP,"Genus2")
ls.bwGen<-as.character(get_taxa_unique(phyT.bwdatPmer, "Genus2")) # compile list of genera to examine (box + whisker genera)

lowA<-1e-3 # set ymin
dt.bwMelt <- psmelt(phyT.bwdatPmer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)

phyT.bwdatP0<-subset_samples(phyT.bwdatP,Da1F=="0")
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
  bwMelt0[,f]<-with(bwMelt0,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Genus2==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Genus2) %>%
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
  dt.bwMelt[,f]<-with(dt.bwMelt,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(dt.bwMelt,Genus2==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (dt.bwMelt). #   
  dt.bwMelt$Q25<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q25,dt.bwMelt$Q25))
  dt.bwMelt$Q50<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q50,dt.bwMelt$Q50))
  dt.bwMelt$Q75<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q75,dt.bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pbox.ProtILF <- ggplot(dt.bwMelt) +
  ggtitle("Macrobdella ILF Taxa") + 
  stat_boxplot(data=dt.bwMelt, aes(x=Da1F, y=Abundance, fill=Da1F)) +    
  scale_y_log10(limits=c(.0001,1)) + 
  facet_grid(~Genus2, scales="free_x",space="free") +
  geom_line(data=dt.bwMelt,aes(x=as.numeric(Da1F),y=Q50),linetype=2) + 
  geom_ribbon(data=dt.bwMelt, aes(x=as.numeric(Da1F), ymin=Q25, ymax=Q75),fill="red", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=rainbow)   
pbox.ProtILF
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbox.ProtILF), size = "last")), filename="Plots/plotBox_Da1FILF.png", width=12,height=8)

#### INTESTINUM #####
### Box + whisker ###
phyT.bwdat<-phyT.mdDaFInt
sample_data(phyT.bwdat)$Da1F = factor(sample_data(phyT.bwdat)$Da1F, levels = c(0,1,4,7,30,100)) # Reorder Da1F
phyT.bwdatP<-prune_taxa(c(ls.coreMdInt),phyT.bwdat) # keep only core ILF taxa (box+whisker data pruned)
phyT.bwdatPmer<-tax_glom(phyT.bwdatP,"Genus2")

lowA<-1e-3 # set ymin
dt.bwMelt <- psmelt(phyT.bwdatPmer) # create data.table from phyloseq object, phyT.bwdatFam (box+whisker melt)
ls.bwGen<-as.character(get_taxa_unique(phyT.bwdatPmer, "Genus2")) # compile list of genera to examine (box + whisker genera)

phyT.bwdatP0<-subset_samples(phyT.bwdatP,Da1F=="0")
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
  bwMelt0[,f]<-with(bwMelt0,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Genus2==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Genus2==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Genus2==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Genus2) %>%
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
  dt.bwMelt[,f]<-with(dt.bwMelt,ifelse(Genus2==eval(f),Abundance,0)) 
  abundPro<-subset(dt.bwMelt,Genus2==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (dt.bwMelt). #   
  dt.bwMelt$Q25<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q25,dt.bwMelt$Q25))
  dt.bwMelt$Q50<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q50,dt.bwMelt$Q50))
  dt.bwMelt$Q75<-with(dt.bwMelt,ifelse(Genus2==eval(f),filter(sumA,Genus2==eval(f))$avg_Q75,dt.bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pbox.ProtInt <- ggplot(dt.bwMelt) +
  ggtitle("Macrobdella Intestinum Taxa") + 
  stat_boxplot(data=dt.bwMelt, aes(x=Da1F, y=Abundance, fill=Da1F)) +    
  scale_y_log10(limits=c(lowA,1)) + 
  facet_grid(~Genus2, scales="free_x",space="free") +
  geom_line(data=dt.bwMelt,aes(x=as.numeric(Da1F),y=Q50),linetype=2) + 
  geom_ribbon(data=dt.bwMelt, aes(x=as.numeric(Da1F), ymin=Q25, ymax=Q75),fill="red", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(1)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=rainbow)   
pbox.ProtInt
##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pbox.ProtInt), size = "last")), filename="Plots/plotBox_Da1FInt.png", width=12,height=8)





###############################################
################ Practice area ################ 
############################################### 


