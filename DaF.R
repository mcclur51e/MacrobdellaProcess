########## Macrobdella Da1F plot with 'Other' category ##########
mac1w<-subset_samples(ctMacro,Da1F%in%c("0","1","2","4","7","31","113","215")) # keep CT samples that were fed blood
mac5w<-subset_samples(maMacro,Da1F%in%c("35")) # keep MA samples that were held for 5w
macDaf<-merge_phyloseq(mac1w,mac5w) # merge fed samples from CT and MA
mdDaF<-sample_data(macDaf)
# combine 31 DaF and 35 DaF into one 30 DaF group
mdDaF$Da1F<-with(mdDaF,
  ifelse(Da1F=="31","30",
  ifelse(Da1F=="35","30",
  as.character(Da1F))))  
macDaf<-phyloseq(otu_table(macDaf),tax_table(macDaf),mdDaF,phy_tree(macDaf))

macDafI<-subset_samples(macDaf,Sample_Type%in%c("ILF","Intestinum")) # Keep only ILF and intestinum samples
macDaFpr<-prune_taxa(taxa_sums(macDafI)>.01,macDafI) # keep taxa with at least 1% of 1 sample
sample_data(macDaFpr)$Da1F = factor(sample_data(macDaFpr)$Da1F, levels = c(0,1,2,4,7,30,113,215)) # Reorder Da1F

### Make plot
dtDaF<-data.table(psmelt(macDaFpr))
dtDaF$Number<-as.character(dtDaF$Number)
dtDaF$Genus<-as.character(dtDaF$Genus)
dtDaF[!dtDaF$Number%in%as.character(coreTot),]$Genus<-" "

pDaF <- ggplot(dtDaF, aes(x=Replicate, y=Abundance, fill=Order)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Sample_Type~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pDaF # print plot

### Macrobdella DaF ILF only
macDafILF<-subset_samples(macDaf,Sample_Type%in%c("ILF")) # Keep only ILF samples
macILFpr<-prune_taxa(taxa_sums(macDafILF)>.01,macDafILF) # keep taxa with at least 1% of 1 sample
macILFfam<-subset_taxa(macILFpr, Family%in%c("Aeromonadaceae","Bacteroidaceae","Ruminococcaceae","unk_Clostridiales"))
sample_data(macILFfam)$Da1F = factor(sample_data(macILFfam)$Da1F, levels = c(0,1,2,4,7,14,28,30,113,215)) # Reorder Da1F
# Convert phyloseq to data table
dtILF<-data.table(psmelt(macILFfam))
dtILF$Number0o<-as.character(dtILF$Number)
dtILF$Genus<-as.character(dtILF$Genus)
dtILF[!dtILF$Number%in%as.character(coreTot),]$Genus<-" "
# Make plot
pDaFilf <- ggplot(dtILF, aes(x=Replicate, y=Abundance, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Order~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pDaFilf # print plot

### Macrobdella DaF Intestinum only
macDafInt<-subset_samples(macDaf,Sample_Type%in%c("Intestinum")) # Keep only Intestinum samples
macIntpr<-prune_taxa(taxa_sums(macDafInt)>.01,macDafInt) # keep taxa with at least 1% of 1 sample
macIntFam<-subset_taxa(macIntpr, Family%in%c("Aeromonadaceae","Bacteroidaceae","Ruminococcaceae","unk_Clostridiales"))
sample_data(macIntFam)$Da1F = factor(sample_data(macIntFam)$Da1F, levels = c(0,1,2,4,7,14,28,30,113,215)) # Reorder Da1F
# Convert phyloseq to data table
dtInt<-data.table(psmelt(macIntFam))
dtInt$Number<-as.character(dtInt$Number)
dtInt$Genus<-as.character(dtInt$Genus)
dtInt[!dtInt$Number%in%as.character(coreTot),]$Genus<-" "
# Make plot
pDaFint <- ggplot(dtInt, aes(x=Replicate, y=Abundance, fill=Genus)) + 
  geom_bar(aes(), stat="identity", position="stack") +
  facet_grid(Order~Da1F, scales="free_x",space="free") +
  theme(text=element_text(size=10), axis.title.x=element_blank()) +
  scale_fill_manual(values=pairBiome)
pDaFint # print plot

##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pDaF), size = "last")), filename="plotDaF.png", width=12,height=8)
##### Figure #####

##### Table #####
write.table(sample_data(macDaFpr), "taxTable.csv", sep=",")
##### Table #####



### Box + whisker ###
library(grid)
library(gtable)

bwdat<-macILFfam

sample_data(macIntFam)$Da1F = factor(sample_data(macIntFam)$Da1F, levels = c(0,1,2,4,7,14,28,30,113,215)) # Reorder Da1F
bwMAT <- sort(taxa_sums(bwdat), TRUE)[1:5] # Identify 5 most abundant taxa
bwdatP<-prune_taxa(names(bwMAT),bwdat) # keep only most abundant taxa (box+whisker data pruned)
bwdatPmer<-merge_taxa(bwdatP,"Genus")

lowA<-1e-3 # set ymin
bwMelt <- psmelt(bwdatPmer) # create data.table from phyloseq object, bwdatFam (box+whisker melt)
bwGen<-as.character(get_taxa_unique(bwdatPmer, "Genus")) # compile list of genera to examine (box + whisker genera)

bwdatP0<-subset_samples(bwdatP,Da1F=="0")
bwMelt0<-psmelt(bwdatP0)

bwMelt$Q25<-0
bwMelt$Q50<-0
bwMelt$Q75<-0
bwMelt0$Q25<-0
bwMelt0$Q50<-0
bwMelt0$Q75<-0
bwMelt0$box2<-0
bwMelt0$box3<-0
bwMelt0$box4<-0
bwMelt0$conf1<-0
bwMelt0$conf2<-0

for (f in c(bwGen)){
  bwMelt0[,f]<-with(bwMelt0,ifelse(Genus==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt0,Genus==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt0). #   
  bwMelt0$conf1<-with(bwMelt0,ifelse(Genus==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[1],bwMelt0$conf1))
  bwMelt0$conf2<-with(bwMelt0,ifelse(Genus==eval(f),boxplot.stats(abundPro[abundPro>0])$conf[2],bwMelt0$conf2))
  bwMelt0$box2<-with(bwMelt0,ifelse(Genus==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[2],bwMelt0$box2))
  bwMelt0$box3<-with(bwMelt0,ifelse(Genus==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[3],bwMelt0$box3))
  bwMelt0$box4<-with(bwMelt0,ifelse(Genus==eval(f),boxplot.stats(abundPro[abundPro>0])$stats[4],bwMelt0$box4))
  bwMelt0$Q25<-with(bwMelt0,ifelse(Genus==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.25)),bwMelt0$Q25))
  bwMelt0$Q50<-with(bwMelt0,ifelse(Genus==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.50)),bwMelt0$Q50))
  bwMelt0$Q75<-with(bwMelt0,ifelse(Genus==eval(f),quantile(abundPro[abundPro>0.001],probs=c(.75)),bwMelt0$Q75))
}

# Compile a table that is the 'average' of Q values from bwMelt0 by Genus
sumA<-bwMelt0 %>% 
  group_by(Genus) %>%
  summarise(avg_Q25 = mean(Q25),
            avg_Q50 = mean(Q50),
            avg_Q75 = mean(Q75),
            hinge1 = mean(box2),
            avg_box3 = mean(box3),
            hinge2 = mean(box4),
            conf1 = mean(conf1),
            conf2 = mean(conf2))

# Enter Abundance in new taxa-specific columns in data.table (bwMelt). 
for (f in c(bwGen)){
  bwMelt[,f]<-with(bwMelt,ifelse(Genus==eval(f),Abundance,0)) 
  abundPro<-subset(bwMelt,Genus==eval(f) | Da1F==0)[,eval(f)]
  # Calculate max, min, mean, median, 25% quartile, 50% quartile, and 75% quartile of taxa-specific Abundance and enter into data.table (bwMelt). #   
  bwMelt$Q25<-with(bwMelt,ifelse(Genus==eval(f),filter(sumA,Genus==eval(f))$avg_Q25,bwMelt$Q25))
  bwMelt$Q50<-with(bwMelt,ifelse(Genus==eval(f),filter(sumA,Genus==eval(f))$avg_Q50,bwMelt$Q50))
  bwMelt$Q75<-with(bwMelt,ifelse(Genus==eval(f),filter(sumA,Genus==eval(f))$avg_Q75,bwMelt$Q75))
}

# VERTICAL boxplot + log y-scale + facet by age + hi/lo bounds
pProt <- ggplot(bwMelt) +
  ggtitle("Macrobdella ILF Taxa") + 
  stat_boxplot(data=bwMelt, aes(x=Da1F, y=Abundance, fill=Da1F)) +    
  scale_y_log10(limits=c(lowA,1)) + 
  facet_grid(~Genus, scales="free_x",space="free") +
  geom_line(data=bwMelt,aes(x=as.numeric(Da1F),y=Q50),linetype=2) + 
  geom_ribbon(data=bwMelt, aes(x=as.numeric(Da1F), ymin=Q25, ymax=Q75),fill="red", alpha=0.1) +
  theme(text=element_text(size=10),strip.text=element_text(size=rel(.5)), axis.title.x=element_blank(), legend.position="none") +
  scale_fill_manual(values=rainbow)   
pProt

##### Figure #####
ggsave(grid.draw(rbind(ggplotGrob(pProt), size = "last")), filename="plotDaFbw.png", width=12,height=8)
##### Figure #####

##### Table #####
write.table(bwMelt0, "meltTable.csv", sep=",")
##### Table #####








