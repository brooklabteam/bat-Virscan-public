rm(list = ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(ggthemes)
library(tidyverse)

#set wd - make a link to homewd for your own folder
homewd = "/Users/carabrook/Developer/bat-VirScan-public"
setwd(homewd)

#load data from Phip-Seq
dat <- read.csv(file = paste0(homewd, "/working-data/merge-data.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)
unique(dat$ID)
dat$ID[dat$ID=="Pa6 - repeat"] <- "Pa6"
length(unique(dat$ID)) #82 unique bats
unique(dat$Species)

#load and merge subfamily data 
subfam <- read.csv(file = paste0(homewd, "/working-data/subfamily-merge.csv"), header = T, stringsAsFactors = F)
head(subfam)

dat <- merge(dat, subfam, by = "virus_species")
head(dat)

length(unique(dat$ID[dat$Species=="Pteropus alecto"])) #77 P. alecto. Some are missing ()
length(unique(dat$ID[dat$Species=="Eonycteris spelaea"])) #5

#load the metadata here to check that the merge was correct
meta <- read.csv(file = paste0(homewd, "/working-data/bat_metadata.csv"), header = T, stringsAsFactors = F)
head(meta)
length(unique(meta$ID)) #97 bats 
#remove those for which no phip-seq id is paired
unique(meta$Run.ID)
unique(meta$ID)
meta$ID[meta$ID=="Pa6 - repeat"] <- "Pa6"
meta$ID[meta$ID=="Pa6 "] <- "Pa6"
meta$ID[meta$ID=="Pa batch 1"] <- "Pa1"

meta <- meta[meta$Run.ID!="",] #95 bats. so there are 11 bats (two Pa1s) for which we are lacking phip-seq data


setdiff(unique(dat$ID), unique(meta$ID)) #none. all phipseq outputs are represented in the metadata
setdiff(unique(meta$ID), unique(dat$ID)) 

dat <- arrange(dat, virus_subfamily, virus_genus, virus_species)
dat$virus_species <- factor(dat$virus_species, levels=unique(dat$virus_species))
head(dat)
 

# we are missing phip-seq data for these bats
#"Pa29"  "Pa50"  "Pa59"  "Pa77"  "Pa81"  "Pa10"  "Pa100" "Pa101" "Pa1"   "Pa2"   "Pa3"  

unique(dat$X) # these only have values on the eonycteris samples, where the columns appear to be shifted
# fix down below


#first, for figure 1, plot a heatmap and compare broad viral exposures for the two species

es.dat = subset(dat, Species=="Eonycteris spelaea") # why are the assigned counts so few here?
head(es.dat) # it appears that the column names are shifted to the right on these guys about halfway through
names(es.dat)[17:22] <- c("X", "Assigned.counts", "Assigned.peptides", "Total.sample.hits", "Total.filtered.sample.hits", "virus_subfamily")
es.dat <- dplyr::select(es.dat, -(X))
head(es.dat)
pa.dat = subset(dat, Species=="Pteropus alecto")
head(pa.dat)
pa.dat <- dplyr::select(pa.dat, -(X))


#first compare the big heatmap across all viruses
pS1 <- ggplot(data=pa.dat)+geom_tile(aes(x=ID, y=virus_species, fill=Assigned.counts)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(), legend.position = c(.95,.1),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family, scales = "free_y")  + ggtitle(bquote(italic("Pteropus alecto")~"virus exposures"))

#and save to supplementary figures


ggsave(file = paste0(homewd,"/supp-figures/figS1.png"),
       plot=pS1,
       units="mm",  
       width=300, 
       height=200, 
       scale=3, 
       dpi=300)


es.dat$Assigned.counts <- as.numeric(es.dat$Assigned.counts)


pS2 <- ggplot(data=es.dat)+geom_tile(aes(x=ID, y=virus_species, fill=Assigned.counts)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(),# legend.position = c(.95,.1),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family, scales = "free_y") + ggtitle(bquote(italic("Eonycteris spelae")~"virus exposures"))

ggsave(file = paste0(homewd,"/supp-figures/figS2.png"),
       plot=pS2,
       units="mm",  
       width=150, 
       height=70, 
       scale=3, 
       dpi=300)


#try at the genus level???

#first compare the big heatmap across all viruses
pS1b <- ggplot(data=pa.dat)+geom_tile(aes(x=ID, y=virus_genus, fill=Assigned.counts)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(), legend.position = c(.95,.1),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family, scales = "free_y")  + ggtitle(bquote(italic("Pteropus alecto")~"virus exposures"))

#and save to supplementary figures


ggsave(file = paste0(homewd,"/supp-figures/figS1b.png"),
       plot=pS1b,
       units="mm",  
       width=300, 
       height=200, 
       scale=3, 
       dpi=300)

es.dat$Assigned.counts <- as.numeric(es.dat$Assigned.counts)
pS2b <- ggplot(data=es.dat)+geom_tile(aes(x=ID, y=virus_genus, fill=Assigned.counts)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(),# legend.position = c(.95,.1),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family, scales = "free_y") + ggtitle(bquote(italic("Eonycteris spelae")~"virus exposures"))

ggsave(file = paste0(homewd,"/supp-figures/figS2b.png"),
       plot=pS2b,
       units="mm",  
       width=150, 
       height=70, 
       scale=3, 
       dpi=300)


# and summarize the data to report in the paper
length(unique(dat$virus_family[dat$Species=="Pteropus alecto"])) #33
length(unique(dat$virus_subfamily[dat$Species=="Pteropus alecto"])) #47
length(unique(dat$virus_genus[dat$Species=="Pteropus alecto"])) #82
length(unique(dat$virus_species[dat$Species=="Pteropus alecto"])) #414


length(unique(dat$virus_family[dat$Species=="Eonycteris spelaea"])) #12
length(unique(dat$virus_subfamily[dat$Species=="Eonycteris spelaea"])) #15
length(unique(dat$virus_genus[dat$Species=="Eonycteris spelaea"])) #23
length(unique(dat$virus_species[dat$Species=="Eonycteris spelaea"])) #78

#### subset out the paramyxovirus data and get counts by viral species 
### and individual to look at co-infection/cross-reactivity ####
### make a heat map of hits

pa.dat1 <- subset(pa.dat, virus_family=="Paramyxoviridae")
#pa.dat1 <- subset(pa.dat, virus_family=="Filoviridae")
#pa.dat1 <- subset(pa.dat, virus_family=="Rhabdoviridae")
#pa.dat1 <- subset(pa.dat, virus_family=="Coronaviridae")
es.dat1 <- subset(es.dat, virus_family=="Paramyxoviridae")
#es.dat1 <- subset(es.dat, virus_family=="Rhabdoviridae")
#es.dat1 <- subset(es.dat, virus_family=="Filoviridae")
#es.dat1 <- subset(es.dat, virus_family=="Coronaviridae")

unique(pa.dat1$virus_subfamily)
unique(es.dat1$virus_subfamily)
pa.dat1$virus_subfamily <- factor(pa.dat1$virus_subfamily, levels = c("Avulavirinae", "Rubulavirinae", "Orthoparamyxovirinae"))
es.dat1$virus_subfamily <- factor(es.dat1$virus_subfamily, levels = c("Avulavirinae",  "Orthoparamyxovirinae"))

Fig1a <- ggplot(data=pa.dat1)+geom_tile(aes(x=ID, y=virus_species, fill=Assigned.counts, group=virus_genus)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank(),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family)  + ggtitle(bquote(italic("Pteropus alecto")~"paramyxoviridae exposures"))


Fig1b <- ggplot(data=es.dat1)+geom_tile(aes(x=ID, y=virus_genus, fill=Assigned.counts)) +
  scale_fill_continuous(low="light blue", high="navy", name="peptide counts") + theme_bw() + 
  theme(axis.text.y=element_text(size=7), axis.text.x = element_blank(), 
        strip.background = element_blank(), strip.text = element_blank(),
        axis.ticks.x=element_blank(), axis.title.x = element_blank()) + ylab("virus species\n")+
  facet_wrap(~virus_family)  + ggtitle(bquote(italic("Pteropus alecto")~"Paramyxoviridae exposures"))


#or, with clustering -- need to reorder data
library(pheatmap)
library(reshape2)
head(pa.dat1)
pa.slim <- dplyr::select(pa.dat1, ID, virus_species, Assigned.counts)

pa.slim <- dcast(melt(pa.slim), formula = virus_species~ID)


pa.mat <- as.matrix(pa.slim[2:ncol(pa.slim)])
rownames(pa.mat) <- pa.slim[,1]

#make all NAs = 0
pa.mat[is.na(pa.mat)]<-0

head(subfam)
#also attach genus
merge.dat <- dplyr::select(dat, virus_species, virus_genus)
merge.dat <- merge.dat[!duplicated(merge.dat),]
subfam <- merge(subfam, merge.dat, by="virus_species", all.x = T)

#and subset the families
pa.merge <- dplyr::select(pa.dat1, virus_species)
pa.merge <- cbind.data.frame(virus_species=pa.merge[!duplicated(pa.merge),])
subclustslim <- merge(subfam, pa.merge, by="virus_species")

subclustplot <- dplyr::select(subclustslim, virus_genus, virus_subfamily)
rownames(subclustplot) <- subclustslim$virus_species

#and choose colors
color.ramp <- viridis::inferno(10)
color.ramp <- viridis::cividis(10)
  

  
pheatmap(pa.mat, 
         color = color.ramp,
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_row = subclustplot)

#and manually add gaps as needed
names(subclustplot) <- c("genus", "subfamily")

mat.colors = list(subfamily = c('Avulavirinae' = 'cornflowerblue', 'Orthoparamyxovirinae' = 'tomato', 'Rubulavirinae' = 'mediumseagreen'),
                  genus = c('Henipavirus' = 'red3', 'Respirovirus' = 'goldenrod1', 'unclassified_Paramyxoviridae' = 'darkorange',
                                  'Rubulavirus' = 'mediumseagreen', 'Avulavirus' = 'cornflowerblue'))


pout <- pheatmap(pa.mat, 
         color = color.ramp,
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_row = subclustplot,
         annotation_colors = mat.colors,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         #gaps_row = c(3, 12))
          gaps_row = c(3, 7, 11, 12))

print(pout)
# can modulate annotation colors with command 'annotation_colors = mat.colors'
# you would make a vector of the right colors
