rm(list = ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(ggthemes)
library(tidyverse)
library(ggplotify)
library(reshape2)
library(pheatmap)


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

#also attach genus to subfam qualifier (for later)
merge.dat <- dplyr::select(dat, virus_species, virus_genus)
merge.dat <- dplyr::select(dat, virus_species, virus_genus, virus_family)
merge.dat <- merge.dat[!duplicated(merge.dat),]
subfam <- merge(subfam, merge.dat, by="virus_species", all.x = T)
head(subfam)

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

# we are missing phip-seq data for these bats
#"Pa29"  "Pa50"  "Pa59"  "Pa77"  "Pa81"  "Pa10"  "Pa100" "Pa101" "Pa1"   "Pa2"   "Pa3"  


dat <- arrange(dat, virus_subfamily, virus_genus, virus_species)
dat$virus_species <- factor(dat$virus_species, levels=unique(dat$virus_species))
head(dat)
 
unique(dat$X) # these only have values on the eonycteris samples, where the columns appear to be shifted
# fix down below


# and summarize the data to report in the paper
length(unique(dat$virus_family[dat$Species=="Pteropus alecto"])) #33
length(unique(dat$virus_subfamily[dat$Species=="Pteropus alecto"])) #47
length(unique(dat$virus_genus[dat$Species=="Pteropus alecto"])) #82
length(unique(dat$virus_species[dat$Species=="Pteropus alecto"])) #414


length(unique(dat$virus_family[dat$Species=="Eonycteris spelaea"])) #12
length(unique(dat$virus_subfamily[dat$Species=="Eonycteris spelaea"])) #15
length(unique(dat$virus_genus[dat$Species=="Eonycteris spelaea"])) #23
length(unique(dat$virus_species[dat$Species=="Eonycteris spelaea"])) #78



#first, for figure 1, plot a heatmap and compare broad viral exposures for the two species
es.dat = subset(dat, Species=="Eonycteris spelaea") # why are the assigned counts so few here?
head(es.dat) # it appears that the column names are shifted to the right on these guys about halfway through
names(es.dat)[17:22] <- c("X", "Assigned.counts", "Assigned.peptides", "Total.sample.hits", "Total.filtered.sample.hits", "virus_subfamily")
es.dat <- dplyr::select(es.dat, -(X))
head(es.dat)
pa.dat = subset(dat, Species=="Pteropus alecto")
head(pa.dat)
pa.dat <- dplyr::select(pa.dat, -(X))


es.dat$Assigned.counts <- as.numeric(es.dat$Assigned.counts)

#mega heatmap
subclustplot <- dplyr::select(subfam,virus_genus, virus_subfamily,  virus_family)
rownames(subclustplot) <- subfam$virus_species
names(subclustplot) <- c("genus", "subfamily", "family")


pa.dat <- arrange(pa.dat, virus_family, virus_subfamily, virus_genus, virus_species)
pa.dat$virus_species <- factor(pa.dat$virus_species, levels = c(unique(pa.dat$virus_species)))

subfam <- arrange(subfam, virus_family, virus_subfamily, virus_genus, virus_species)
subfam$virus_species <- factor(subfam$virus_species, levels = c(unique(pa.dat$virus_species)))

#remove those viral families with only one entry
vsum <- ddply(pa.dat, .(virus_family), summarise, N=length(virus_species))
vsum <- vsum[vsum$N>1,]
pa.merge <-merge(pa.dat, vsum, by = "virus_family")
pa.merge = subset(pa.merge, virus_family!="unclassified_family")
pa.slim <- dplyr::select(pa.merge, virus_species)
pa.slim <- cbind.data.frame(virus_species=pa.slim [!duplicated(pa.slim ),])

subfam.plot <- merge(subfam, pa.slim, by ="virus_species")

col.dat <- cbind.data.frame(ID = unique(pa.merge$ID), xlab =1:length(unique(pa.merge$ID)))
pa.merge <- merge(pa.merge, col.dat, by="ID")
pa.merge <- arrange(pa.merge, xlab)
pa.merge$xlab <- pa.merge$xlab + 2.5

figS1 <- ggplot(pa.merge) + theme_bw() + theme(panel.grid = element_blank(),
         strip.background = element_rect(fill="white"),
         axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
         axis.text.x = element_blank(), axis.title.y = element_blank(),
         legend.position = "bottom", legend.title = element_blank(),
         legend.direction = "horizontal")+
         geom_tile(aes(x=xlab, y=virus_species, fill=Assigned.counts), width=1) +
         #add in the bars
         scale_fill_viridis_c() + scale_y_discrete(position = "right") +
         #facet_grid(virus_family~., scales="free_y") +
         facet_wrap(virus_family~., scales="free_y") + coord_cartesian(xlim=c(.5, (max(pa.merge$xlab)+.5)), expand = F)+
         ggnewscale::new_scale_fill() + 
         geom_tile(data= subfam.plot, aes(x=2, y=virus_species, fill=virus_genus), width=1, color="black") +
         geom_tile(data= subfam.plot, aes(x=1, y=virus_species, fill=virus_subfamily), width=1, color="black") +
         guides(fill=guide_legend(nrow=4,byrow=TRUE))
         #geom_tile(data= subfam.plot, aes(x=1, y=virus_species, fill=virus_family),color="black",show.legend = F) 
         #facet_wrap(~virus_family, scales = "free_y")
    

ggsave(file = paste0(homewd,"/supp-figures/figS1.png"),
       plot=figS1,
       units="mm",  
       width=400, 
       height=200, 
       scale=3, 
       dpi=300)


#repeat for Eonycteris
#remove those viral families with only one entry
vsum <- ddply(es.dat, .(virus_family), summarise, N=length(virus_species))
vsum <- vsum[vsum$N>1,]
es.merge <-merge(es.dat, vsum, by = "virus_family")
es.merge = subset(es.merge, virus_family!="unclassified_family")
es.slim <- dplyr::select(es.merge, virus_species)
es.slim <- cbind.data.frame(virus_species=es.slim [!duplicated(es.slim ),])

subfam.plot <- merge(subfam, es.slim, by ="virus_species")



col.dat <- cbind.data.frame(ID = unique(es.merge$ID), xlab =1:length(unique(es.merge$ID)))
es.merge <- merge(es.merge, col.dat, by="ID")
es.merge <- arrange(es.merge, xlab)
es.merge$xlab <- es.merge$xlab + 1.5


figS2 <- ggplot(es.merge) + theme_bw() + theme(panel.grid  = element_blank(),
                                               strip.background = element_rect(fill="white"),
                                                axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                                axis.text.x = element_blank(), axis.title.y = element_blank(),
                                                legend.position = "bottom", #legend.title = element_blank(),
                                                legend.direction = "horizontal")+
  geom_tile(aes(x=xlab, y=virus_species, fill=Assigned.counts), width=1) +
  #add in the bars
  scale_fill_viridis_c(name="peptide counts") + scale_y_discrete(position = "right") +
  #facet_grid(virus_family~., scales="free_y") +
  facet_wrap(virus_family~., scales="free_y") + 
  ggnewscale::new_scale_fill() + coord_cartesian(xlim=c(.75, (max(es.merge$xlab)+.5)), expand = F)+
  geom_tile(data= subfam.plot, aes(x=1.5, y=virus_species, fill=virus_genus), width=.5, color="black") +#,width=.5,
  geom_tile(data= subfam.plot, aes(x=1, y=virus_species, fill=virus_subfamily),width=.5,  color="black") +
  scale_fill_discrete(name="virus genus\n& subfamily") +
  guides(fill=guide_legend(nrow=4,byrow=TRUE,)) 
#geom_tile(data= subfam.plot, aes(x=1, y=virus_species, fill=virus_family),color="black",show.legend = F) 
#facet_wrap(~virus_family, scales = "free_y")

ggsave(file = paste0(homewd,"/supp-figures/figS2.png"),
       plot=figS2,
       units="mm",  
       width=150, 
       height=80, 
       scale=3, 
       dpi=300)

#then, use the pheatmap package for the main text figure (just paramyxovirus)
#can comment on picorna and rhabdo
pa.para = subset(pa.merge, virus_family=="Paramyxoviridae")
es.para = subset(es.merge, virus_family=="Paramyxoviridae")
subfam.para = subset(subfam.plot, virus_family=="Paramyxoviridae")

pbar1 <- ggplot(subfam.para) + theme_bw() + theme(panel.background = element_blank(),
         axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())+
  geom_tile(aes(x=1, y=virus_species, fill=virus_genus), show.legend = F) +
  geom_tile(aes(x=2, y=virus_species, fill=virus_subfamily), show.legend = F) + coord_cartesian(expand = F)
  
pbar1.leg <- ggplot(subfam.para) + theme_bw() + theme(panel.background = element_blank(), legend.position = "horizontal",
                                                  axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())+
  geom_tile(aes(x=1, y=virus_species, fill=virus_genus)) +
  geom_tile(aes(x=2, y=virus_species, fill=virus_subfamily)) + coord_cartesian(expand = F)

pbar1.leg <- cowplot::get_legend(pbar1.leg)

pheat <- ggplot(pa.para) + theme_bw() + theme(panel.grid = element_blank(),
                                               strip.background = element_rect(fill="white"),
                                               axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                               axis.text.x = element_blank(), axis.title.y = element_blank())+
  geom_tile(aes(x=xlab, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
  #add in the bars
  scale_fill_viridis_c() + scale_y_discrete(position = "right")

pheat.leg <- ggplot(pa.para) + theme_bw() + theme(panel.grid = element_blank(),
                                              strip.background = element_rect(fill="white"),
                                              axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                              axis.text.x = element_blank(), axis.title.y = element_blank(),
                                              legend.position = "bottom", legend.title = element_blank(),
                                              legend.direction = "horizontal")+
  geom_tile(aes(x=xlab, y=virus_species, fill=Assigned.counts), width=1) +
  #add in the bars
  scale_fill_viridis_c(name="peptide counts") + scale_y_discrete(position = "right")
  
pheat.leg <- cowplot::get_legend(pheat.leg)

pboth <- cowplot::plot_grid(pbar1, pheat, ncol=2, nrow=1, rel_widths = c(.1,1))

pleg <- cowplot::plot_grid(pbar1.leg, pheat.leg, ncol=2, nrow = 1, rel_widths = c(1, .2))

Fig1 <- cowplot::plot_grid(pboth, pleg, ncol=1, nrow = 2, rel_heights = c(1, .1))


ggsave(file = paste0(homewd,"/final-figures/fig1draft.png"),
       plot=Fig1,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)

#and mega heatmap
#convert to matrix
df <- dplyr::select(pa.para, ID, virus_species, Assigned.counts)
df <- dcast(melt(df), formula = virus_species~ID)


df.mat <- as.matrix(df[2:ncol(df)])
rownames(df.mat) <- df[,1]

#make all NAs = 0
df.mat[is.na(df.mat)]<-0

color.ramp <- viridis::cividis(100)

pheatmap(df.mat, 
         color=color.ramp,
         cluster_cols = F, 
         cluster_rows = F, 
         gaps_row = c(3,7, 11,12),
         annotation_row = subclustplot)

df2 <- dplyr::select(es.para, ID, virus_species, Assigned.counts)
df2 <- dcast(melt(df2), formula = virus_species~ID)


df2.mat <- as.matrix(df2[,2:ncol(df2)])
rownames(df2.mat) <- df2[,1]

#make all NAs = 0
df2.mat[is.na(df2.mat)]<-0

df2.mat <- log10(df2.mat)
df2.mat[df2.mat=="-Inf"] <- 0

pheatmap(df2.mat, 
         color=color.ramp,
         gaps_row = c(1,5,6),
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_row = subclustplot)




subclustplot <- dplyr::select(subfam, virus_genus, virus_subfamily)
rownames(subclustplot) <- subfam$virus_species
names(subclustplot) <- c("genus", "subfamily")

#then, split the data into matrices by family and populate a grid
# first, remove the one virus family with only one entry

pa.dat.list <- dlply(pa.dat[pa.dat$virus_family!="Picobirnaviridae" & 
                              pa.dat$virus_family!="Podoviridae",], .(virus_family))


plot.function <- function(df, clusters){
  
  print(unique(df$virus_family))
  
  
  #convert to matrix
  df1 <- dplyr::select(df, ID, virus_species, Assigned.counts)
  
  df1 <- dcast(melt(df1), formula = virus_species~ID)
  
  
  df.mat <- as.matrix(df1[2:ncol(df1)])
  rownames(df.mat) <- df1[,1]
  
  #make all NAs = 0
  df.mat[is.na(df.mat)]<-0
  
  color.ramp <- viridis::cividis(10)
  
  p.out <- as.ggplot(pheatmap(df.mat, 
           color = color.ramp,
           cluster_cols = F, 
           cluster_rows = F, 
           name = unique(df$virus_family),
           annotation_row = clusters))
  
  return(p.out)
}

pa.plot.list <- lapply(pa.dat.list, plot.function, clusters=subclustplot)



pS1b <- cowplot::plot_grid(plotlist = pa.plot.list, ncol = 6)


ggsave(file = paste0(homewd,"/supp-figures/figS1b.png"),
       #plot=pS1b,
       units="mm",  
       width=300, 
       height=200, 
       scale=3, 
       dpi=300)




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
  

  
gtable::gtable(pheatmap(pa.mat, 
         color = color.ramp,
         cluster_cols = F, 
         cluster_rows = F, 
         annotation_row = subclustplot))

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
