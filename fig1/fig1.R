rm(list = ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(ggthemes)
library(tidyverse)
library(ggplotify)
library(reshape2)
library(pheatmap)
library(ggh4x)
library(metafolio)
library(RColorBrewer)


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

sort(unique(dat$virus_genus))
dat$virus_genus[dat$virus_genus=="unclassified_Paramyxoviridae"] <- "unclassified Paramyxoviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Anelloviridae"] <- "unclassified Anelloviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Arenaviridae"] <- "unclassified Arenaviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Astroviridae"] <- "unclassified Astroviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Cyclovirus"] <- "unclassified Cyclovirus"
dat$virus_genus[dat$virus_genus=="unclassified_Hepeviridae"] <- "unclassified Hepeviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Papillomaviridae"] <- "unclassified Papillomaviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Polyomaviridae"] <- "unclassified Polyomaviridae"
dat$virus_genus[dat$virus_genus=="unclassified_Retroviridae"] <- "unclassified Retroviridae"



length(unique(dat$virus_subfamily)) #47 subfamilies with nested colors within them

# first, for figure 1, plot a heatmap and compare broad viral exposures for the two species
# include a supplementary figure that is focused on just one or two viral families (maybe just paramyxoviruses)
es.dat = subset(dat, Species=="Eonycteris spelaea") # why are the assigned counts so few here?
head(es.dat) # it appears that the column names are shifted to the right on these guys about halfway through
names(es.dat)[17:22] <- c("X", "Assigned.counts", "Assigned.peptides", "Total.sample.hits", "Total.filtered.sample.hits", "virus_subfamily")
es.dat <- dplyr::select(es.dat, -(X))
head(es.dat)
pa.dat = subset(dat, Species=="Pteropus alecto")
head(pa.dat)
pa.dat <- dplyr::select(pa.dat, -(X))

es.dat$Assigned.counts <- as.numeric(es.dat$Assigned.counts)



pa.dat <- arrange(pa.dat, virus_family, virus_subfamily, virus_genus, virus_species)
pa.dat$virus_species <- factor(pa.dat$virus_species, levels = c(unique(pa.dat$virus_species)))

es.dat <- arrange(es.dat, virus_family, virus_subfamily, virus_genus, virus_species)
es.dat$virus_species <- factor(es.dat$virus_species, levels = c(unique(es.dat$virus_species)))


#and get rank numbers for plotting for each bat
rank.bat = cbind.data.frame(ID = unique(pa.dat$ID), rank = seq(1, length(unique(pa.dat$ID))))
pa.dat <- merge(pa.dat, rank.bat, by="ID")


rank.bat2 = cbind.data.frame(ID = unique(es.dat$ID), rank = seq(1, length(unique(es.dat$ID))))
es.dat <- merge(es.dat, rank.bat2, by="ID")


#remove those viral families with only one entry
vsum <- ddply(pa.dat, .(virus_family), summarise, N=length(virus_species))
vsum <- vsum[vsum$N>1,]
pa.merge <-merge(pa.dat, vsum, by = "virus_family")
pa.merge = subset(pa.merge, virus_family!="unclassified_family")


#make a list by virus family
pa.list <- dlply(pa.merge, .(virus_family))

#testing a few things to include in plotting function below
dat.sum <- ddply(dat,.(virus_family), summarise, N_subfamily =length(unique(virus_subfamily)))
max(dat.sum$N_subfamily) #at most 5 subfamilies within a family
max(pa.merge$Assigned.counts) #326

#plotting function for one virus family -- apply over lists and populate grid
plot.heat <- function(df, all.bat, leg.name){
  
  
  
  #first, attach blanks for the bats with no hits to this virus
  #all.bat = setdiff(all.bat, unique(df$rank))
  rank.df <- cbind.data.frame(rank=rep(all.bat, length(unique(df$virus_species))), virus_species=rep(unique(df$virus_species), each=length(all.bat)))
  subfam.sum <- ddply(df, .(virus_subfamily, virus_genus, virus_species), summarise)
  rank.df <- merge(rank.df, subfam.sum, by ="virus_species", all.x = T)
  #rank.df$Assigned.counts <- 0
  head(rank.df)
  
  df1 <- merge(df, rank.df, by = c("rank", "virus_subfamily", "virus_genus", "virus_species"), all = T)
  df1 <- arrange(df1, rank)
  df1$ID <- factor(df1$ID, levels=unique(df1$ID))
  
  # #need list of 5 possible palette
   subfam.col.list = c("Blues", "Reds", "Purples", "Oranges", "Greens")
  # 
  # #nested colors for strip labels
   ncol = length(unique(df$virus_subfamily)) #number distinct subfamily
   subfam.col.list <- subfam.col.list[1:ncol] 
  # 
  # #and get the sub-palettes
   df.sum <- ddply(df1, .(virus_subfamily), summarise, N_genus = length(unique(virus_genus)))
  # 
   df.sum$N_genus <- df.sum$N_genus+1 #add one for the subfamily
  # 
   palette.list <- list()
   for(i in 1:length(df.sum$virus_subfamily)){
     palette.list[[i]] <- rev(brewer.pal(n=df.sum$N_genus[i], name=subfam.col.list[i]))[1:df.sum$N_genus[i]]
   }
  # 
  #this gives a nested list of colors
  
  palette.vector <- c(unlist(palette.list))
  names(palette.vector) <- c(unlist(mapply(rep, as.list(df.sum$virus_subfamily), each=as.list(df.sum$N_genus), SIMPLIFY = F)))
  
  for(i in 1:length(df.sum$virus_subfamily)){
    names(palette.vector)[names(palette.vector)==df.sum$virus_subfamily[i]] <- c(df.sum$virus_subfamily[i], unique(df$virus_genus[df$virus_subfamily==df.sum$virus_subfamily[i]]))
  }
  
  #then, reorder for plotting to take the virus_subfamily entries first
  new.palette <- list()
  for(i in 1:length(df.sum$virus_subfamily)){
    index <- which(names(palette.vector)==df.sum$virus_subfamily[i]) 
    new.palette[[i]] <- palette.vector[index]
  }
  new.palette <- c(unlist(new.palette))
  new.palette  <- c(new.palette, setdiff(palette.vector, new.palette))
  
  
  #now set custom strips function to take in your colors
  custom_strips <- strip_nested(
    background_y = elem_list_rect(
      fill = new.palette),
      by_layer_y = FALSE
  )
  
  p1leg <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                         panel.grid.major.x =  element_blank(),
                                         panel.grid.major.y = element_line(color="black"),
                                         panel.grid.minor.y = element_line(color="black"),
                                          axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                         axis.text.x = element_blank(), axis.title.y = element_blank(),
                                         legend.position = "bottom", #legend.title = element_blank(),
                                          legend.direction = "horizontal")+
    geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
    scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
    facet_nested(virus_subfamily + virus_genus~., scales = "free_y", 
    switch = "y", strip = custom_strips) + coord_cartesian(expand=F) +
    geom_vline(xintercept = seq(1.5, 77.5,1), color = "white")
    
  
  p1leg <- cowplot::get_legend( p1leg )
  
  if(unique(df$virus_family !=leg.name)){
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free_y", 
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
    
  }else{
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free_y", 
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
  }
  return(p1)
  
}

pa.plot.list <- lapply(pa.list, plot.heat, all.bat=sort(unique(pa.dat$rank)),  leg.name="Betaherpesvirinae")
#legend <-  plot.heat(pa.list[[1]],all.bat=sort(unique(pa.dat$rank)), return.leg=T)

figS1 <- cowplot::plot_grid(plotlist = pa.plot.list, ncol = 6, align = "v", axis="lr")    

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

#make a list by virus family
es.list <- dlply(es.merge, .(virus_family))

es.plot.list <- lapply(es.list, plot.heat, all.bat=sort(unique(es.dat$rank)),  leg.name="Togaviridae")
figS2 <- cowplot::plot_grid(plotlist = es.plot.list, ncol = 4, align = "v", axis="lr")    

ggsave(file = paste0(homewd,"/supp-figures/figS2.png"),
       plot=figS2,
       units="mm",  
       width=200, 
       height=100, 
       scale=3, 
       dpi=300)


#then, make subplot just for paramyxoviruses