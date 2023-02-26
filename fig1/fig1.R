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
library(ggpubr)


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
dat$virus_genus[dat$virus_genus=="Pneumovirus"] <- "Orthopneumovirus"
dat$virus_species[dat$virus_species=="Middle_East_respiratory_syndrome_related_coronavirus"] <- "MERS_related_coronavirus"


#load and merge subfamily data 
subfam <- read.csv(file = paste0(homewd, "/working-data/subfamily-merge.csv"), header = T, stringsAsFactors = F)
head(subfam)
subfam$virus_species[subfam$virus_species=="Middle_East_respiratory_syndrome_related_coronavirus"] <- "MERS_related_coronavirus"

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

# we are missing phip-seq data for these bats - is this because there were no hits?
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
dat$virus_genus[dat$virus_genus=="unclassified_Paramyxoviridae"] <- "unc. Paramyxo"
dat$virus_genus[dat$virus_genus=="unclassified_Anelloviridae"] <- "unc. Anello"
dat$virus_genus[dat$virus_genus=="unclassified_Arenaviridae"] <- "unc. Arena"
dat$virus_genus[dat$virus_genus=="unclassified_Astroviridae"] <- "unc. Astro"
dat$virus_genus[dat$virus_genus=="unclassified_Cyclovirus"] <- "unc. Cyclo"
dat$virus_genus[dat$virus_genus=="unclassified_Hepeviridae"] <- "unc. Hepe"
dat$virus_genus[dat$virus_genus=="unclassified_Papillomaviridae"] <- "unc. Papilloma"
dat$virus_genus[dat$virus_genus=="unclassified_Polyomaviridae"] <- "unc. Polyoma"
dat$virus_genus[dat$virus_genus=="unclassified_Retroviridae"] <- "unc. Retro"
dat$virus_species <- sub("_", " ", dat$virus_species)
dat$virus_species <- sub("_", " ", dat$virus_species)
dat$virus_species <- sub("_", " ", dat$virus_species)
dat$virus_species <- sub("_", " ", dat$virus_species)


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

  
  if(unique(df$virus_family !=leg.name)){
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           #strip.text = element_text(size=12),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
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
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
  }
  return(p1)
  
}
plot.heat.bigger <- function(df, all.bat, leg.name){
  
  
  
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
  
  
  if(unique(df$virus_family !=leg.name)){
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           strip.text = element_text(size=12),
                                           title = element_text(size=14),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
    
    
  }else{
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           strip.text = element_text(size=12),
                                           title = element_text(size=14),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family), ) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
  }
  return(p1)
  
}

pa.plot.list <- lapply(pa.list, plot.heat, all.bat=sort(unique(pa.dat$rank)),  leg.name="Betaherpesvirinae")


plot.list.a <- list()
plot.list.b <- list()
for(i in 1:16){
  plot.list.a[[i]] <- pa.list[[i]]
}
for(i in 1:14){
  plot.list.b[[i]] <- pa.list[[16+i]]
}


plot.list.a1 <- lapply(plot.list.a, plot.heat.bigger, all.bat=sort(unique(pa.dat$rank)),  leg.name="Nairoviridae")
plot.list.a2 <- lapply(plot.list.b, plot.heat.bigger, all.bat=sort(unique(pa.dat$rank)),  leg.name="Togaviridae")
#legend <-  plot.heat(pa.list[[1]],all.bat=sort(unique(pa.dat$rank)), return.leg=T)

figS1a <- cowplot::plot_grid(plotlist = plot.list.a1, ncol = 4, nrow=4, align = "v", axis="lr")    


ggsave(file = paste0(homewd,"/supp-figures/figS1a.pdf"),
       plot=figS1a,
       units="mm",  
       width=200, 
       height=200, #280 is a full page
       scale=3.5, 
       dpi=300)

figS1b <- cowplot::plot_grid(plotlist = plot.list.a2, ncol = 4, nrow=4, align = "v", axis="lr")    


ggsave(file = paste0(homewd,"/supp-figures/figS1b.pdf"),
       plot=figS1b,
       units="mm",  
       width=200, 
       height=200, #280 is a full page
       scale=3, 
       dpi=300)

#and again,


#then combine as 2 pages


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


#then, make subplot just for the main bat viruses, plus a few purely human libraries

#make a list by virus family

pa.list <- dlply(pa.merge, .(virus_family))

#testing a few things to include in plotting function below
dat.sum <- ddply(dat,.(virus_family), summarise, N_subfamily =length(unique(virus_subfamily)))
max(dat.sum$N_subfamily) #at most 5 subfamilies within a family
max(pa.merge$Assigned.counts) 



pa.list.bat <- list()
#list.fam = c("Paramyxoviridae", "Filoviridae", "Rhabdoviridae", "Coronaviridae", "Hantaviridae", "Adenoviridae")
#list.fam = c("Paramyxoviridae", "Filoviridae", "Rhabdoviridae", "Coronaviridae", "Flaviviridae", "Adenoviridae")
list.fam = c( "Paramyxoviridae", "Picornaviridae")
for (i in 1:length(list.fam)){
  pa.list.bat[[i]] <-  pa.list[[list.fam[i]]]
}




pa.list.bat.df <- data.table::rbindlist(pa.list.bat)
unique(pa.list.bat.df$virus_subfamily)
pa.list.bat.df$virus_subfamily[pa.list.bat.df$virus_subfamily=="Heptrevirinae"] <- "Hep."
pa.list.bat.df$virus_subfamily[pa.list.bat.df$virus_subfamily=="Picornaviridae"] <- "Pic."

unique(pa.list.bat.df$virus_genus)
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="unc. Paramyxo"] <- "unc.\nPar."
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Ledantevirus"] <- "Ledante\nvirus"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Pegivirus"] <- "Pegi\nvirus"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Avulavirus"]<- "Avulavirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Henipavirus"] <- "Henipavirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Respirovirus"]<- "Respirovirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Rubulavirus"] <- "Rubulavirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Cosavirus"]<- "Cosavirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Marburgvirus"] <- "Marburgvirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Ebolavirus"] <- "Ebolavirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Lyssavirus"] <- "Lyssavirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Vesiculovirus"] <- "Vesiculovirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Betacoronavirus"] <- "Betacoronavirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Flavivirus"] <- "Flavivirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Betacoronavirus"] <- "Betacoronavirus\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Hepacivirus"] <- "Hepacivirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Hepatovirus"] <- "Hep.\n"
#pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Mastadenovirus"] <- "Mastadenovirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Rosavirus"] <- "Ros.\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Enterovirus"] <- "Enterovirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Kobuvirus"] <- "Kob.\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Salivirus"] <- "Salivirus\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Cardiovirus"] <- "Cardio.\n"
pa.list.bat.df$virus_genus[pa.list.bat.df$virus_genus=="Parechovirus"] <- "Par.\n"




#pa.list.bat.df$virus_family <- factor(pa.list.bat.df$virus_family, levels = c("Paramyxoviridae", "Filoviridae",  "Coronaviridae",  "Rhabdoviridae", "Flaviviridae", "Adenoviridae"))
pa.list.bat.df$virus_family <- factor(pa.list.bat.df$virus_family, levels = c("Paramyxoviridae", "Picornaviridae"))
pa.list.bat <- dlply(pa.list.bat.df, .(virus_family))

#function with bigger titles
plot.heat.big <- function(df, all.bat, leg.name){
  
  
  
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
  
  
  if(unique(df$virus_family !=leg.name)){
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           plot.title = element_text(size=18),
                                           panel.grid.major.x =  element_blank(),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           plot.margin = unit(c(.5,.5,.5,.1), "cm"),
                                           strip.text = element_text(size=13),
                                           axis.text.y = element_text(size=16), 
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
    
    
  }else{
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           plot.title = element_text(size=18),
                                           strip.text = element_text(size=13),
                                           axis.text.y = element_text(size=16), 
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
  }
  return(p1)
  
}



#fig1.list  <-  lapply(pa.list.bat, plot.heat.big, all.bat=sort(unique(pa.dat$rank)),  leg.name="Coronaviridae")
fig1.list  <-  lapply(pa.list.bat, plot.heat.big, all.bat=sort(unique(pa.dat$rank)), leg.name="Picornaviridae")
#fig1 <- cowplot::plot_grid(plotlist = fig1.list, ncol = 3, align = "v", axis="lr", labels = c("A", "B", "C", "D", "E", "F"), label_size = 20, label_y = 1.005, rel_heights = c(1,1,1,1.2,1.2,1.2))    
fig1right <- cowplot::plot_grid(plotlist = fig1.list, nrow = 2, align = "v", axis="lr", labels = c("B", "C"), label_size = 20, label_y = 1.005, rel_heights = c(.667,1))    


#and add in panel A


#and make the mock dataset for the cartoon - do it on a subset of 3 bats
pa.mock <- cbind.data.frame(virus_family = rep("Virus-Family-Ex", 10),
                            virus_subfamily = c(rep("Virus-Subfamily-A", 6), rep("Virus-Subfamily-B", 4)),
                            virus_genus=c(rep("Virus-Genus-a", 3), rep("Virus-Genus-b", 3), rep("Virus-Genus-c", 2), rep("Virus-Genus-d", 2)),
                            virus_species= paste0("Virus-Species-", seq(1,10,1)))

#and add in the rank of the bat
#bat 1 is exposed across the entire family
pa.mock$rank <- NA
pa.mock$Assigned.counts <- NA
pa.mock$rank[1:4] <- pa.mock$rank[7:10] <- 1 
pa.mock$Assigned.counts[1:4] <-  c(20,25,15,10)
pa.mock$Assigned.counts[7:10] <- c(5,22,27,21)

#bat 2 is just subfamily B - but multiple genera
pa.mock <- rbind(pa.mock, pa.mock[8,])
pa.mock$rank[9:11] <- 2
pa.mock$Assigned.counts[9:11] <- c(30, 50, 70)
#bat 3 is just genera a but multiple species
pa.mock$rank[5:6] <- 3
pa.mock$Assigned.counts[5:6]<- c(205, 150)

pa.mock <- merge(pa.mock, pa.merge, by=c("virus_species", "virus_family", "virus_subfamily", "virus_genus", "rank", "Assigned.counts"), all.x = T)
head(pa.mock)
pa.mock <- dplyr::select(pa.mock, names(pa.merge))
head(pa.mock)

pa.mock$virus_genus[pa.mock$virus_genus=="Virus-Genus-a"] <- "Virus-Genus-a\n"
pa.mock$virus_genus[pa.mock$virus_genus=="Virus-Genus-b"] <- "Virus-Genus-b\n"
pa.mock$virus_genus[pa.mock$virus_genus=="Virus-Genus-c"] <- "Virus-Genus-c\n"
pa.mock$virus_genus[pa.mock$virus_genus=="Virus-Genus-d"] <- "Virus-Genus-d\n"


pa.mock <- arrange(pa.mock, rank)

plot.mock.heat <- function(df, all.bat, leg.name){
  
  
  
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
  
  df1$virus_species <- factor(df1$virus_species, levels = c("Virus-Species-1", "Virus-Species-2", "Virus-Species-3", "Virus-Species-4", "Virus-Species-5", "Virus-Species-6", "Virus-Species-7", "Virus-Species-8", "Virus-Species-9", "Virus-Species-10"))
  df1$virus_subfamily <- factor(df1$virus_subfamily, levels = c("Virus-Subfamily-B", "Virus-Subfamily-A"))
  df1$virus_genus <- factor(df1$virus_genus, levels = c("Virus-Genus-d\n", "Virus-Genus-c\n", "Virus-Genus-b\n", "Virus-Genus-a\n"))
  
  df1$label <- NA
  df1$label[df1$rank==1] <- "Bat seropositive at\n'Example-Virus-Family'\nlevel only"
  df1$label[df1$rank==2] <- "Bat seropositive at\n'Example-Virus-Family'\nand 'Virus-Subfamily-B'"
  df1$label[df1$rank==3] <- "Bat seropositive at\n'Example-Virus-Family',\n'Virus-Subfamily-A',\nand 'Virus-Genus-b'"
  
  df1$rank <- df1$label
  df1$rank <- factor(df1$rank, levels=c("Bat seropositive at\n'Example-Virus-Family'\nlevel only",
                                        "Bat seropositive at\n'Example-Virus-Family'\nand 'Virus-Subfamily-B'",
                                        "Bat seropositive at\n'Example-Virus-Family',\n'Virus-Subfamily-A',\nand 'Virus-Genus-b'"))
  if(unique(df$virus_family !=leg.name)){
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           plot.title = element_text(size=18),
                                           panel.grid.major.x =  element_line(color="black", size=1),
                                           panel.grid.minor.x =  element_line(color="black", size=1),
                                           panel.grid.major.y = element_line(color="black", size=1),
                                           panel.grid.minor.y = element_line(color="black", size=1),
                                           plot.margin = unit(c(.5,.1,.5,.5), "cm"),
                                           strip.text = element_text(size=13),
                                           axis.text.x = element_text(size=11), 
                                           axis.text.y = element_text(size=12), 
                                           axis.ticks.x = element_blank(), 
                                           axis.title.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1, show.legend = F) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = "Example-Virus-Family") +
      geom_vline(xintercept = seq(1.5, (3+.5),1), color = "black")
    
    
  }else{
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           panel.grid.major.x =  element_blank(),
                                           plot.title = element_text(size=18),
                                           panel.grid.major.y = element_line(color="black"),
                                           panel.grid.minor.y = element_line(color="black"),
                                           axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                                           axis.text.x = element_blank(), axis.title.y = element_blank()) +
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = unique(df$virus_family)) +
      geom_vline(xintercept = seq(1.5, (max(df1$rank)+.5),1), color = "black")
  }
  return(p1)
  
}

fig1left <- plot.mock.heat(df=pa.mock, all.bat = 1:3, leg.name = "")


fig1 <- cowplot::plot_grid(fig1left, fig1right, ncol=2, rel_widths = c(1,1.1), labels = c("A", ""), label_size = 20)



ggsave(file = paste0(homewd,"/final-figures/fig1.png"),
       plot=fig1,
       units="mm",  
       width=140, 
       height=110, 
       scale=3, 
       dpi=300)

 
plot.mock.heat.leg <- function(df, all.bat){
  
  
  
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
  
  df1$virus_species <- factor(df1$virus_species, levels = c("Virus-Species-1", "Virus-Species-2", "Virus-Species-3", "Virus-Species-4", "Virus-Species-5", "Virus-Species-6", "Virus-Species-7", "Virus-Species-8", "Virus-Species-9", "Virus-Species-10"))
  df1$virus_subfamily <- factor(df1$virus_subfamily, levels = c("Virus-Subfamily-B", "Virus-Subfamily-A"))
  df1$virus_genus <- factor(df1$virus_genus, levels = c("Virus-Genus-d\n", "Virus-Genus-c\n", "Virus-Genus-b\n", "Virus-Genus-a\n"))
  
  df1$label <- NA
  df1$label[df1$rank==1] <- "Bat seropositive\nat 'Example-\nVirus-Family'\nlevel only"
  df1$label[df1$rank==2] <- "Bat seropositive\nat 'Example-\nVirus-Family'\nand 'Virus\n-Subfamily-B'"
  df1$label[df1$rank==3] <- "Bat seropositive\nat 'Example-\nVirus-Family',\n'Virus-Subfamily\n-A', and 'Virus-\nGenus-b'"
  
  df1$rank <- df1$label
  df1$rank <- factor(df1$rank, levels=c("Bat seropositive\nat 'Example-\nVirus-Family'\nlevel only",
                                        "Bat seropositive\nat 'Example-\nVirus-Family'\nand 'Virus\n-Subfamily-B'",
                                        "Bat seropositive\nat 'Example-\nVirus-Family',\n'Virus-Subfamily\n-A', and 'Virus-\nGenus-b'"))
  
    p1 <- ggplot(df1) + theme_bw() + theme(panel.grid.minor = element_blank(), 
                                           plot.title = element_text(size=18),
                                           panel.grid.major.x =  element_line(color="black", size=1),
                                           panel.grid.minor.x =  element_line(color="black", size=1),
                                           panel.grid.major.y = element_line(color="black", size=1),
                                           panel.grid.minor.y = element_line(color="black", size=1),
                                           plot.margin = unit(c(.5,.1,.5,.5), "cm"),
                                           strip.text = element_text(size=13),
                                           legend.position = "bottom",
                                           axis.text.x = element_text(size=11), 
                                           axis.text.y = element_text(size=16), 
                                           axis.ticks.x = element_blank(), 
                                           axis.title.x = element_blank(), axis.title.y = element_blank()) +#,
      #legend.position = "bottom", #legend.title = element_blank(),
      #legend.direction = "horizontal")+
      geom_tile(aes(x=rank, y=virus_species, fill=Assigned.counts), width=1) +
      scale_fill_viridis_c(values=c(0,.05,.1,1), limits=c(0,330), na.value = "gray", name="peptide\ncounts") + scale_y_discrete(position = "right") +
      facet_nested(virus_subfamily + virus_genus~., scales = "free", space = "free_y",
                   switch = "y", strip = custom_strips) + coord_cartesian(expand=F) + ggtitle(label = "Example-Virus-Family") +
      geom_vline(xintercept = seq(1.5, (3+.5),1), color = "black")
  
  return(p1)
  
}


fig1a <- plot.mock.heat.leg(df=pa.mock, all.bat = 1:3)

pa.para = subset(pa.list.bat.df, virus_family=="Paramyxoviridae")
pa.picorna = subset(pa.list.bat.df, virus_family=="Picornaviridae")
#and try a version with B and C in separate columns
fig1b = plot.heat.big(df=pa.para, all.bat=sort(unique(pa.dat$rank)), leg.name="")
fig1c = plot.heat.big(df=pa.picorna, all.bat=sort(unique(pa.dat$rank)), leg.name="")

fig1 <- cowplot::plot_grid(fig1a, fig1b, fig1c, ncol=3, nrow = 1,  
                           labels = c("A", "B", "C"), label_size = 20, rel_widths = c(1,1.1,1.1))

ggsave(file = paste0(homewd,"/final-figures/fig1.png"),
       plot=fig1,
       units="mm",  
       width=180, 
       height=90, 
       scale=3, 
       dpi=300)
