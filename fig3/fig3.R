rm(list=ls())

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


#load the data
homewd = "/Users/carabrook/Developer/bat-VirScan-public"
setwd(homewd)

#load the metadata here to check that the merge was correct
meta <- read.csv(file = paste0(homewd, "/working-data/bat_metadata.csv"), header = T, stringsAsFactors = F)
head(meta)
names(meta)
meta.merge <- dplyr::select(meta, ID, Sex, Age_cat, Age_scale, Mass_g, Forearm_mm, residuals, Condition)

exp.dat <- read.csv(file = paste0(homewd, "/working-data/all_bat_exposures.csv"), header = T, stringsAsFactors = F)
names(exp.dat)

dat <- merge(exp.dat, meta.merge, by="ID", all.x = T)
head(dat)

dat$seropos=1

dat.sub.y = subset(dat, virus_genus=="Henipavirus" & bat_species=="Pteropus alecto")
dat.sub.n = subset(dat, virus_genus!="Henipavirus" & bat_species=="Pteropus alecto")
dat.sub.n$seropos<- 0
dat.sub.n$virus_family <- "Paramyxoviridae"
dat.sub.n$virus_subfamily <- "Orthoparamyxovirinae"
dat.sub.n$virus_genus <- "Henipavirus"
dat.sub.n <- dat.sub.n[!duplicated(dat.sub.n),]
head(dat.sub.n)
pos.ID <- unique(dat.sub.y$ID)
dat.sub.n$keep = 1
for(i in 1:length(pos.ID)){
  dat.sub.n$keep[dat.sub.n$ID==pos.ID[i]]<- 0  
}
dat.sub.n = subset(dat.sub.n, keep==1)
length(setdiff(dat.sub.n$ID, dat.sub.y$ID)) #66
length(setdiff(dat.sub.y$ID, dat.sub.n$ID)) #0 # not sure where the other 5 bats are

dat.sub.y$keep = 1
dat.sub <- rbind(dat.sub.n, dat.sub.y)
dat.sub$seropos <- factor(dat.sub$seropos)
ggplot(data=dat.sub) + geom_boxplot(aes(x=seropos, y=Mass_g))
