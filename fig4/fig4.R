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


#load the data and make age-seroprevalence curves
homewd = "/Users/carabrook/Developer/bat-VirScan-public"
setwd(homewd)

exp.dat <- read.csv(file = paste0(homewd, "/working-data/all_bat_exposures.csv"), header = T, stringsAsFactors = F)
names(exp.dat)
length(unique(exp.dat$ID[exp.dat$bat_species=="Pteropus alecto"]))#77 all bats here


#reshape data with all pathogens and the proportion seropositive vs. negative above
#need to select some focal pathogens to look at, let's say at the genus level

#first get the proportion seropositive for each genus to get a sense of sample size for comparison
unique(exp.dat$virus_genus)
sum.genus <- ddply(exp.dat, .(virus_genus), summarise, N=length(ID))
sum.genus <- sum.genus[complete.cases(sum.genus),]
sum.genus$prev <- sum.genus$N/77

sum.genus <- arrange(sum.genus, desc(prev))
head(sum.genus) #this matches Fig. 2



dat=subset(exp.dat, bat_species=="Pteropus alecto")


length(unique(dat$ID)) #77
head(dat)

dat.long <- dplyr::select(dat, rank, ID, virus_genus)
length(unique(dat.long$ID)) #77
unique(dat.long$virus_genus)

dat.wide <- dcast(melt(dat.long), formula= ID~virus_genus)
head(dat.wide)
dat.wide <- dplyr::select(dat.wide, -("NA"))
names(dat.wide)
nrow(dat.wide)#77. one row per bat

#slim to just the pathogens of interest
#let's take the same ones as in Fig 1 for consistency
#the main four bat families + adeno + flavi but all at the genus level
#dat.wide <- dplyr::select(dat.wide, -("NA"))
#dat.wide <- dplyr::select(dat.wide, ID, Betacoronavirus, Henipavirus, Mastadenovirus, Lyssavirus, Marburgvirus, Ebolavirus)
head(dat.wide)
names(dat.wide)

unique(unlist(dat.wide[,2:ncol(dat.wide)]))
#which(dat.wide==2, arr.ind=TRUE) #all row 43
#dat.wide[43,] #Pa6
#subset(dat.long, ID=="Pa6") #oh! it has multiple entries for the same genera...
#only going to count that as seropositive once
#dat.wide[43,][dat.wide[43,]==2] <- 1
#unique(unlist(dat.wide[,2:ncol(dat.wide)])) #just 1s and 0s now



#now convert back to long
dat.long <- melt(dat.wide)
names(dat.long) <- c("ID", "virus_genus", "serostatus")
nrow(dat.long) #
unique(dat.long$serostatus)


#and merge with metadata
exp.dat.merge <- exp.dat #dplyr::select(exp.dat, -(virus_genus))
dat.all <- merge(dat.long, exp.dat.merge, by=c("ID", "virus_genus"), all.x = T)
head(dat.all)
unique(dat.all$serostatus)

unique(dat.all$age_cat)

head(dat.all)

#plot age-seroprevalence by boxplot
unique(dat.all$age_cat)
age.cat.prev <- ddply(subset(dat.all, !is.na(age_cat)), .(virus_genus, age_cat), summarise, N=sum(serostatus))
age.cat.base<- ddply(subset(dat.all, !is.na(age_cat)), .(age_cat), summarise, Ntot =length(unique(ID)))

age.cat.prev <- age.cat.prev[complete.cases(age.cat.prev),]
age.cat.base <- age.cat.base[complete.cases(age.cat.base),]



age.cat.prev <- merge(age.cat.prev, age.cat.base, by="age_cat", all.x = T)

age.cat.prev$seroprevalence <- age.cat.prev$N/age.cat.prev$Ntot

age.cat.prev$age_cat <- factor(age.cat.prev$age_cat, levels=c("Juvenile", "Sub Adult", "Adult"))

p1 <- ggplot(data = age.cat.prev) + 
      geom_bar(aes(x=age_cat, y=seroprevalence), stat = "identity") +
      facet_wrap(~virus_genus)

p2 <- ggplot(data = age.cat.prev) + 
  geom_point(aes(x=age_cat, y=seroprevalence, size=Ntot)) +
  geom_line(aes(x=age_cat, y=seroprevalence, group=virus_genus)) +
  facet_wrap(~virus_genus) + coord_cartesian(ylim=c(0,1))

#betaCoV and Flavi are highest in juveniles -- this is what is known for CoVs in the literature

unique(dat.all$age_tooth[dat.all$age_cat=="Juvenile"]) # NA,  1
unique(dat.all$age_tooth[dat.all$age_cat=="Sub Adult"]) # NA, 2, 4, 3, 5, 1
unique(dat.all$age_tooth[dat.all$age_cat=="Adult"]) # NA, 1-10

#and plot by age bin
#plot age-seroprevalence by boxplot
unique(dat.all$age_tooth)

dat.all$age_bin <- NA
dat.all$age_bin[dat.all$age_tooth=="1"] <- 1
dat.all$age_bin[dat.all$age_tooth==2 | dat.all$age_tooth==3] <- 2.5
dat.all$age_bin[dat.all$age_tooth==4 | dat.all$age_tooth==5] <- 4.5
dat.all$age_bin[dat.all$age_tooth==6 | dat.all$age_tooth==7] <- 6.5
dat.all$age_bin[dat.all$age_tooth==8 | dat.all$age_tooth==9] <- 8.5
dat.all$age_bin[dat.all$age_tooth==10 | dat.all$age_tooth==11] <- 10.5
dat.all$age_bin[dat.all$age_tooth==12 | dat.all$age_tooth==13] <- 12.5

dat.all <- subset(dat.all, !is.na(age_tooth)) #97 entries with age
length(unique(dat.all$ID)) #39 aged bats at all

#write function to get age-seroprevalence for each of the viruses
#split by virus genus
dat.split <- dlply(dat.all, .(virus_genus)) #39 genera

get.age.seroprev <- function(df, df.all){
  age.scale.base<- ddply(df.all, .(age_bin), summarise, N =length(unique(ID)))
  age.scale.prev <- ddply(df, .(virus_genus, age_bin), summarise, Npos=sum(serostatus))
  
  age.sero = age.scale.base
  
  age.sero <- merge(age.sero, age.scale.prev, by="age_bin", all.x = T)
  age.sero$virus_genus <- unique(df$virus_genus)
  age.sero$Npos[is.na(age.sero$Npos)] <- 0
  
  age.sero$seroprevalence <- age.sero$Npos/age.sero$N
  return(age.sero)
}

dat.split.sero <- lapply(dat.split, get.age.seroprev, df.all=dat.all)

dat.sero <- data.table::rbindlist(dat.split.sero)
head(dat.sero)


FigS7 <- ggplot(data = dat.sero) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16),
        axis.text = element_text(size=14), strip.text = element_text(size=12),
        legend.position = c(.8,.04), legend.direction = "horizontal")+
  geom_point(aes(x=age_bin, y=seroprevalence, size=N)) + xlab("age (yrs)") +
  ylab("seroprevalence\n") +
  geom_line(aes(x=age_bin, y=seroprevalence, group=virus_genus)) +
  facet_wrap(~virus_genus) + coord_cartesian(ylim=c(0,1), xlim=c(0,13)) +
  scale_x_continuous(breaks=c(0,4,8,12))



ggsave(file = paste0(homewd,"/supp-figures/figS7.png"),
       plot=FigS7,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


#and just the ones of interest
age.scale.plot = subset(dat.sero, 
                          virus_genus=="Henipavirus" | virus_genus=="Marburgvirus"|
                          virus_genus=="Betacoronavirus"| virus_genus=="Lyssavirus" | 
                          virus_genus=="Flavivirus" | virus_genus=="Mastadenovirus")


Fig4 <- ggplot(data = age.scale.plot) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16),
        axis.text = element_text(size=14), strip.text = element_text(size=12),
        legend.position = c(.95,.87), legend.direction = "vertical", 
        legend.background = element_rect(color="black"))+
  geom_point(aes(x=age_bin, y=seroprevalence, size=N)) + xlab("age (yrs)") +
  ylab("seroprevalence\n") +
  geom_line(aes(x=age_bin, y=seroprevalence, group=virus_genus)) +
  facet_wrap(~virus_genus) + coord_cartesian(ylim=c(0,1), xlim=c(0,13)) +
  scale_x_continuous(breaks=c(0,4,8,12))

#Fig4


ggsave(file = paste0(homewd,"/final-figures/fig4.png"),
       plot=Fig4,
       units="mm",  
       width=90, 
       height=60, 
       scale=3, 
       dpi=300)

#for model fitting... which genera have the most aged positives?
dat.sero.sum <- ddply(dat.sero, .(virus_genus), summarise, totPos = sum(Npos))
dat.sero.sum <- arrange(dat.sero.sum, desc(totPos)) #flavivirus and mastadenovirus
