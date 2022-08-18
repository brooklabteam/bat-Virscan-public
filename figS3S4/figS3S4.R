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

exp.dat <- read.csv(file = paste0(homewd, "/working-data/all_bat_exposures.csv"), header = T, stringsAsFactors = F)
names(exp.dat)
length(unique(exp.dat$ID[exp.dat$bat_species=="Pteropus alecto"]))#77 all bats here


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

meta.merge <- dplyr::select(meta, ID, Sex, Age_cat, Age_scale, Mass_g, Forearm_mm, Condition)
names(meta.merge) <- c("ID", "sex", "age_cat", "age_scale", "mass_g", "forearm_mm", "condition")

#reshape data with all pathogens and the proportion seropositive vs. negative above
#need to select some focal pathogens to look at, let's say at the genus level

#first get the proportion seropositive for each genus to get a sense of sample size for comparison
unique(exp.dat$virus_genus)
sum.genus <- ddply(exp.dat, .(virus_genus), summarise, N=length(ID))
sum.genus <- sum.genus[complete.cases(sum.genus),]
sum.genus$prev <- sum.genus$N/77

sum.genus <- arrange(sum.genus, desc(prev))
head(sum.genus) #this matches Fig. 2



dat <- merge(exp.dat, meta.merge, by="ID", all.x = T)

dat=subset(dat, bat_species=="Pteropus alecto")

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
which(dat.wide==2, arr.ind=TRUE) #all row 43
dat.wide[43,] #Pa6
subset(dat.long, ID=="Pa6") #oh! it has multiple entries for the same genera...
#only going to count that as seropositive once
dat.wide[43,][dat.wide[43,]==2] <- 1
unique(unlist(dat.wide[,2:ncol(dat.wide)])) #just 1s and 0s now



#now convert back to long
dat.long <- melt(dat.wide)
names(dat.long) <- c("ID", "virus_genus", "serostatus")
nrow(dat.long) #462 = 77*6
unique(dat.long$serostatus)

#and merge with metadata
dat.all <- merge(dat.long, meta.merge, by="ID", all.x = T)
head(dat.all)
unique(dat.all$serostatus)

#now get residuals
unique(dat.all$age_cat)
dat.all$age_cat[dat.all$age_cat=="Sub adult"] <- "Sub Adult"
dat.all$age_cat[dat.all$age_cat=="Mature Adult"] <- "Adult"

#and plot
p1 <- ggplot(data=dat.all) + geom_point(aes(x=forearm_mm, y=mass_g, color=age_cat, shape=sex)) + facet_grid(age_cat~.)
      
      

#fit separate regression to adult and sub-adult data to get residuals
m1 <- lm(mass_g~forearm_mm, data=subset(dat.all, age_cat=="Adult"))
summary(m1)

m2 <- lm(mass_g~forearm_mm, data=subset(dat.all, age_cat=="Sub Adult"))
summary(m2)

#and add prediction to the dataset
dat.all$predicted_mass_g <- NA
dat.all$predicted_mass_g[dat.all$age_cat=="Adult" &!is.na(dat.all$mass_g)] <- predict(m1)
dat.all$predicted_mass_g[dat.all$age_cat=="Sub Adult" &!is.na(dat.all$mass_g)] <- predict(m2)



#and plot
# This can be supplementary FigS3 of the paper
FigS3 <- ggplot(data=subset(dat.all, !is.na(age_cat) & age_cat!="Juvenile")) + 
      geom_point(aes(x=forearm_mm, y=mass_g, color=age_cat, shape=sex)) +
      facet_grid(age_cat~.) + ylab("mass (g)") + xlab("forearm (mm)") +
      geom_line(aes(x=forearm_mm, y=predicted_mass_g)) + theme_bw() +
      theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))

ggsave(file = paste0(homewd,"/supp-figures/figS3.png"),
       plot=FigS3,
       units="mm",  
       width=50, 
       height=60, 
       scale=2.8, 
       dpi=300)


#now calculate the residuals
dat.all$mass_residuals <- dat.all$mass_g-dat.all$predicted_mass_g


#and plot with the serological data
dat.all$serostatus[dat.all$serostatus==0] <- "no"
dat.all$serostatus[dat.all$serostatus==1] <- "yes"
dat.all$serostatus <- as.factor(dat.all$serostatus)

#all
p3 <- ggplot(data=dat.all) +
      geom_boxplot(aes(x=serostatus, y=mass_residuals, fill=serostatus)) +
      facet_wrap(~virus_genus)

#and overall - not really helpful
p4 <- ggplot(data=dat.all) +
  geom_boxplot(aes(x=serostatus, y=mass_residuals, fill=serostatus)) 
  #probably no different

#check do stats on p3
library(mgcv)
library(lme4)
dat.all$virus_genus <- as.factor(dat.all$virus_genus)

m3 <- lmer(mass_residuals~serostatus + (1|virus_genus), data = dat.all)
summary(m3) #not sig

dat.all$serostatus_num <- as.character(dat.all$serostatus )
dat.all$serostatus_num[dat.all$serostatus_num=="no"] <- 0
dat.all$serostatus_num[dat.all$serostatus_num=="yes"] <- 1
dat.all$serostatus_num <- as.numeric(dat.all$serostatus_num)
m3 <- glmer(serostatus~mass_residuals + (1|virus_genus), family="binomial", data = dat.all)
summary(m3) #not sig

#try as gam
m3gam <- gam(mass_residuals~s(serostatus, bs="re") +
                            s(virus_genus, bs="re"), data = dat.all)


m3gam <- gam(serostatus_num~s(mass_residuals, bs="tp") +
               s(virus_genus, bs="re"), data = dat.all)

summary(m3gam)
plot(m3gam)

#try some individual families
dat.sub = subset(dat.all, virus_genus=="Betacoronavirus" |
                          virus_genus=="Henipavirus" | 
                          virus_genus=="Mastadenovirus" |
                          virus_genus=="Lyssavirus" |
                          virus_genus=="Marburgvirus" |
                          virus_genus=="Flavivirus")#|
                          #virus_genus=="Ebolavirus") 
                          #removed b/c no mass resid for seropositives here

#and this could be a supplementary plot


dat.sub$serostatus[dat.sub$serostatus =="no"] 
dat.sub$serostatus <- as.character(dat.sub$serostatus)
dat.sub$serostatus[dat.sub$serostatus =="no"]  <- "seronegative"
dat.sub$serostatus[dat.sub$serostatus =="yes"]  <- "seropositive"
dat.sub$serostatus <- as.factor(dat.sub$serostatus)

#Fig S4
FigS4 <- ggplot(data=dat.sub) + theme_bw() +
  geom_boxplot(aes(x=serostatus, y=mass_residuals, fill=serostatus), show.legend = F) +
  facet_wrap(~virus_genus) + ylab("mass:forearm residuals") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        strip.background = element_rect(fill="white")) +
  geom_hline(data=dat.sub, aes(yintercept=0), linetype=2)

ggsave(file = paste0(homewd,"/supp-figures/figS4.png"),
       plot=FigS4,
       units="mm",  
       width=80, 
       height=50, 
       scale=2.8, 
       dpi=300)

#and stats by family

m4 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Betacoronavirus"))
summary(m4)


m5 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Henipavirus"))
summary(m5)


m6 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Lyssavirus"))
summary(m6)

m7 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Marburgvirus"))
summary(m7)

m8 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Mastadenovirus"))
summary(m8)

m9 <- glm(serostatus_num ~ mass_residuals, family="binomial", data=subset(dat.sub, virus_genus=="Flavivirus"))
summary(m9)
