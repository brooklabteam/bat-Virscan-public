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



#now convert back to long
dat.long <- melt(dat.wide)
names(dat.long) <- c("ID", "virus_genus", "serostatus")
nrow(dat.long) #462 = 77*6
unique(dat.long$serostatus)

#and merge with metadata
head(exp.dat)
head(dat.long)
merge.dat <- dplyr::select(exp.dat, -(virus_genus))
dat.all <- merge(dat.long, merge.dat, by=c("ID"), all.x = T)
head(dat.all)
unique(dat.all$serostatus)

#now get residuals
unique(dat.all$age_cat)

#check by condition
unique(dat.all$condition)
dat.all$condition <- factor(dat.all$condition, levels=c("Poor", "Fair", "Good", "Excellent"))

p1 <- ggplot(data=subset(dat.all, !is.na(condition))) + theme_bw() + 
      theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.text=element_text(size=14)) +
      geom_boxplot(aes(x=condition, y=mass_residuals, color=condition), show.legend = F) +
      geom_hline(aes(yintercept=0), linetype=2) + ylab("mass residuals") + xlab("\nbat body condition")
#great support that mass residuals actually tell us something about mortality



#and plot with the serological data


plot.dat <- subset(dat.all, !is.na(mass_residuals))
# plot.sum <- ddply(plot.dat, .(virus_genus), summarise, Npos = sum(serostatus), Ntot=length(serostatus))
# plot.sum$Nneg = plot.sum$Ntot-plot.sum$Npos
# 
# plot.sum = subset(plot.sum, Npos>0 & Nneg>0)
# plot.merge <- dplyr::select(plot.sum, virus_genus)
# plot.merge$keep = 1
# 
# plot.dat <- merge(plot.dat, plot.merge, by="virus_genus", all.x = T)
# plot.dat <- subset(plot.dat, keep ==1)

plot.dat$serostatus[plot.dat$serostatus==0] <- "seroneg"
plot.dat$serostatus[plot.dat$serostatus==1] <- "seropos"
plot.dat$serostatus <- factor(plot.dat$serostatus, levels=c("seropos", "seroneg"))


#all
FigS4 <- ggplot(data=plot.dat) + theme_bw()+
      theme(panel.grid = element_blank(), axis.title.x = element_blank(),
            legend.position = c(.9,.04), legend.title = element_blank(),
            legend.background = element_rect(color="black")) +
      geom_violin(aes(x=serostatus, y=mass_residuals, fill=serostatus), 
                  draw_quantiles = c(0,.25,.5,.75)) +
      geom_jitter(aes(x=serostatus, y=mass_residuals), width=.1, size=.1) +
      geom_hline(yintercept = 1, linetype=2) +
      facet_wrap(~virus_genus) + ylab("mass residuals")

ggsave(file = paste0(homewd,"/supp-figures/figS4.png"),
       plot=FigS4,
       units="mm",  
       width=100, 
       height=90, 
       scale=2.8, 
       dpi=300)



#check do stats on FigS4
library(mgcv)
library(sjPlot)
plot.dat$virus_genus <- as.factor(plot.dat$virus_genus)

plot.dat$serostatus_num <- as.character(plot.dat$serostatus )
plot.dat$serostatus_num[plot.dat$serostatus_num=="seroneg"] <- 0
plot.dat$serostatus_num[plot.dat$serostatus_num=="seropos"] <- 1
plot.dat$serostatus_num <- as.numeric(plot.dat$serostatus_num)

plot.dat$ID <- as.factor(plot.dat$ID)

#is poor body condition associated with positive serostatus?
m1 <- glmer(serostatus_num~mass_residuals + (1|virus_genus) + (1|ID), family="binomial", data = plot.dat)
summary(m1) #not sig but trends negative
plot_model(m1, "std") #here's the fixed effects visualized
plot_model(m1, "pred") #here's the marginal effects if how mass resid relates to serostatus
plot_model(m1, "re") #modulating effects by genera

#sig negative association between mass-residuals and serostatus


#try as gam as well
plot.dat$condition
plot.dat$sex <- as.factor(plot.dat$sex)
plot.dat$age_cat <- as.factor(plot.dat$age_cat)
plot.dat$forearm_mm

m1gam <- gam(serostatus_num~s(mass_residuals, bs="tp") +
               #s(forearm_mm, bs="tp") +
               #s(mass_g, bs="tp") +
               #s(condition, bs="re") +
               #s(age_cat, bs="re") +
               #s(sex, bs="re") +
               s(ID, bs="re") +
               s(virus_genus, bs="re"), 
              family = "binomial",
              data = plot.dat)

summary(m1gam)
source(paste0(homewd,"/prep-scripts/mollentze-streicker-2020-functions.R"))
genus.df <- get_partial_effects(fit=m1gam, var="virus_genus")
resid.df <- get_partial_effects_continuous(gamFit=m1gam, var="mass_residuals")

#partial effects of mass residuals on serostatus
#not sig but trends negative
plot.partial.cont(df=resid.df, var="mass_residuals", response_var = "serostatus", log=FALSE, alt_var = "mass:forearm residuals")
#and modulating effects by genus
plot.partial(df=genus.df, var="virus_genus", response_var = "serostatus")
# I trust above more than the glm. No clear sig interactions for a single virus.
# Likly we need much more data
length(unique(plot.dat$ID)) #only 59 bats here

# Mass residuals are not significantly associated serostatus,
# Slight trend of lower residuals trending towards higher
# serostatus but not significant. All effects are modulated
# by differences by individual and relationship varies
# significantly across disparate viral genera.

# This is shown by two statistical methods.

##########################################
##########################################

# Now look at predictors of (a) number of hits (total peptide hits)
# And (b) number of exposures

# These should all be summarized by individual

ind.dat <- ddply(dat.all, .(ID, rank, sex, condition, age_cat, age_tooth, mass_g, 
                             forearm_mm, mass_residuals, tot_hits, tot_filter_hits), summarise, N_exposures=sum(serostatus))
head(ind.dat)

#positive condition associated with more total hits
ggplot(subset(ind.dat, !is.na(condition))) + 
  geom_boxplot(aes(x=condition, y=tot_hits, color=condition),show.legend = F)

#filter hits just maps total hits
ggplot(subset(ind.dat, !is.na(condition))) + 
  geom_boxplot(aes(x=condition, y=tot_filter_hits, color=condition),show.legend = F) +
  geom_jitter(aes(x=condition, y=tot_filter_hits, color=condition), width = .1)

#filter hits just maps total hits
ggplot(subset(ind.dat, !is.na(condition))) + 
  geom_boxplot(aes(x=condition, y=N_exposures, color=condition),show.legend = F) +
  geom_jitter(aes(x=condition, y=N_exposures, color=condition), width = .1)



ind.dat$condition_simp <- NA
ind.dat$condition_simp[ind.dat$condition=="Poor" | ind.dat$condition=="Fair"] <- "neg"
ind.dat$condition_simp[ind.dat$condition=="Good" | ind.dat$condition=="Excellent"] <- "pos"

ggplot(subset(ind.dat, !is.na(condition_simp))) + 
  geom_boxplot(aes(x=condition_simp, y=tot_hits, color=condition_simp),show.legend = F)


ind.dat$age_cat <- factor(ind.dat$age_cat, levels=c("Juvenile", "Sub Adult", "Adult" ))

#definitely more hits in older bats
ggplot(subset(ind.dat, !is.na(age_cat))) + 
  geom_boxplot(aes(x=age_cat, y=tot_hits, color=age_cat),show.legend = F) +
  geom_jitter(aes(x=age_cat, y=tot_hits, color=age_cat),width=.1,show.legend = F) 

#nominal differences in exposures
#makes since b/c our exposures count is missing bat-specific exposures
ggplot(subset(ind.dat, !is.na(age_cat))) + 
  geom_boxplot(aes(x=age_cat, y=N_exposures, color=age_cat),show.legend = F)+
  geom_jitter(aes(x=age_cat, y=N_exposures, color=age_cat),width=.1,show.legend = F) 

#no differences by sex
ggplot(subset(ind.dat, !is.na(sex))) + 
  geom_boxplot(aes(x=sex, y=tot_hits, color=sex),show.legend = F)
ggplot(subset(ind.dat, !is.na(sex))) + 
  geom_boxplot(aes(x=sex, y=N_exposures, color=sex),show.legend = F)

#vaguely positice association with hits and tooth
ggplot(subset(ind.dat, !is.na(age_tooth))) + 
  geom_point(aes(x=age_tooth, y=tot_hits, color=age_cat),show.legend = F)

#falls apart with exposures
ggplot(subset(ind.dat, !is.na(age_tooth))) + 
  geom_point(aes(x=age_tooth, y=N_exposures, color=age_cat),show.legend = F)


#exposure by body size? yes
ggplot(subset(ind.dat, !is.na(mass_g))) + 
  geom_point(aes(x=mass_g, y=tot_hits, color=condition),show.legend = F)

ggplot(subset(ind.dat, !is.na(forearm_mm))) + 
  geom_point(aes(x=forearm_mm, y=tot_hits, color=condition),show.legend = F)


ggplot(subset(ind.dat, !is.na(mass_residuals))) + 
  geom_point(aes(x=mass_residuals, y=tot_hits, color=condition),show.legend = F)


ind.dat$MF_split <- NA
ind.dat$MF_split[ind.dat$mass_residuals<=0] <- "low"
ind.dat$MF_split[ind.dat$mass_residuals>1] <- "high"


ggplot(subset(ind.dat, !is.na(MF_split))) + 
  geom_violin(aes(x=MF_split, y=tot_hits, color=MF_split),
              show.legend = F, draw_quantiles = c(0,.25,.5,.75)) +
  geom_jitter(aes(x=MF_split, y=tot_hits, color=MF_split), width = .1)



ggplot(subset(ind.dat, !is.na(MF_split))) + 
  geom_violin(aes(x=MF_split, y=N_exposures, color=MF_split),
              show.legend = F, draw_quantiles = c(0,.25,.5,.75)) +
  geom_jitter(aes(x=MF_split, y=N_exposures, color=MF_split), width = .1)



## try some stats - predict tot hits - can't do it or N_exposures with any of these variables
length(ind.dat$age_tooth[!is.na(ind.dat$age_tooth)]) #40
length(ind.dat$age_cat[!is.na(ind.dat$age_cat)]) #63
length(ind.dat$sex[!is.na(ind.dat$sex)]) #72
length(ind.dat$mass_residuals[!is.na(ind.dat$mass_residuals)]) #59

ind.dat$age_tooth[ind.dat$age_cat=="Juvenile"] <- 1 # not working
ind.dat$age_tooth <- as.numeric(ind.dat$age_tooth)
ind.dat$sex <- as.factor(ind.dat$sex)
ind.dat$age_cat <- as.factor(ind.dat$age_cat)
ind.dat$condition <- as.factor(ind.dat$condition)

gam2 <- gam(tot_hits ~ s(age_tooth, bs="tp"), #s(age_cat, bs="re") +
                      #s(mass_g, bs="tp"),# +
                      #s(mass_residuals, bs="tp"),# +
                      #s(sex, bs="re"),# +
                      #s(condition, bs="re"), 
                      data = ind.dat)
summary(gam2) #nothing sig


gam3 <- gam(N_exposures ~ s(tot_hits, bs="tp"),
              #s(age_tooth, bs="tp"),#s(age_cat, bs="re") +
              #s(mass_g, bs="tp") +
              #s(mass_residuals, bs="tp") +
            #s(sex, bs="re"),# +
            #s(condition, bs="re"), 
            data = ind.dat)
summary(gam3) #nothing sig


ggplot(ind.dat,) + 
  geom_point(aes(x=tot_hits, y=N_exposures),show.legend = F)

#exposures is realted to total hits. this is predictive. 
#beyond that, can't tell what is the driver

