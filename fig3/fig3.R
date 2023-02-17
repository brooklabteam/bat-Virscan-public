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
library(sjPlot)


#load the data
#comment out when not for your system
homewd = "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/bat-VirScan-public"
#homewd = "/Users/carabrook/Developer/bat-Virscan-public"
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

#now, slim to not repeat the same individuals
dat.mass <- dplyr::select(dat.all, ID, rank, condition, mass_residuals)
dat.mass <- dat.mass[!duplicated(dat.mass),]
dat.mass <- subset(dat.mass, !is.na(condition))

#great support that mass residuals actually tell us something about body condition
Fig3a1 <- ggplot(data=dat.mass) + theme_bw() + 
      theme(panel.grid = element_blank(), axis.title = element_text(size=14), axis.text=element_text(size=12),
            plot.margin = unit(c(.5,.1,2.2,.1), "cm")) +
            #plot.margin = unit(c(.1,.1,.1,.1), "cm")) +
      geom_boxplot(aes(x=condition, y=mass_residuals, color=condition), show.legend = F) +
      geom_jitter(aes(x=condition, y=mass_residuals, color=condition), width = .1, show.legend = F) +
      geom_hline(aes(yintercept=0), linetype=2) + ylab("mass : forearm residuals") + xlab("bat body condition")


m1 <- lm(mass_residuals ~ condition,data = dat.mass)
summary(m1)
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          -85.12      31.83  -2.674  0.00985 ** 
#   conditionFair         41.41      34.87   1.187  0.24018    
# conditionGood         87.52      33.10   2.644  0.01065 *  
#   conditionExcellent   250.00      42.11   5.936 2.04e-07 ***


#is it supported statistically? (yes)
Fig3a <- Fig3a1 + 
        annotate("text", label = "**",x=1, y=210, size = 8)+
        annotate("text", label = "*",x=3, y=210, size = 8)+
        annotate("text", label = "***",x=4, y=210, size = 8) 

# (note that I also ran this with a GAM, which showed the same results - but 
# seems plenty clean and simple to just report the glm)

####################################################################################
####################################################################################

# Now, let's look at the relationship between mass residuals and serostatus,
# as a marker of the impact of host "health" on virus exposure

plot.dat <- subset(dat.all, !is.na(mass_residuals))
plot.dat$serostatus_num <- plot.dat$serostatus
plot.dat$serostatus[plot.dat$serostatus==0] <- "seroneg"
plot.dat$serostatus[plot.dat$serostatus==1] <- "seropos"
plot.dat$serostatus <- factor(plot.dat$serostatus, levels=c("seropos", "seroneg"))


plot.dat$virus_genus <- as.factor(plot.dat$virus_genus)
plot.dat$serostatus_num <- as.numeric(plot.dat$serostatus_num)
plot.dat$age_tooth<- as.numeric(plot.dat$age_tooth)
plot.dat$ID <- as.factor(plot.dat$ID)


### First, for the SuppMat, plot all

# add which were significant as stars
plot.dat$sig_label <- ""
plot.dat$sig_label[plot.dat$virus_genus=="Alphavirus"| plot.dat$virus_genus=="Enterovirus" | 
                   plot.dat$virus_genus== "Gammaretrovirus" | plot.dat$virus_genus=="Betainfluenzavirus" | plot.dat$virus_genus=="Mastadenovirus" | 
                   plot.dat$virus_genus=="Orthopoxvirus"| plot.dat$virus_genus=="Rotavirus" | plot.dat$virus_genus=="Rubivirus"] <- "***"


# plot all in the SuppMat
FigS4 <- ggplot(data=plot.dat) + theme_bw()+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        legend.position = c(.9,.04), legend.title = element_blank(),
        legend.background = element_rect(color="black")) +
  geom_violin(aes(x=serostatus, y=mass_residuals, fill=serostatus), 
              draw_quantiles = c(0,.25,.5,.75)) +
  geom_label(aes(x=1.5, y=170, label=sig_label), size=5, label.size = 0) +
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



# Next, assess this relationship statistically
# We ask: is poor body condition associated with positive serostatus?
# And is this relationship different for different viruses? 

#skipping version here with no random effects
#m2 <-  glm(serostatus_num~mass_residuals:virus_genus, family = "binomial", data = plot.dat)
#summary(m2)

# version here with random effects
# need the nAGQ=0 to avoid convergence issues.
# This was fine, but we appear to have lost the confidence
# intervals on the interaction terms
library(lme4)
m2 <- glmer(serostatus_num~ mass_residuals:virus_genus + (1|ID), nAGQ=0, family="binomial", plot.dat)
summary(m2)

#now wdjust for multiple comparisons by bonferonni:
adjusted.P <- p.adjust(summary(m2)$coefficients[,ncol(summary(m2)$coefficients)][-1], n= length(summary(m2)$coefficients[,ncol(summary(m2)$coefficients)][-1]))

adjusted.P <- cbind.data.frame(virus=sort(as.character(unique(plot.dat$virus_genus))), estimate=round(summary(m2)$coefficients[,1][-1],2), raw_P=summary(m2)$coefficients[,ncol(summary(m2)$coefficients)][-1], adj_P=adjusted.P)
rownames(adjusted.P) <- c()
adjusted.P$isSigRaw <- "no"
adjusted.P$isSigRaw[adjusted.P$raw_P<0.01] <- "yes"
adjusted.P$isSigCorrect <- "no"
adjusted.P$isSigCorrect[adjusted.P$adj_P<0.01] <- "yes"

#and look at those that are significant
subset(adjusted.P, isSigCorrect=="yes")

## This is your subset for the main text

#and Fig 3C side by side

# plot.sub.df = subset(plot.dat, virus_genus =="Gammaretrovirus" | virus_genus == "Alphavirus"| virus_genus =="Avulavirus" | 
#                        virus_genus == "Parechovirus" | virus_genus=="Betacoronavirus" | virus_genus == "Rotavirus" | 
#                        virus_genus == "Orthopoxvirus"| virus_genus == "Flavivirus" | virus_genus =="Lyssavirus" )
# 
# plot.sub.df = subset(plot.dat,virus_genus == "Alphavirus"| virus_genus =="Avulavirus" | 
#                        virus_genus == "Orthohantavirus" | virus_genus=="Betacoronavirus" | virus_genus == "Rotavirus" | 
#                        virus_genus == "Orthopoxvirus"| virus_genus == "Influenzavirus_B" | virus_genus =="Lyssavirus" )


plot.sub.df = subset(plot.dat,virus_genus == "Alphavirus"| virus_genus =="Enterovirus" | 
                       virus_genus == "Gammaretrovirus" | virus_genus=="Betainfluenzavirus" | virus_genus == "Mastadenovirus" | 
                       virus_genus == "Orthopoxvirus"| virus_genus == "Rotavirus" | virus_genus =="Rubivirus" )



plot.sub.df$interaction <- "negative slope"
#plot.sub.df$interaction[plot.sub.df$virus_genus=="Betacoronavirus" | plot.sub.df$virus_genus=="Orthohantavirus" | plot.sub.df$virus_genus=="Lyssavirus"] <- "positive slope"
plot.sub.df$interaction[plot.sub.df$virus_genus=="Enterovirus" | plot.sub.df$virus_genus=="Rubivirus"] <- "positive slope"


Fig3b <- ggplot(plot.sub.df, aes(x=virus_genus, y=mass_residuals)) + theme_bw()+
  theme(legend.position = c(.08,.9), legend.title = element_blank(),
        legend.background = element_rect(color="black"))+
  geom_violin(aes(fill=serostatus), 
              draw_quantiles = c(0,.25,.50,.75,1), show.legend = T) + 
  geom_hline(yintercept = 0, linetype=2) + ylab("mass : forearm residuals") +
  geom_point(aes(fill=serostatus), position = position_jitterdodge(), 
             size=.2, show.legend = F) + ggh4x::facet_nested(~interaction, scales= "free", space = "free") +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"),
        axis.title.x = element_blank(), axis.title.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.y = element_text(size=12), axis.text.x = element_text(size=12, angle = 290, hjust = 0),
        plot.margin = unit(c(.5,.1,.1,.1), "cm"))



# Note that we also tried including age (tooth) in this analysis but it was not 
# significant, so we ultimately just investigated the interaction of virus genus 
# and mass residual

#visualize interactions quickly
plot_model(m2, "int", mdrt.values="minmax")
plot_model(m2, "int", mdrt.values="meansd")

#and the random effects
plot_model(m2, "re")

plot_model(m2, "pred", terms=c("mass_residuals", "virus_genus"))

#and extract important features
int.dat <- plot_model(m2, "int", mdrt.values="meansd", ci.lvl = .9)$data
int.dat <- plot_model(m2, "int")$data

#this gets you the same
#int.dat.2 <- plot_model(m2, "pred", terms=c("mass_residuals", "virus_genus"), ci.lvl=.9)$data

int.df <- cbind.data.frame(int.dat)
head(int.df)

#now add which were significant as stars
int.df$sig_label <- ""
int.df$sig_label[int.df$group=="Alphavirus"| int.df$group=="Enterovirus" | 
                   int.df$group== "Gammaretrovirus" | int.df$group=="Betainfluenzavirus" | int.df$group=="Mastadenovirus" | 
                   int.df$group=="Orthopoxvirus"| int.df$group=="Rotavirus" | int.df$group=="Rubivirus"] <- "***"



#test code for adding back in confidence intervals - taken from summary(m2)
int.df$conf.low <- 0
int.df$conf.high <- 0
int.df$conf.low[int.df$group=="Alphavirus"] <- (int.df$predicted[int.df$group=="Alphavirus"] - (0.0032452))
int.df$conf.low[int.df$group=="Enterovirus"] <- (int.df$predicted[int.df$group=="Enterovirus"] - (0.0038048))
int.df$conf.low[int.df$group=="Gammaretrovirus"] <- (int.df$predicted[int.df$group=="Gammaretrovirus"] - (0.0036252))
int.df$conf.low[int.df$group=="Betainfluenzavirus"] <- (int.df$predicted[int.df$group=="Betainfluenzavirus"] - (0.0035187))
int.df$conf.low[int.df$group=="Mastadenovirus"] <- (int.df$predicted[int.df$group=="Mastadenovirus"] - (0.0035884))
int.df$conf.low[int.df$group=="Orthopoxvirus"] <- (int.df$predicted[int.df$group=="Orthopoxvirus"] - (0.0032728))
int.df$conf.low[int.df$group=="Rotavirus"] <- (int.df$predicted[int.df$group=="Rotavirus"] - (0.0033401))
int.df$conf.low[int.df$group=="Rubivirus"] <- (int.df$predicted[int.df$group=="Rubivirus"] - (0.0040620))

int.df$conf.high[int.df$group=="Alphavirus"] <- (int.df$predicted[int.df$group=="Alphavirus"] + (0.0032452))
int.df$conf.high[int.df$group=="Enterovirus"] <- (int.df$predicted[int.df$group=="Enterovirus"] + (0.0038048))
int.df$conf.high[int.df$group=="Gammaretrovirus"] <- (int.df$predicted[int.df$group=="Gammaretrovirus"] + (0.0036252))
int.df$conf.high[int.df$group=="Betainfluenzavirus"] <- (int.df$predicted[int.df$group=="Betainfluenzavirus"] + (0.0035187))
int.df$conf.high[int.df$group=="Mastadenovirus"] <- (int.df$predicted[int.df$group=="Mastadenovirus"] + (0.0035884))
int.df$conf.high[int.df$group=="Orthopoxvirus"] <- (int.df$predicted[int.df$group=="Orthopoxvirus"] + (0.0032728))
int.df$conf.high[int.df$group=="Rotavirus"] <- (int.df$predicted[int.df$group=="Rotavirus"] + (0.0033401))
int.df$conf.high[int.df$group=="Rubivirus"] <- (int.df$predicted[int.df$group=="Rubivirus"] + (0.0040620))


#Also for the Supp Matt, plot it all
FigS5 <- ggplot(data=int.df) + 
         geom_line(aes(x=x, y=predicted, color=group), show.legend = F)+
         #geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=.3, show.legend = F)+
         facet_wrap(group~.) + ylab("serostatus") + xlab("mass : forearm residuals") +
         scale_y_continuous(breaks=c(0,1)) + geom_vline(xintercept = 0, linetype=2) +
         coord_cartesian(ylim=c(0,1.1)) + theme_bw() + 
         geom_label(aes(x=0,y=.85, label=sig_label), label.size = 0, size=5) +
         theme(panel.grid = element_blank(), strip.background = element_blank(), 
               strip.text = element_text(size=12), axis.title = element_text(size=12), axis.text = element_text(size=9))


ggsave(file = paste0(homewd,"/supp-figures/figS5.png"),
       plot=FigS5,
       units="mm",  
       width=120, 
       height=90, 
       scale=2.8, 
       dpi=300)

#and subset to just a few examples that are significant for the main text

#these are the *** sig genera
# sub.df = subset(int.df, group =="Gammaretrovirus" | group == "Alphavirus"| group =="Avulavirus" | 
#                   group == "Parechovirus" | group=="Betacoronavirus" | group == "Rotavirus" | 
#                   group == "Orthopoxvirus"| group == "Flavivirus" | group =="Lyssavirus" )



sub.df = subset(int.df,group == "Alphavirus"| group =="Enterovirus" | 
                       group == "Gammaretrovirus" | group=="Betainfluenzavirus" | group == "Mastadenovirus" | 
                       group == "Orthopoxvirus"| group == "Rotavirus" | group =="Rubivirus" )

sub.df.test = subset(int.df,group == "Alphavirus")

# sub.df = subset(int.df, group == "Alphavirus"| group =="Avulavirus" | 
#                   group == "Orthohantavirus" | group=="Betacoronavirus" | group == "Rotavirus" | 
#                   group == "Orthopoxvirus"| group == "Influenzavirus_B" | group =="Lyssavirus" )


# @ Emily, not sure if you can find a way to get the confidence intervals back...?
# Fig3cflip <- ggplot(data=sub.df) + theme_bw() +
#   theme(panel.grid = element_blank(), legend.title = element_blank(),
#         legend.position = "bottom",legend.background = element_rect(color="black"),
#         legend.direction = "horizontal",
#         axis.title = element_text(size=16), axis.text = element_text(size=14),
#         plot.margin = unit(c(.5,.1,.1,.1), "cm")) +
#   geom_line(aes(x=x, y=predicted, color=group))+ 
#   ylab("serostatus") + xlab("mass : forearm residual") +
#   scale_y_continuous(breaks=c(0,1)) + geom_vline(xintercept = 0, linetype=2) +
#   coord_flip(ylim=c(0,1.1)) #+ geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=.1)

# @ Emily, not sure if you can find a way to get the confidence intervals back...?
Fig3c <- ggplot(data=sub.df, aes(x=x, y=predicted, color=group)) + theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = "bottom",legend.background = element_rect(color="black"),
        legend.direction = "horizontal",
        axis.title = element_text(size=16), axis.text = element_text(size=14),
        plot.margin = unit(c(.5,.1,.1,.1), "cm")) +
  geom_line()+
  ylab("serostatus") + xlab("mass : forearm residual") +
  scale_y_continuous(breaks=c(0,1)) + geom_vline(xintercept = 0, linetype=2) +
  coord_cartesian(ylim=c(0,1.1)) + 
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=.5, color=NA)

#test code doesn't work
#plot(predicted ~ x, data=sub.df.test, ylim=c(0,1.1),type='n')
#polygon(c(sub.df.test$x,rev(sub.df.test$x)),c(sub.df.test$conf.low,rev(sub.df.test$conf.high)),col = "grey75", border = FALSE)
#plotCI(x=sub.df.test$x, y=sub.df.test$predicted, li=sub.df.test$conf.low, ui=sub.df.test$conf.high)
#Fig3c + geom_errorbar(aes(x=x, ymin=conf.low, ymax=conf.high), color="gray")


#and combine these together:

# Fig3flip <- cowplot::plot_grid(Fig3a, Fig3b, Fig3cflip, ncol=3, labels = c("A", "B", "C"), rel_widths = c(1,1.2,1), label_size = 22)
# 
# ggsave(file = paste0(homewd,"/final-figures/fig3flip.png"),
#        plot=Fig3flip,
#        units="mm",  
#        width=150, 
#        height=50, 
#        scale=3, 
#        dpi=300)

Fig3 <- cowplot::plot_grid(Fig3a, Fig3b, Fig3c, ncol=3, labels = c("A", "B", "C"), rel_widths = c(1,1.2,1), label_size = 22)

ggsave(file = paste0(homewd,"/final-figures/fig3.png"),
       plot=Fig3,
       units="mm",  
       width=150, 
       height=50, 
       scale=3, 
       dpi=300)

