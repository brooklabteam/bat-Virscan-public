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
homewd = "/Users/carabrook/Developer/bat-Virscan-public"
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


# Finally, look at predictors of (a) number of hits (total peptide hits)
# And (b) number of exposures

# These should all be summarized by individual

ind.dat <- ddply(dat.all, .(ID, rank, sex, condition, age_cat, age_tooth, mass_g, 
                            forearm_mm, mass_residuals, tot_hits, tot_filter_hits), summarise, N_exposures=sum(serostatus))
head(ind.dat)
ind.dat$age_cat <- factor(ind.dat$age_cat, levels = c("Juvenile", "Sub Adult", "Adult"))

# Age cat is merely pulled from tooth age:
p1 <- ggplot(data=subset(ind.dat, !is.na(age_cat))) + 
  geom_boxplot(aes(x=age_cat, y=age_tooth, fill=age_cat)) +
  geom_jitter(aes(x=age_cat, y=age_tooth), width = .1)
# Since tooth age is more quantitative, we'll use it as a predictor


# What are the predictors of total hits?

# This old stuff is incorrect -- skip over commented out
# # Is mass residual a predictor?
# m3 <- lm(log10(tot_hits)~mass_residuals, data = ind.dat)
# summary(m3) #nope
# 
# 
# # Age?
# m4 <- lm(log10(tot_hits)~age_tooth, data = ind.dat)
# summary(m4) # yes, weakly
# 
# m4alt <- glm(tot_hits~age_tooth, data = ind.dat, family="poisson")
# summary(m4alt) # yes, positive, strongly
# plot_model(m4alt, type="pred")

# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.48639    0.06142  40.484   <2e-16 ***
#   age_tooth    0.02234    0.01228   1.819   0.0768 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2046 on 38 degrees of freedom
# (37 observations deleted due to missingness)
# Multiple R-squared:  0.08008,	Adjusted R-squared:  0.05587 
# F-statistic: 3.308 on 1 and 38 DF,  p-value: 0.07683

#gam shows the same curve
#m4 <- gam(log10(tot_hits)~s(age_tooth, bs="tp"), data = ind.dat)
#summary(m4) # yes, weakly
# 
# ind.dat$N_exposures_scale <- ind.dat$N_exposures + 1
# m5gam <- gam(log10(N_exposures_scale)~s(age_tooth, bs="tp"), data = ind.dat)
# summary(m5gam) #not sig
# 
# m5lm <- lm(log10(N_exposures_scale)~age_tooth, data = ind.dat)
# summary(m5lm) #nope
# 
# m5alt <- glm(N_exposures~age_tooth, data = ind.dat, family="poisson")
# summary(m5alt) #weak positive association
# plot_model(m5alt, type="pred")
# 
# 
# m6 <- lm(log10(N_exposures_scale)~mass_residuals, data = ind.dat)
# summary(m6) #nope
# m6alt <- glm(N_exposures~mass_residuals, data = ind.dat, family="poisson")
# summary(m6alt) # yes, negative slope
# plot_model(m6alt, type="pred") #fewer exposures in healthier bats


# Look at age and mass residuals together as predictors of total hits and 
# overall viral exposures :

m7a <- glm(tot_hits~mass_residuals + age_tooth, data = ind.dat, family="poisson")
summary(m7a) #both sig but mass resid negative now

p2 <- plot_model(m7a, type="pred")$mass_residuals
p1 <- plot_model(m7a, type="pred")$age_tooth

m7b <- glm(N_exposures~mass_residuals + age_tooth, data = ind.dat, family="poisson")
summary(m7b) #both sig. mirrors above

p4 <-  plot_model(m7b, type="pred")$mass_residuals
p3 <- plot_model(m7b, type="pred")$age_tooth

cowplot::plot_grid(p1,p2,p3,p4, nrow=2, ncol=2)


# And pull the data out and put into ggplot- this is the new figure 4.

dat1 <- p1$data
head(dat1)
dat1$response <- "total peptide hits"
dat1$pred <- "age"
dat1$label <- "A"

dat1$label_x = 1
dat1$label_y = 1050

dat2 <- p2$data
head(dat2)
dat2$response <- "total peptide hits"
dat2$pred <- "mass : forearm residuals"
dat2$label <- "C"

dat2$label_x = -140
dat2$label_y = 1050



dat3 <- p3$data
head(dat3)
dat3$response <- "viral exposures"
dat3$pred <- "age"
dat3$label <- "B"

dat3$label_x = 1
dat3$label_y = 34

dat4 <- p4$data
head(dat4)
dat4$response <-  "viral exposures"
dat4$pred <- "mass : forearm residuals"
dat4$label <- "D"
dat4$label_x = -140
dat4$label_y = 34

# join the datasets
combine.df <- rbind(dat1,dat2,dat3,dat4)
head(combine.df)

unique(combine.df$response)
unique(combine.df$pred)

#and summarise data to include
library(reshape2)
ind.slim <- dplyr::select(ind.dat, ID, age_tooth, mass_residuals, tot_hits, N_exposures)
head(ind.slim)
ind.sum.age <- melt(ind.slim, id.vars = c("ID", "age_tooth"))
head(ind.sum.age)
ind.sum.age = subset(ind.sum.age, variable!="mass_residuals")
ind.sum.mass <- melt(ind.slim, id.vars = c("ID", "mass_residuals"))
ind.sum.mass = subset(ind.sum.mass, variable!="age_tooth")
names(ind.sum.age) <- c("ID", "x", "response", "y")
head(ind.sum.age)
ind.sum.age$pred <- "age"
ind.sum.age$response <-  as.character(ind.sum.age$response)
ind.sum.age$response[ind.sum.age$response=="tot_hits"] <- "total peptide hits"
ind.sum.age$response[ind.sum.age$response=="N_exposures"] <- "viral exposures"

head(ind.sum.mass)
names(ind.sum.mass) <- c("ID", "x", "response", "y")
head(ind.sum.mass)
ind.sum.mass$pred <- "mass : forearm residuals"
ind.sum.mass$response <-  as.character(ind.sum.mass$response)
ind.sum.mass$response[ind.sum.mass$response=="tot_hits"] <- "total peptide hits"
ind.sum.mass$response[ind.sum.mass$response=="N_exposures"] <- "viral exposures"

head(ind.sum.mass)

ind.sum <- rbind(ind.sum.age, ind.sum.mass)
head(ind.sum)
max(ind.sum$y[ind.sum$response=="viral exposures"]) #88
max(combine.df$predicted[combine.df$response=="viral exposures"]) #16

ind.sum$flag=0
ind.sum$flag[ind.sum$response=="viral exposures" & ind.sum$y>40] <- 1
ind.sum$flag[ind.sum$pred=="mass : forearm residuals" & ind.sum$x< -150] <- -1
ind.sum <- subset(ind.sum, flag==0)
#plot with ggplot:

Fig4 <- ggplot(data=combine.df) + 
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high), alpha=.3) +
  geom_line(aes(x=x, y=predicted)) + 
  geom_label(aes(x=label_x, y=label_y, label=label), label.size = 0, size=6, fontface="bold") +
  geom_point(data=ind.sum, (aes(x=x, y=y)), color="gray40") +
  facet_grid(response~pred, scales="free", switch = "both") + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        axis.text = element_text(size=12),
        strip.placement = "outside")


ggsave(file = paste0(homewd,"/final-figures/fig4.png"),
       plot=Fig4,
       units="mm",  
       width=70, 
       height=60, 
       scale=3, 
       dpi=300)

# 
# 
# #include the first model in the plot
# ind.dat$predicted_hits <- ind.dat$predicted_hits_lci <- ind.dat$predicted_hits_uci <-NA
# ind.dat$predicted_hits[!is.na(ind.dat$age_tooth)] <- 10^predict(m4)
# ind.dat$predicted_hits_lci <- ind.dat$predicted_hits_uci <- NA
# 
# ind.dat$predicted_hits_lci[!is.na(ind.dat$age_tooth)] <- 10^(predict(m4) -1.96*predict(m4, type = "response", se.fit = T)$se)
# ind.dat$predicted_hits_uci[!is.na(ind.dat$age_tooth)] <- 10^(predict(m4) +1.96*predict(m4, type = "response", se.fit = T)$se)
# # plotting shows us evidence of a weak association between
# # tot_hits and age
# max(ind.dat$tot_hits[!is.na(ind.dat$age_tooth)]) #872
# 
# Fig3b <- ggplot(subset(ind.dat, !is.na(age_tooth))) + 
#   geom_point(aes(x=age_tooth, y=tot_hits, color=age_cat), size=3) +
#   geom_line(aes(x=age_tooth, y= predicted_hits)) +
#   geom_ribbon(aes(x=age_tooth, ymin= predicted_hits_lci,  ymax= predicted_hits_uci), alpha=.3) +
#   coord_cartesian(ylim=c(0,900), xlim=c(0,13)) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_text(size=14), 
#         plot.margin = unit(c(.1,.1,.1,.1), "cm"),
#         axis.text = element_text(size=12), legend.position = c(.85,.15),
#         legend.title = element_blank(), legend.background = element_rect(color="black")) +
#   ylab("total peptide hits") + xlab("age (yrs, by teeth)")
# 
# 
# 
# Fig3top<- cowplot::plot_grid(Fig3a, Fig3b, labels = c("A", "B"), label_size = 20, nrow=1, ncol=2)
# Fig3bottom <- cowplot::plot_grid(Fig3c, Fig3d, labels = c("C", "D"), label_size = 20, nrow=1, ncol=2)
# 
# dev.new()
# Fig3 <- cowplot::plot_grid(Fig3top, Fig3bottom, ncol=1, nrow = 2, rel_heights = c(1,1.2))
# Fig3
# 
# ggsave(file = paste0(homewd,"/final-figures/fig3.png"),
#        plot=Fig3,
#        units="mm",  
#        width=90, 
#        height=80, 
#        scale=3, 
#        dpi=300)
# 
# 
# #and, as supplementary:
# 
# max(ind.dat$N_exposures[!is.na(ind.dat$age_tooth)]) #35
# 
# FigS6 <- ggplot(subset(ind.dat, !is.na(age_tooth))) + 
#   geom_point(aes(x=age_tooth, y=N_exposures, color=age_cat), size=3) +
#   #geom_line(aes(x=age_tooth, y= predicted_hits)) +
#   #geom_ribbon(aes(x=age_tooth, ymin= predicted_hits_lci,  ymax= predicted_hits_uci), alpha=.3) +
#   coord_cartesian(ylim=c(0,40), xlim=c(0,13)) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_text(size=14), 
#         plot.margin = unit(c(.1,.1,.1,.1), "cm"),
#         axis.text = element_text(size=12), legend.position = c(.87,.15),
#         legend.title = element_blank(), legend.background = element_rect(color="black")) +
#   ylab("total exposures") + xlab("age (yrs, by teeth)")
# 
# 
# 
# ggsave(file = paste0(homewd,"/supp-figures/figS6.png"),
#        plot=FigS6,
#        units="mm",  
#        width=50, 
#        height=40, 
#        scale=3, 
#        dpi=300)
