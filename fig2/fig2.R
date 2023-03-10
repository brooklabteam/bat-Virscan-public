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
library(paletteer)


#for fig2, report the seroprevalence of all viruses based on these rules:

## if an individual bat is seropositive to multiple subfamilies in the same
## viral family at once, then we only report to family, as this is likely cross-reactivity

#if an individual bat is seropositive to multiple genera in the same viral subfamily at once,
#then we only report to subfamily, as this is likely cross-reactivity

#if an individual bat is seropositive to multiple species in the same viral genus
#at once, then we only report to genus, as this is likely cross-reactivity


########################################################################
########################################################################

#set wd - make a link to homewd for your own folder
homewd = "/Users/carabrook/Developer/bat-VirScan-public"
homewd = "/Users/emilyruhs/Desktop/UChi_Brook_Lab/GitHub_repos/bat-VirScan-public"

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

dat$Age_tooth[dat$Age_cat=="Juvenile"] <- 0

length(unique(dat$ID[dat$Species=="Pteropus alecto"])) #77 P. alecto. Some are missing ()
length(unique(dat$ID[dat$Species=="Eonycteris spelaea"])) #5


dat$Age_cat[dat$Age_cat=="Mature Adult"] <- "Adult"

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

head(dat)
head(meta)

meta$Age_cat[meta$Age_cat=="Mature Adult"] <- "Adult"
meta$Age_cat[meta$Age_cat=="Sub adult"] <- "Sub Adult"
meta$Age_cat <- factor(meta$Age_cat, levels=rev(c("Juvenile", "Sub Adult", "Adult")))
ggplot(data=subset(meta, Species=="Pteropus alecto" & !is.na(Age_cat))) + geom_point(aes(x=Forearm_mm, y=Mass_g, color=Age_cat, shape=Sex)) + facet_grid(Age_cat~.)

#fit separate regression to adult and sub-adult data to get residuals
m1 <- lm(Mass_g~Forearm_mm, data=subset(meta, Age_cat=="Adult" & Species=="Pteropus alecto"))
summary(m1)

sink("m1_adult_massforearm.txt")
summary(m1)
sink()

m2 <- lm(Mass_g~Forearm_mm, data=subset(meta, Age_cat=="Sub Adult"& Species=="Pteropus alecto"))
summary(m2)

sink("m2_subadult_massforearm.txt")
summary(m2)
sink()

#and add prediction to the dataset
meta$predicted_mass_g <- NA
meta$predicted_mass_g[meta$Age_cat=="Adult" &!is.na(meta$Mass_g) & meta$Species=="Pteropus alecto"] <- predict(m1)
meta$predicted_mass_g[meta$Age_cat=="Sub Adult" &!is.na(meta$Mass_g)& meta$Species=="Pteropus alecto"] <- predict(m2)



#and plot
# This can be supplementary FigS3 of the paper
# This can be supplementary FigS3 of the paper
library(ggtext)
ann_text <- data.frame(Age_cat = c("Adult",  "Sub Adult" ), 
                       xpos = c(155,155), ypos = c(955,955), 
                       lab = c("0.386","0.453"))

FigS3 <- ggplot(data=subset(meta, !is.na(Age_cat) & Age_cat!="Juvenile" & Species=="Pteropus alecto")) + 
  geom_point(aes(x=Forearm_mm, y=Mass_g, color=Age_cat, shape=Sex)) +
  facet_grid(Age_cat~.) + ylab("mass (g)") + xlab("forearm (mm)") +
  geom_line(aes(x=Forearm_mm, y=predicted_mass_g)) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"))+
  geom_richtext(data = ann_text, size=3, vjust=0, hjust=0,aes(x = xpos,  y = ypos, label = paste0("adj R<sup>2</sup> = ", lab))
  )


ggsave(file = paste0(homewd,"/supp-figures/figS3_annotated.png"),
       plot=FigS3,
       units="mm",  
       width=50, 
       height=60, 
       scale=2.8, 
       dpi=300)


#now calculate the residuals
meta$mass_residuals <- meta$Mass_g-meta$predicted_mass_g

setdiff(unique(dat$ID), unique(meta$ID)) #none. all phipseq outputs are represented in the metadata
setdiff(unique(meta$ID), unique(dat$ID)) 

# we are missing phip-seq data for these bats
#"Pa29"  "Pa50"  "Pa59"  "Pa77"  "Pa81"  "Pa10"  "Pa100" "Pa101" "Pa1"   "Pa2"   "Pa3"  

#want to merge in the body condition scores
names(meta)
meta.merge <- dplyr::select(meta, ID, Sex, Age_cat, Age_scale, Condition, Mass_g, Forearm_mm, mass_residuals)

names(dat)
dat <- dplyr::select(dat, -(MF_resid), -(Sex), -(Age_cat), -(Mass), -(Forearm), -(Age_tooth))
dat <- merge(dat, meta.merge, by="ID", all.x = T)
head(dat)
dat <- arrange(dat, virus_subfamily, virus_genus, virus_species)
dat$virus_species <- factor(dat$virus_species, levels=unique(dat$virus_species))
head(dat)

unique(dat$X) # these only have values on the eonycteris samples, where the columns appear to be shifted
# fix down below


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


length(unique(dat$virus_subfamily)) #44 subfamilies with nested colors within them

# first, for figure 1, plot a heatmap and compare broad viral exposures for the two species
# include a supplementary figure that is focused on just one or two viral families (maybe just paramyxoviruses)
es.dat = subset(dat, Species=="Eonycteris spelaea") # why are the assigned counts so few here?
head(es.dat) # it appears that the column names are shifted to the right on these guys about halfway through
names(es.dat)[11:23] <- c("X", "Assigned.counts", "Assigned.peptides", "Total.sample.hits", "Total.filtered.sample.hits", "virus_subfamily", "Sex", "Age_cat", "Age_scale", "Condition","Mass_g", "Forearm_mm",  "mass_residuals")
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
all.dat <- rbind(pa.dat, es.dat)

rank.bat = cbind.data.frame(ID = unique(all.dat$ID), rank = seq(1, length(unique(all.dat$ID))))
all.dat <- merge(all.dat, rank.bat, by="ID")

names(all.dat)[16:21] <- c("sex", "age_cat", "age_tooth", "condition", "mass_g", "forearm_mm")
head(all.dat)
#all.dat$Age_tooth[all.dat$Age_cat=="Juvenile"] <- 0
#plot body condition by number of peptides


########################################################################
########################################################################

#now split by individual bat and score the seropositivity for each

all.bat.id <- dlply(all.dat, .(rank))

score.virus <- function(df1){
  n.viral.subfam <- length(unique(df1$virus_subfamily))
  
  ## if an individual bat is seropositive to multiple subfamilies in the same
  ## viral family at once, then we only report to family, as this is likely cross-reactivity
  
  if (n.viral.subfam>1){
    out.df <- cbind.data.frame(rank = unique(df1$rank), ID = unique(df1$ID), 
                               virus_family = unique(df1$virus_family), 
                               virus_subfamily = NA, virus_genus=NA,
                               cat="virus_family")
  }else{
    
    #if an individual bat is seropositive to multiple genera in the same viral subfamily at once,
    #then we only report to subfamily, as this is likely cross-reactivity
    
    n.viral.genera <- length(unique(df1$virus_genus))
    
    if (n.viral.genera>1){
      out.df <- cbind.data.frame(rank = unique(df1$rank), ID = unique(df1$ID), 
                                 virus_family = unique(df1$virus_family),  
                                 virus_subfamily =  unique(df1$virus_subfamily), 
                                 virus_genus=NA,
                                 cat="virus_subfamily")
    }else{
      
      #if an individual bat is seropositive to multiple species in the same viral genus
      #at once, then we only report to genus, as this is likely cross-reactivity
      
      ## once we get to this level, we are scoring only to genus regardless
      
      n.viral.species <- length(unique(df1$virus_species))
      
      out.df <- cbind.data.frame(rank = unique(df1$rank), ID = unique(df1$ID),
                                 virus_family = unique(df1$virus_family),  
                                 virus_subfamily =  unique(df1$virus_subfamily), 
                                 virus_genus=  unique(df1$virus_genus), 
                                 cat="virus_genus")
    }}
  
  return(out.df)
}
score.sero <- function(df){
  
  #now split by viral family and apply
  df.split <- dlply(df, .(virus_family))
  
  #and apply
  df.split.out <- lapply(df.split, score.virus)
  
  out.df <- data.table::rbindlist(df.split.out) 
  #nrow(out.df) is the number of exposures for this bat
  
  out.df$bat_species <- unique(df$Species)
  
  return(out.df)
}

#and apply
all.bat.sero <- lapply(all.bat.id, score.sero)

all.bat.df <- data.table::rbindlist(all.bat.sero)

head(all.bat.df)

#and bind the total hits by bat id
head(all.dat)
dat.add <- dplyr::select(all.dat, ID, sex, age_cat, age_tooth, condition,  mass_g, forearm_mm, mass_residuals, Total.sample.hits, Total.filtered.sample.hits)
names(dat.add) <- c("ID", "sex", "age_cat", "age_tooth", "condition", "mass_g", "forearm_mm", "mass_residuals", "tot_hits", "tot_filter_hits")
dat.add <- dat.add[!duplicated(dat.add),]
all.bat.df <- merge(all.bat.df, dat.add, by="ID", all.x = T)
head(all.bat.df)

#check and clean as needed
unique(all.bat.df$age_cat)
unique(all.bat.df$condition)
all.bat.df$age_tooth[all.bat.df$age_cat=="Juvenile"] #one juvenile bat is scored as age 3
all.bat.df[all.bat.df$age_cat=="Juvenile",] 
#he has forearm == 145. this is erroneous. 
#we are going to reclass him as 1 year like all the other juveniles, as well as the juvenile who is labelled as NA for age
all.bat.df$age_tooth[all.bat.df$ID=="Pa96"] <- 1
all.bat.df$age_tooth[all.bat.df$ID=="Pa37"] <- 1

head(all.bat.df)


#and save this - for use in Fig 3 and 4
write.csv(all.bat.df, file = paste0(homewd, "/working-data/all_bat_exposures.csv"), row.names = F)

length(unique(all.bat.df$rank[all.bat.df$bat_species=="Pteropus alecto"])) #77
length(unique(all.bat.df$rank[all.bat.df$bat_species=="Eonycteris spelaea"])) #5

#and summarise - bv virus family, subfamily, genus
sum.family <- ddply(all.bat.df, .(bat_species, virus_family), summarise, N=length(unique(rank)))
sum.subfamily <- ddply(subset(all.bat.df, !is.na(virus_subfamily)), .(bat_species, virus_family, virus_subfamily), summarise, N=length(unique(rank)))
sum.genus <- ddply(subset(all.bat.df, !is.na(virus_genus)), .(bat_species, virus_family, virus_subfamily, virus_genus), summarise, N=length(unique(rank)))

#and seroprev
sum.family$seroprev <- NA
sum.family$seroprev[sum.family$bat_species=="Pteropus alecto"] <- sum.family$N[sum.family$bat_species=="Pteropus alecto"]/77
sum.family$seroprev[sum.family$bat_species=="Eonycteris spelaea"] <- sum.family$N[sum.family$bat_species=="Eonycteris spelaea"]/5


sum.subfamily$seroprev <- NA
sum.subfamily$seroprev[sum.subfamily$bat_species=="Pteropus alecto"] <- sum.subfamily$N[sum.subfamily$bat_species=="Pteropus alecto"]/77
sum.subfamily$seroprev[sum.subfamily$bat_species=="Eonycteris spelaea"] <- sum.subfamily$N[sum.subfamily$bat_species=="Eonycteris spelaea"]/5


sum.genus$seroprev <- NA
sum.genus$seroprev[sum.genus$bat_species=="Pteropus alecto"] <- sum.genus$N[sum.genus$bat_species=="Pteropus alecto"]/77
sum.genus$seroprev[sum.genus$bat_species=="Eonycteris spelaea"] <- sum.genus$N[sum.genus$bat_species=="Eonycteris spelaea"]/5


#and plot them all
pa.family <- subset(sum.family, bat_species=="Pteropus alecto")
pa.subfamily <- subset(sum.subfamily, bat_species=="Pteropus alecto")
pa.genus <- subset(sum.genus, bat_species=="Pteropus alecto")

es.family <- subset(sum.family, bat_species=="Eonycteris spelaea")
es.subfamily <- subset(sum.subfamily, bat_species=="Eonycteris spelaea")
es.genus <- subset(sum.genus, bat_species=="Eonycteris spelaea")


pa.family <- arrange(pa.family, bat_species, desc(seroprev))
pa.family$virus_family <- factor(pa.family$virus_family, levels = c(unique(pa.family$virus_family)))

pa.subfamily <- arrange(pa.subfamily, bat_species, desc(seroprev))
pa.subfamily$virus_subfamily <- factor(pa.subfamily$virus_subfamily, levels = c(unique(pa.subfamily$virus_subfamily)))

pa.genus <- arrange(pa.genus, bat_species, desc(seroprev))
pa.genus$virus_genus <- factor(pa.genus$virus_genus, levels = c(unique(pa.genus$virus_genus)))

pa.family$cat <- "virus family"
pa.subfamily$cat <- "virus subfamily"
pa.genus$cat <- "virus genus"


es.family <- arrange(es.family, bat_species, desc(seroprev))
es.family$virus_family <- factor(es.family$virus_family, levels = c(unique(es.family$virus_family)))

es.subfamily <- arrange(es.subfamily, bat_species, desc(seroprev))
es.subfamily$virus_subfamily <- factor(es.subfamily$virus_subfamily, levels = c(unique(es.subfamily$virus_subfamily)))

es.genus <- arrange(es.genus, bat_species, desc(seroprev))
es.genus$virus_family <- factor(es.genus$virus_family, levels = c(unique(es.genus$virus_family)))


es.family$cat <- "virus family"
es.subfamily$cat <- "virus subfamily"
es.genus$cat <- "virus genus"

#are there any families within the genus and subfamily datasets that are skipped?
setdiff(unique(pa.genus$virus_family), unique(pa.family$virus_family)) #no
setdiff(unique(pa.subfamily$virus_family), unique(pa.family$virus_family)) #no
setdiff(unique(pa.genus$virus_subfamily), unique(pa.subfamily$virus_subfamily)) #no

setdiff(unique(es.genus$virus_family), unique(es.family$virus_family)) #no
setdiff(unique(es.subfamily$virus_family), unique(es.family$virus_family)) #no
setdiff(unique(es.genus$virus_subfamily), unique(es.subfamily$virus_subfamily)) #no


# and summarize the data to report in the paper
length(unique(pa.family$virus_family)) #33
length(unique(pa.subfamily$virus_subfamily)) #41
length(unique(pa.genus$virus_genus)) #57


length(unique(es.family$virus_family)) #12
length(unique(es.subfamily$virus_subfamily)) #14
length(unique(es.genus$virus_genus)) #9

sort(as.character(unique(pa.genus$virus_genus)))

# [1] "Alphainfluenzavirus" "Alphatorquevirus"    "Alphavirus"          "Avulavirus"          "Betacoronavirus"    
# [6] "Betainfluenzavirus"  "Betapapillomavirus"  "Betapolyomavirus"    "Circovirus"          "Deltavirus"         
# [11] "Ebolavirus"          "Enterovirus"         "Erythroparvovirus"   "Flavivirus"          "Gammainfluenzavirus"
# [16] "Gammaretrovirus"     "Henipavirus"         "Hepacivirus"         "Hepatovirus"         "Husavirus"          
# [21] "Lentivirus"          "Lymphocryptovirus"   "Lyssavirus"          "Mammarenavirus"      "Marburgvirus"       
# [26] "Mastadenovirus"      "Metapneumovirus"     "Molluscipoxvirus"    "Norovirus"           "Orthobunyavirus"    
# [31] "Orthohantavirus"     "Orthohepadnavirus"   "Orthohepevirus"      "Orthonairovirus"     "Orthopneumovirus"   
# [36] "Orthopoxvirus"       "Orthoreovirus"       "Parapoxvirus"        "Parechovirus"        "Pegivirus"          
# [41] "Phlebovirus"         "Picobirnavirus"      "Respirovirus"        "Roseolovirus"        "Rotavirus"          
# [46] "Rubivirus"           "Salivirus"           "Sapovirus"           "Seadornavirus"       "Simplexvirus"       
# [51] "Spumavirus"          "Tl2011virus"         "unc. Arena"          "unc. Cyclo"          "unc. Papilloma"     
# [56] "unc. Polyoma"        "Vesiculovirus" 

#sort(as.character(unique(es.genus$virus_genus)))
# [1] "Alphavirus"        "Avulavirus"        "Enterovirus"       "Marburgvirus"      "Orthohepadnavirus"
# [6] "Orthopneumovirus"  "Sapovirus"         "Simplexvirus"      "Spumavirus"     

## ideally, would make nested color ramp with virus genus colors nested in virus subfamily and family colors

Fig2a <- ggplot(data=pa.family) + geom_bar(aes(x=virus_family, y=seroprev, fill=virus_family), stat="identity", position = "dodge") + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=9, vjust=-2,  color="black"), strip.text = element_text(size=16),
                      axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)+scale_fill_grey()


Fig2b <- ggplot(data=pa.subfamily) + geom_bar(aes(x=virus_subfamily, y=seroprev, fill=virus_subfamily), stat="identity", position = "dodge")  + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7),  strip.background = element_rect(fill="white"), strip.placement = "outside",
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=9, vjust=-2,  color="black"), strip.text = element_text(size=14),
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)+scale_fill_grey()

Fig2c <- ggplot(data=pa.genus) + geom_bar(aes(x=virus_genus, y=seroprev, fill=virus_genus), stat="identity", position = "dodge")  + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=9, vjust=-2,  color="black"), strip.text = element_text(size=16),
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)+scale_fill_grey()

#Fig2abc <- cowplot::plot_grid(Fig2a, Fig2b, Fig2c, ncol=1, nrow=3, labels = c("A", "B", "C"), label_size = 20)
Fig2ac <- cowplot::plot_grid(Fig2a, Fig2c,ncol=1, nrow=2, labels = c("A", "B"), label_size = 20)


head(pa.genus)
# make color ramp to nest the genera within the subfamily
# should just be able to order them and assign a color ramp

# first, get a data table of just the virus families and their genera
pa.df1 <- dplyr::select(pa.family, virus_family)
pa.df2 <- dplyr::select(pa.genus, virus_family)
pa.sum <- ddply(pa.genus, .(virus_family, virus_genus), summarise)
head(pa.sum)
pa.sum$virus_genus <- as.character(pa.sum$virus_genus)
pa.sum$virus_family <- as.character(pa.sum$virus_family)

# are there families not captured here?

setdiff(unique(pa.sum$virus_family), unique(pa.family$virus_family))
setdiff(unique(pa.family$virus_family), unique(pa.sum$virus_family)) 
#"Astroviridae" - needs to be added at the end

pa.add <- c("Astroviridae", NA)
pa.sum <-rbind(pa.sum, pa.add)

# split by family
pa.split <- dlply(pa.sum, .(virus_family))

# and make a list of palettes in R that can be subsampled
# add in a few that you made on your own
col.pals <- as.list(rep(c("Blues", "Reds", "Purples", "Oranges", "Greens", "Greys", "Turquoises", "LightPinks", "Yellows", "Seagreens", "Magentas", "Cyans", "Chartreuses", "SlateBlues"), length.out=length(pa.split)))

# within, each family, collapse into a vector with family name
# first (only once), followed by the names of different genera
nested.colors <- function(df, col.pal){
  
  fam.name = unique(df$virus_family)
  #print(fam.name)
  gen.names = c(df$virus_genus)
  ordered.vector <- c(fam.name, gen.names)
  ordered.vector <- ordered.vector[!is.na(ordered.vector)]
  
  #if(length(ordered.vector)>2){
  # colorz = rev(brewer.pal(n=length(ordered.vector), name = col.pal))  
  #}else{
  brewer.palettes <- c("Blues", "Reds", "Purples", "Oranges", "Greens", "Greys")
  brewer.pal = brewer.palettes[brewer.palettes==col.pal]
  if (length(brewer.pal)>0){
    colorz = rev(brewer.pal(n=4, name = col.pal))  
    colorz <- colorz[1:length(ordered.vector)]
  }else{
    #here are the manual color palettes
    LightPinks <- rev(c("lightpink", "lightpink1", "lightpink2", "lightpink3"))
    Yellows <- rev(c("yellow1", "yellow2", "yellow3", "yellow4"))
    Seagreens <- rev(c("seagreen1", "seagreen2", "seagreen3", "seagreen4"))
    Turquoises <- rev(c("turquoise1", "turquoise2", "turquoise3", "turquoise4"))
    Magentas <- rev(c("violetred1", "violetred2", "violetred3", "violetred4"))
    Cyans <- rev(c("cyan1", "cyan2", "cyan3", "cyan4"))
    Chartreuses <- rev(c("chartreuse1", "chartreuse2", "chartreuse3", "chartreuse4"))
    SlateBlues <- rev(c("slateblue1", "slateblue2", "slateblue3", "slateblue4"))
    
    palette.list <- list(LightPinks, Yellows, Seagreens, Turquoises, Magentas, Cyans, Chartreuses, SlateBlues)
    names(palette.list) <- c("LightPinks",  "Yellows", "Seagreens", "Turquoises", "Magentas", "Cyans", "Chartreuses", "SlateBlues")
    palette = palette.list[[col.pal]]
    
    colorz <- palette[1:length(ordered.vector)]
    
  }
  
  #}
  
  names(colorz) <- ordered.vector
  
  
  
  return(colorz)
}

# or you can use this function to give the same color to the genera and families:
same.nested.colors <- function(df, col.pal.num){
  
  fam.name = unique(df$virus_family)
  #print(fam.name)
  gen.names = c(df$virus_genus)
  ordered.vector <- c(fam.name, gen.names)
  ordered.vector <- ordered.vector[!is.na(ordered.vector)]
  
  all.colors = paletteer_c("grDevices::rainbow", 33)
  
  
  colorz <- rep(all.colors[col.pal.num], length(ordered.vector))
  
  names(colorz) <- ordered.vector  
  
  return(colorz)
}

#here's the output with nested colors - you could add more manual palettes
pa.out <- mapply(nested.colors, df=pa.split, col.pal=col.pals, SIMPLIFY = FALSE)

#and here's the one for the same colors
pa.out <- mapply(same.nested.colors, df=pa.split, col.pal.num=as.list(1:33), SIMPLIFY = FALSE)
names(pa.out) <- c()
colorz <- c(unlist(pa.out))

# now you have a vector of colors in rank order, named by 
# nested virus family and genera

# rerun the plots (leaving off the subfamily plot)

# first, virus family:
Fig2a <- ggplot(data=pa.family) + geom_bar(aes(x=virus_family, y=seroprev, fill=virus_family), stat="identity", position = "dodge") + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=9, vjust=-2,  color="black"), strip.text = element_text(size=16),
                      axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16))+
  scale_fill_manual(values=colorz) +
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)#+scale_fill_manual(values=cols1) 

#then, go straight to genera:
Fig2b <- ggplot(data=pa.genus) + geom_bar(aes(x=virus_genus, y=seroprev, fill=virus_genus), stat="identity", position = "dodge")  + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=9, vjust=-1,  color="black"), strip.text = element_text(size=16),
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_text(size=16))+
  scale_fill_manual(values=colorz) +
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)#+scale_fill_manual(values=cols1) 

Fig2ab <- cowplot::plot_grid(Fig2a, Fig2b,ncol=1, nrow=2, labels = c("A", "B"), label_size = 20)


#and compare average virus exposures per bat

head(all.bat.df) # each line per ID is a different exposure

bat.ID.sum <- ddply(all.bat.df, .(bat_species, rank, ID), summarise, N_exposures=length(ID))
names(bat.ID.sum)[names(bat.ID.sum)=="bat_species"] <- "species"
#and the elledge data
ell.dat <- c(rep(0,12), rep(2,34), rep(4,46), rep(6,77), rep(8,112), rep(10,94), rep(12,72), rep(14,55), rep(16,30), rep(18, 14), rep(20,17), rep(22,3), rep(24,3), rep(26,1))
hist(ell.dat)
hum.dat <- cbind.data.frame(species="Homo sapiens", rank=NA, ID=NA, N_exposures=ell.dat)
head(bat.ID.sum)
bat.ID.sum <- rbind(bat.ID.sum, hum.dat)
tail(bat.ID.sum)
bat.ID.sum$species <- factor(bat.ID.sum$species, levels =rev(c("Eonycteris spelaea", "Pteropus alecto", "Homo sapiens")))

#test elledge data
test<- ggplot(subset(bat.ID.sum, species=="Homo sapiens"), aes(x=N_exposures, fill=species)) + 
  geom_histogram(binwidth=2, color="black") + theme_bw() +
  theme(legend.title = element_blank(), legend.text=element_text(face="italic"),
        legend.position = c(.85,.85)) + xlab("number of virus exposures") + ylab("frequency") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100,120)) + scale_x_continuous(breaks=c(0,5,10,15,20,25,30))
#this replicates the plot from Xu et al. 2015



Fig2d1 <- ggplot(data=bat.ID.sum) + theme_bw() +
  geom_violin(aes(x=species, y=N_exposures, fill=species), scale = "width",
              draw_quantiles=c(0.025,0.5,0.975), show.legend = F) + ylab("number viruses detected") +
  geom_jitter(aes(x=species, y=N_exposures),size=1,width=.1,height=0) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size=14, face="italic", color="black"),
        plot.margin = unit(c(.1,.1,1.5,.1), "cm"),
        axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_blank())+scale_fill_grey()

# Fig2d2 <- ggplot(data=bat.ID.sum) + theme_bw() +
#   geom_boxplot(aes(x=species, y=N_exposures, fill=species), scale = "width",
#                draw_quantiles=c(0.025,0.5,0.975), show.legend = F) + ylab("number viruses detected") +
#   #geom_jitter(aes(x=species, y=N_exposures),size=1,width=.1,height=0) +
#   theme(panel.grid = element_blank(), axis.text.x = element_text(size=14, face="italic", color="black"),
#         plot.margin = unit(c(.1,.1,1.5,.1), "cm"),
#         axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_blank())


#bat.ID.sum$species <- factor(bat.ID.sum$species, levels =rev(c("Eonycteris spelaea", "Pteropus alecto", "Homo sapiens")))


Fig2d3 <- ggplot(bat.ID.sum, aes(x=N_exposures, fill=species)) + 
  geom_histogram(binwidth=2, color="black") + theme_bw() +
  theme(legend.title = element_blank(), legend.text=element_text(face="italic"), panel.grid = element_blank(),
        axis.title = element_text(size=16), axis.text = element_text(size=14),
        plot.margin = unit(c(.2,.1,1.5,.1), "cm"),
        legend.position = c(.85,.85)) + xlab("number of viruses detected") + ylab("frequency") +
  scale_y_continuous(breaks = c(0,20,40,60,80,100,120))+scale_fill_grey()

#compare
sum.dat <- ddply(bat.ID.sum, .(species), summarise, mean=mean(N_exposures), median = median(N_exposures))


#and plot
#Fig2def <- cowplot::plot_grid(Fig2d3, Fig2d1, ncol=1, nrow=3, labels = c("D", "E"), label_size = 20)
Fig2cd <- cowplot::plot_grid(Fig2d3, Fig2d1, ncol=1, nrow=2, labels = c("C", "D"), label_size = 20)



#Fig2 <- cowplot::plot_grid(Fig2abc, Fig2def , ncol=2, nrow=1, rel_widths = c(1,.6))
Fig2 <- cowplot::plot_grid(Fig2ac, Fig2cd , ncol=2, nrow=1, rel_widths = c(1,.6))


ggsave(file = paste0(homewd,"/final-figures/fig2_blackwhite.png"),
       plot=Fig2,
       units="mm",  
       width=90, 
       height=70, 
       scale=4.5, 
       dpi=300)

# #or the boxplot version
# 
# Fig2def <- cowplot::plot_grid(Fig2d3, Fig2d2, ncol=1, nrow=3, labels = c("D", "E"), label_size = 20)
# 
# Fig2 <- cowplot::plot_grid(Fig2abc, Fig2def , ncol=2, nrow=1, rel_widths = c(1,.6))
# ggsave(file = paste0(homewd,"/final-figures/fig2box.png"),
#        plot=Fig2,
#        units="mm",  
#        width=90, 
#        height=90, 
#        scale=4.5, 
#        dpi=300)


