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
all.dat <- rbind(pa.dat, es.dat)

rank.bat = cbind.data.frame(ID = unique(all.dat$ID), rank = seq(1, length(unique(all.dat$ID))))
all.dat <- merge(all.dat, rank.bat, by="ID")

# rank.bat = cbind.data.frame(ID = unique(pa.dat$ID), rank = seq(1, length(unique(pa.dat$ID))))
# pa.dat <- merge(pa.dat, rank.bat, by="ID")
# 
# rank.bat2 = cbind.data.frame(ID = unique(es.dat$ID), rank = seq(1, length(unique(es.dat$ID))))
# es.dat <- merge(es.dat, rank.bat2, by="ID")


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

## ideally, would make nested color ramp with virus genus colors nested in virus subfamily and family colors

Fig2a <- ggplot(data=pa.family) + geom_bar(aes(x=virus_family, y=seroprev, fill=virus_family), stat="identity", position = "dodge") + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=8, vjust=-2,  color="black"), strip.text = element_text(size=14),
                      axis.text.y=element_text(size=12,color="black"),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)#+scale_fill_manual(values=cols1) 


Fig2b <- ggplot(data=pa.subfamily) + geom_bar(aes(x=virus_subfamily, y=seroprev, fill=virus_subfamily), stat="identity", position = "dodge")  + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7),  strip.background = element_rect(fill="white"), strip.placement = "outside",
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=8, vjust=-2,  color="black"), strip.text = element_text(size=14),
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      axis.text.y=element_text(size=12,color="black"),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)#+scale_fill_manual(values=cols1) 

Fig2c <- ggplot(data=pa.genus) + geom_bar(aes(x=virus_genus, y=seroprev, fill=virus_genus), stat="identity", position = "dodge")  + facet_grid(~cat)+
  theme_bw() +  theme(legend.title = element_blank(), legend.text =element_text(size=7), strip.background = element_rect(fill="white"), strip.placement = "outside",
                      panel.grid = element_blank(), axis.text.x = element_text(angle=300, size=8, vjust=-2,  color="black"), strip.text = element_text(size=14),
                      plot.margin = unit(c(.2,.8,.1,.2), "cm"),
                      axis.text.y=element_text(size=12,color="black"),axis.title.y=element_text(size=14),axis.title.x=element_text(size=14))+
  labs(x="", y="seroprevalence\n")+theme(legend.position="none")+ylim(0,0.35)#+scale_fill_manual(values=cols1) 

Fig2abc <- cowplot::plot_grid(Fig2a, Fig2b, Fig2c, ncol=1, nrow=3, labels = c("A", "B", "C"), label_size = 20)


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
                     draw_quantiles=c(0.025,0.5,0.975), show.legend = F) + ylab("number virus exposures") +
         geom_jitter(aes(x=species, y=N_exposures),size=1,width=.1,height=0) +
        theme(panel.grid = element_blank(), axis.text.x = element_text(size=14, face="italic", color="black"),
              plot.margin = unit(c(.1,.1,1.5,.1), "cm"),
        axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_blank())

Fig2d2 <- ggplot(data=bat.ID.sum) + theme_bw() +
  geom_boxplot(aes(x=species, y=N_exposures, fill=species), scale = "width",
              draw_quantiles=c(0.025,0.5,0.975), show.legend = F) + ylab("number virus exposures") +
  #geom_jitter(aes(x=species, y=N_exposures),size=1,width=.1,height=0) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(size=14, face="italic", color="black"),
        plot.margin = unit(c(.1,.1,1.5,.1), "cm"),
        axis.text.y=element_text(size=14,color="black"),axis.title.y=element_text(size=16),axis.title.x=element_blank())


#bat.ID.sum$species <- factor(bat.ID.sum$species, levels =rev(c("Eonycteris spelaea", "Pteropus alecto", "Homo sapiens")))


Fig2d3 <- ggplot(bat.ID.sum, aes(x=N_exposures, fill=species)) + 
          geom_histogram(binwidth=2, color="black") + theme_bw() +
          theme(legend.title = element_blank(), legend.text=element_text(face="italic"), panel.grid = element_blank(),
                axis.title = element_text(size=16), axis.text = element_text(size=14),
                plot.margin = unit(c(.2,.1,1.5,.1), "cm"),
                legend.position = c(.85,.85)) + xlab("number of virus exposures") + ylab("frequency") +
          scale_y_continuous(breaks = c(0,20,40,60,80,100,120))

#compare
sum.dat <- ddply(bat.ID.sum, .(species), summarise, mean=mean(N_exposures), median = median(N_exposures))


#and plot
Fig2def <- cowplot::plot_grid(Fig2d3, Fig2d1, ncol=1, nrow=3, labels = c("D", "E"), label_size = 20)



Fig2 <- cowplot::plot_grid(Fig2abc, Fig2def , ncol=2, nrow=1, rel_widths = c(1,.6))


ggsave(file = paste0(homewd,"/final-figures/fig2.png"),
       plot=Fig2,
       units="mm",  
       width=90, 
       height=90, 
       scale=4.5, 
       dpi=300)

#or the boxplot version

Fig2def <- cowplot::plot_grid(Fig2d3, Fig2d2, ncol=1, nrow=3, labels = c("D", "E"), label_size = 20)

Fig2 <- cowplot::plot_grid(Fig2abc, Fig2def , ncol=2, nrow=1, rel_widths = c(1,.6))
ggsave(file = paste0(homewd,"/final-figures/fig2box.png"),
       plot=Fig2,
       units="mm",  
       width=90, 
       height=90, 
       scale=4.5, 
       dpi=300)



########
########
# and repeat in the supplement for eonycteris