library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)
library(Hmisc)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/local_ancestry")


## Load the different populations/superpopulations
populations=read.table("../population-variant-stats-per-sample/populations.txt", header=TRUE)
head(populations)

## load the summary data for all populations
d <- read.table("summary.all.txt")
names(d)=c("ref", "pop", "local", "mean", "+/-", "stdev")
dd = inner_join(d, populations, by=c("pop"="Population_code"))
head(dd)



## just focus on grch38 for now
###############################################################################
hg38 = dd %>% filter(ref=="grch38")
head(hg38)

## tally EUR populations

hg38_eur = hg38 %>% filter(Superpopulation_code=="EUR")
hg38_eur

eur<-ggplot(hg38_eur, aes(x=local, y=mean, fill=local)) + 
     geom_bar(stat='identity') + facet_grid(pop ~ .) + 
     theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
     theme(legend.position = "none") +
     ggtitle("EUR Samples/GRCh38")

eur

## tally AFR populations 

hg38_afr = hg38 %>% filter(Superpopulation_code=="AFR")
hg38_afr

afr<-ggplot(hg38_afr, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  ggtitle("AFR Samples/GRCh38")

afr


## tally AMR populations 

hg38_amr = hg38 %>% filter(Superpopulation_code=="AMR")
hg38_amr

amr<-ggplot(hg38_amr, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  theme(legend.position = "none") +
  ggtitle("AMR Samples/GRCh38")

amr


## tally EAS populations 

hg38_eas = hg38 %>% filter(Superpopulation_code=="EAS")
hg38_eas

eas<-ggplot(hg38_eas, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  theme(legend.position = "none") +
  ggtitle("EAS Samples/GRCh38")

eas

## tally SAS populations 

hg38_sas = hg38 %>% filter(Superpopulation_code=="SAS")
hg38_sas

sas<-ggplot(hg38_sas, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  theme(legend.position = "none") +
  ggtitle("SAS Samples/GRCh38")

sas



## make a composite plot
grid.arrange(eur, afr, nrow=1, widths=c(1,1.25))


grid.arrange(amr, eas, eur, sas,  afr, nrow=1, widths=c(1,1,1,1,1.25))


## close look at CEU

grch38_ceu <- read.table("data/grch38.CEU.ancestry.txt")
names(grch38_ceu) = c("chr", "start", "end", "cnt", "pop")

ggplot(grch38_ceu, aes(x=pop, y=cnt, fill=pop)) + geom_violin() + stat_summary(fun.y=mean, geom="point", size=2, color="black") +
  coord_cartesian(ylim=c(0,20)) + ggtitle("CEU/GRCh38 variant density")


grch38_ceu_afr = grch38_ceu %>% filter(pop=="AFR")
grch38_ceu_eur = grch38_ceu %>% filter(pop=="EUR")

t.test(grch38_ceu_afr$cnt, grch38_ceu_eur$cnt)
summary(grch38_ceu_afr$cnt)
summary(grch38_ceu_eur$cnt)





## Now CHM13
###############################################################################

chm13 = dd %>% filter(ref=="chm13")
head(chm13)

## tally EUR populations

chm13_eur = chm13 %>% filter(Superpopulation_code=="EUR")
chm13_eur

eur<-ggplot(chm13_eur, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  theme(legend.position = "none") +
  ggtitle("EUR Samples/CHM13")

eur

## tally AFR populations 

chm13_afr = chm13 %>% filter(Superpopulation_code=="AFR")
chm13_afr

afr<-ggplot(chm13_afr, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  ggtitle("AFR Samples/CHM13")

afr


## make a composite plot
grid.arrange(eur, afr, nrow=1, widths=c(1,1.25))

