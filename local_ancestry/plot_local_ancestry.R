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
populations


d <- read.table("summary.all.txt")
names(d)=c("ref", "pop", "local", "mean", "+/-", "stdev")
head(d)


## just focus on grch38 for now

hg38 = d %>% filter(ref=="grch38")
hg38_sup = inner_join(hg38, populations, by=c("pop"="Population_code"))
head(hg38_sup)

## tally EUR populations

hg38_sup_eur = hg38_sup %>% filter(Superpopulation_code=="EUR")
hg38_sup_eur

eur<-ggplot(hg38_sup_eur, aes(x=local, y=mean, fill=local)) + 
     geom_bar(stat='identity') + facet_grid(pop ~ .) + 
     theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
     theme(legend.position = "none") +
     ggtitle("EUR Samples")

eur

## tally AFR populations 

hg38_sup_afr = hg38_sup %>% filter(Superpopulation_code=="AFR")
hg38_sup_afr

afr<-ggplot(hg38_sup_afr, aes(x=local, y=mean, fill=local)) + 
  geom_bar(stat='identity') + facet_grid(pop ~ .) + 
  theme(axis.text.x = element_text(angle = 90)) + geom_text(aes(label=mean), vjust=+5) + 
  ggtitle("AFR Samples")

afr


## make a composite plot
grid.arrange(eur, afr, nrow=1, widths=c(1,1.25))



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
