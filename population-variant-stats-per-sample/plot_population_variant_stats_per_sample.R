library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/population-variant-stats-per-sample")


## Summarize the Superpopulation/population map, only needs to be done once
###############################################################################

#stats <- as.data.frame(t(read.table("../samtools_stats_chm13/2021.04.22.samtools.stats.all.txt", header=TRUE, row.names=1)))
#stats$Population_code      = as.factor(stats$Population_code)
#stats$Superpopulation_code = as.factor(stats$Superpopulation_code)
#populations = data.frame(Population_code=stats$Population_code, Superpopulation_code=stats$Superpopulation_code)
#populations = populations %>% distinct() %>% arrange(Population_code) %>% arrange(Superpopulation_code)
#populations
#write.table(populations, "populations.txt", row.names=FALSE, quote=FALSE)

## Load the different populations/superpopulations
populations=read.table("populations.txt", header=TRUE)
populations


## Load the genome/population variant counts
###############################################################################
cnts = data.frame()

for (genome in c("chm13", "hg38"))
{
  print(genome)
  for (pop in populations$Population_code)
  {
    filename = paste0(genome,"/population/",pop,".pop.psc")
    print(paste("loading", filename))
    tt <- read.table(filename, header=TRUE)

    tt_long = tt %>% 
      rename(sample=X.3.sample, nNonRefHom=X.5.nNonRefHom, nHets=X.6.nHets, nIndels=X.9.nIndels, nHapAlt=X.13.nHapAlt) %>%
      select(chr, sample, nNonRefHom, nHets, nIndels, nHapAlt) %>%
      pivot_longer(cols=c(nNonRefHom, nHets, nIndels, nHapAlt), names_to="type", values_to="count")
    
    tt_long$genome = genome
    tt_long$population = pop
    
    cnts <- rbind(cnts, tt_long)
  }
}

cnts$population = factor(cnts$population, populations$Population_code)
cnts$genome = factor(cnts$genome, levels=c("hg38", "chm13"), labels=c("GRCh38", "CHM13"))


## Plot all
###############################################################################

cntall <- cnts %>% group_by(sample, genome, population) %>% summarize(variants=sum(count))
cntall = inner_join(cntall, populations, by=c("population"="Population_code"))
cntall$population = factor(cntall$population, populations$Population_code)

## by individual population
plotall = ggplot(cntall, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~population) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

plotall

## by superpopulation
plotall_super = ggplot(cntall, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~Superpopulation_code) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

plotall_super


## Just hets
###############################################################################

cnthet <- cnts %>% filter(type=="nHets") %>% group_by(sample, genome, population) %>% summarize(variants=sum(count))
cnthet = inner_join(cnthet, populations, by=c("population"="Population_code"))
cnthet$population = factor(cntall$population, populations$Population_code)

plothet = ggplot(cnthet, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~population) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Heterozygous Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

plothet

plothet_super = ggplot(cnthet, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~Superpopulation_code) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Heterozygous Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

plothet_super



## Just homozygous
###############################################################################

cnthom <- cnts %>% filter(type=="nNonRefHom") %>% group_by(sample, genome, population) %>% summarize(variants=sum(count))
cnthom = inner_join(cnthom, populations, by=c("population"="Population_code"))
cnthom$population = factor(cntall$population, populations$Population_code)

plothom = ggplot(cnthom, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~population) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Homozygous Variant Counts") + theme(plot.title = element_text(hjust = 0.5)) 

plothom

plothom_super = ggplot(cnthom, aes(x=genome, y=variants)) + geom_boxplot(aes(fill=Superpopulation_code)) + facet_grid(~Superpopulation_code) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Homozygous Variant Counts") + theme(plot.title = element_text(hjust = 0.5)) 

plothom_super



## Final plot with all panels together
###############################################################################

grid.arrange(plotall + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()),
             plothet + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()),
             plothom + theme(legend.position="bottom", legend.title=element_blank()) + guides(colour = guide_legend(nrow = 1)), 
             ncol=1,heights=c(1,1,1.4))


grid.arrange(plotall_super + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()),
             plothet_super + theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()),
             plothom_super + theme(legend.position="bottom", legend.title=element_blank()) + guides(colour = guide_legend(nrow = 1)), 
             ncol=1,heights=c(1,1,1.4))