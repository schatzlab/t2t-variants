library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/population-variant-stats")

## Summarize the populations, only needs to be done once
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
cnts = data.frame()

for (genome in c("chm13", "hg38"))
{
  print(genome)
  for (pop in populations$Population_code)
  {
    filename = paste0(genome,"/population/",pop,".chr.stats")
    print(paste("loading", filename))
    tt <- read.table(filename, header=TRUE)
    
    tt$chr = factor(tt$chr, tt$chr)
    tt_long = pivot_longer(tt, cols=2:3, names_to="type", values_to="count")
    tt_long$genome = genome
    tt_long$population = pop
    
    cnts <- rbind(cnts, tt_long)
  }
}

cnts
cnts$population = factor(cnts$population, populations$Population_code)


## Plot chr1

chr1 <- cnts %>% filter(chr=="chr1")
chr1

pchr1 = ggplot(chr1, aes(x=genome, y=count)) + 
  geom_col(aes(fill=type)) + facet_grid(~population) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Chr1 Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

pchr1


## Summarize the whole genome

cnts %>% filter(genome=="chm13", population=="MXL") %>% summarize(sum(count))

cntall <- cnts %>% group_by(genome, population, type) %>% summarize(variants = sum(count))
cntall %>% filter(population=="MXL")

pgenome = ggplot(cntall, aes(x=genome, y=variants)) + 
  geom_col(aes(fill=type)) + facet_grid(~population) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genomewide Variant Counts") + theme(plot.title = element_text(hjust = 0.5))

pgenome

grid.arrange(pchr1, pgenome, ncol=1)

png("pop_variants.png", width=23, height=13, units="in", res=300)
grid.arrange(pchr1, pgenome, ncol=1)
dev.off()



### todo

cnts

ggplot(cnts, aes(x=genome, y=count)) + 
  geom_col(aes(fill=type)) + facet_grid(~population + chr) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Per Chromosome Variant Counts") + theme(plot.title = element_text(hjust = 0.5))
