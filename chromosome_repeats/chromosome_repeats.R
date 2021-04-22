library(dplyr)   # For data manipulation
library(ggplot2) # For data visualization
library(tidyr)
library(scales)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/chromosome_repeats")
dir.create("plot")

## chm13 
###############################################################################

chm_uniq  <- read.table("chm13v1_75mer.txt")
chm_n     <- read.table("chm13v1_N.txt")

chm <- data.frame(chr=chm_uniq[[1]], unique=chm_uniq[[2]], n=chm_n[[2]], repeats=chm_uniq[[4]]-chm_uniq[[2]]-chm_n[[2]], total=chm_uniq[[4]])
chm$chr = factor(chm$chr, chm$chr)
chm_long <- pivot_longer(chm, cols=2:4, names_to="type", values_to="len")
chm_long$sample <- "CHM13"

pdf("plot/chm13_chromosome_uniqueness.pdf", width=11, height=8.5)
ggplot(chm_long, aes(x=chr, y=len)) + geom_col(aes(fill=type)) + ggtitle("CHM13 Chromosome Uniquenes (75-mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## hg38
###############################################################################

hg38_uniq <- read.table("GRCh38p13_75mer.txt")
hg38_n    <- read.table("GRCh38p13_N.txt")

hg38 <- data.frame(chr=hg38_uniq[[1]], unique=hg38_uniq[[2]], n=hg38_n[[2]], repeats=hg38_uniq[[4]]-hg38_uniq[[2]]-hg38_n[[2]], total=hg38_uniq[[4]])
hg38$chr = factor(hg38$chr, hg38$chr)
hg38_long <- pivot_longer(hg38, cols=2:4, names_to="type", values_to="len")
hg38_long$sample <- "GRCh38"

pdf("plot/grch38_chromosome_uniqueness.pdf", width=11, height=8.5)
ggplot(hg38_long, aes(x=chr, y=len)) + geom_col(aes(fill=type)) + ggtitle("GRCh38 Chromosome Uniquenes (75-mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## genome comparison
###############################################################################

all_long <- rbind(hg38_long, chm_long)

pdf("plot/cmp_chromosome_uniqueness.pdf", width=11, height=8.5)
ggplot(all_long, aes(x=sample, y=len)) + geom_col(aes(fill=type)) + facet_grid(~chr) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Chromosome Uniquenes (75-mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()



## extra unique bases per chromosome
###############################################################################

colors=hue_pal()(3)

extra_unique = data.frame(chr=hg38$chr, len=chm$unique-hg38$unique)
sum(extra_unique$len)

pdf("plot/cmp_chromosome_extra_unique.pdf", width=11, height=8.5)
ggplot(extra_unique, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Additional Unique Bases per Chromosome (75 mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## extra total bases per chromosome
###############################################################################

extra_total = data.frame(chr=hg38$chr, len=chm$total-hg38$total)
sum(extra_total$len)

pdf("plot/cmp_chromosome_extra_bases.pdf", width=11, height=8.5)
ggplot(extra_total, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (all bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

## extra non-n
###############################################################################

extra_nonn = data.frame(chr=hg38$chr, len=(chm$total-chm$n)-(hg38$total-hg38$n))
sum(extra_nonn$len)

pdf("plot/cmp_chromosome_extra_nonn.pdf", width=11, height=8.5)
ggplot(extra_nonn, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (non-n bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
