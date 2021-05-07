library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/chromosome_variants")


## chm13 on its own
###############################################################################

chm <- read.table("chm13.chr.stats", header=TRUE)

chm$chr = factor(chm$chr, chm$chr)
chm_long = pivot_longer(chm, cols=2:3, names_to="type", values_to="count")
chm_long$sample = "CHM13"

chm_long

ggplot(chm_long, aes(x=chr, y=count)) + geom_col(aes(fill=type)) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("CHM13 Variants") + 
  theme(plot.title = element_text(hjust = 0.5))



## hg38 on its own
###############################################################################

hg38 <- read.table("hg38.chr.stats", header=TRUE)

hg38$chr = factor(hg38$chr, hg38$chr)
hg38_long = pivot_longer(hg38, cols=2:3, names_to="type", values_to="count")
hg38_long$sample = "GRCh38"

hg38_long

ggplot(hg38_long, aes(x=chr, y=count)) + geom_col(aes(fill=type)) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("GRCh38 Variants") + 
  theme(plot.title = element_text(hjust = 0.5))

## combined
###############################################################################

all_long <- rbind(hg38_long, chm_long)

all_long

ggplot(all_long, aes(x=sample, y=count)) + 
  geom_col(aes(fill=type)) + facet_grid(~chr) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Variant Counts") + theme(plot.title = element_text(hjust = 0.5))


png("chr_variants.png", width=23, height=13, units="in", res=300)
ggplot(all_long, aes(x=sample, y=count)) + 
  geom_col(aes(fill=type)) + facet_grid(~chr) +
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Variant Counts") + theme(plot.title = element_text(hjust = 0.5))
dev.off()
