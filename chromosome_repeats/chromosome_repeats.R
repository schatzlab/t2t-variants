library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/chromosome_repeats")
dir.create("plot")

plot_order=c("N", "all", "250", "100", "50", "25")
colors = c("darkgrey", brewer_pal(5,"Spectral")(5))

## chm13 uniqueness
###############################################################################
chm = read.table("chm13v1.uniqueness", header=TRUE)

N_row = which(grepl("N", chm$k))
for (row in (N_row-1):2)
{
  for (col in 2:dim(chm)[2])
  {
    chm[row,col] = chm[row,col] - chm[row-1,col]
  }
}

all_row = which(grepl("all", chm$k))
for (col in 2:dim(chm)[2])
{
  chm[all_row,col] = chm[all_row,col] - sum(chm[1:(all_row-1), col])
}

chm$k = factor(chm$k, plot_order)

chm_long = pivot_longer(chm, cols=2:24, names_to="chr", values_to="len")
chm_long$sample = "CHM13v1.0"
chm_long$chr = factor(chm_long$chr, chm_long$chr[1:23])

chm_long

ggplot(chm_long, aes(x=chr, y=len)) + geom_col(aes(fill=k)) +
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("CHM13v1 Chromosome Uniquenes") + theme(plot.title = element_text(hjust = 0.5))


## hg38 uniqueness
###############################################################################

hg38 = read.table("hg38.uniqueness", header=TRUE)

N_row = which(grepl("N", hg38$k))
for (row in (N_row-1):2)
{
  for (col in 2:dim(hg38)[2])
  {
    hg38[row,col] = hg38[row,col] - hg38[row-1,col]
  }
}

all_row = which(grepl("all", hg38$k))
for (col in 2:dim(hg38)[2])
{
  hg38[all_row,col] = hg38[all_row,col] - sum(hg38[1:(all_row-1), col])
}

hg38$k = factor(hg38$k, plot_order)

hg38_long = pivot_longer(hg38, cols=2:24, names_to="chr", values_to="len")
hg38_long$sample = "GRCh38"
hg38_long$chr = factor(hg38_long$chr, chm_long$chr[1:23])

hg38_long

#ggplot(hg38_long, aes(x=chr, y=len)) + geom_col(aes(fill=k)) +
#  scale_fill_brewer(palette="Spectral") +
#  theme(axis.text.x = element_text(angle = 90)) + ggtitle("GRCh38 Chromosome Uniquenes") + theme(plot.title = element_text(hjust = 0.5))

ggplot(hg38_long, aes(x=chr, y=len)) + geom_col(aes(fill=k)) +
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(angle = 90)) + ggtitle("GRCh38 Chromosome Uniquenes") + theme(plot.title = element_text(hjust = 0.5))


## Compare CHM13 with GRCh38
###############################################################################

all_long <- rbind(hg38_long, chm_long)

all_long

chr_plot = ggplot(all_long, aes(x=sample, y=len)) + 
  geom_col(aes(fill=k)) + facet_grid(~chr) +
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Chromosome") + theme(plot.title = element_text(hjust = 0.5))
  
chr_plot  

chm_total  = chm  %>% rowwise(k) %>% mutate(total=sum(c_across(2:23)))
hg38_total = hg38 %>% rowwise(k) %>% mutate(total=sum(c_across(2:23)))

chm_total
hg38_total

total=data.frame(k=chm_total$k, "CHM13v1.0" = chm_total$total, "GRCh38" = hg38_total$total)

total

total_long = pivot_longer(total, cols=2:3, names_to="sample", values_to="len")

total_long

genome_plot = ggplot(total_long, aes(x=sample, y=len)) + 
  geom_col(aes(fill=k)) +
  scale_fill_manual(values=colors) +
  theme(axis.text.x = element_text(hjust=1, angle=90)) + theme(axis.title.x = element_blank()) +
  ggtitle("Genome") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none")

genome_plot

grid.arrange(genome_plot, chr_plot, nrow=1, widths=c(4,27))







total_change=data.frame(k=chm_total$k, "Gains"=chm_total$total - hg38_total$total)

total_change

ggplot(total_change, aes(x=k, y=Gains)) + 
  geom_col(aes(fill=k)) +
  scale_fill_manual(values=colors) +
  ggtitle("Uniqueness Gains") + theme(plot.title = element_text(hjust = 0.5))




## old stuff
###############################################################################



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

png("plot/chm13_chromosome_uniqueness.png", width=11, height=8.5, units="in", res=300)
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

hg38_long

png("plot/grch38_chromosome_uniqueness.png", width=11, height=8.5, units="in", res=300)
ggplot(hg38_long, aes(x=chr, y=len)) + geom_col(aes(fill=type)) + ggtitle("GRCh38 Chromosome Uniquenes (75-mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


## genome comparison
###############################################################################

all_long <- rbind(hg38_long, chm_long)

pdf("plot/cmp_chromosome_uniqueness.pdf", width=11, height=8.5)
ggplot(all_long, aes(x=sample, y=len)) + geom_col(aes(fill=type)) + facet_grid(~chr) + theme(axis.text.x = element_text(angle = 90)) + ggtitle("Chromosome Uniquenes (75-mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("plot/cmp_chromosome_uniqueness.png", width=11, height=8.5, units="in", res=300)
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

png("plot/cmp_chromosome_extra_unique.png", width=11, height=8.5, units="in", res=300)
ggplot(extra_unique, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Additional Unique Bases per Chromosome (75 mers)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


## extra total bases per chromosome
###############################################################################

extra_total = data.frame(chr=hg38$chr, len=chm$total-hg38$total)
sum(extra_total$len)

pdf("plot/cmp_chromosome_extra_bases.pdf", width=11, height=8.5)
ggplot(extra_total, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (all bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("plot/cmp_chromosome_extra_bases.png", width=11, height=8.5, units="in", res=300)
ggplot(extra_total, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (all bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


## extra non-n
###############################################################################

extra_nonn = data.frame(chr=hg38$chr, len=(chm$total-chm$n)-(hg38$total-hg38$n))
sum(extra_nonn$len)

pdf("plot/cmp_chromosome_extra_nonn.pdf", width=11, height=8.5)
ggplot(extra_nonn, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (non-n bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

png("plot/cmp_chromosome_extra_nonn.png", width=11, height=8.5, units="in", res=300)
ggplot(extra_nonn, aes(x=chr, y=len)) + geom_col(fill=colors[3]) + labs(title="Chromosome length differences (non-n bases)") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

