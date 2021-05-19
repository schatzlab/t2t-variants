library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/superpopulation-af")

alpha_level = 0.6
size_level  = 1.5


## AFR
###############################################################################################

chm_afr <- read.table("singleton_chm13_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_afr <- read.table("singleton_grch38_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

afr <- data.frame(af=chm_afr$X.3.allele.frequency, chm13=chm_afr$X.4.number.of.SNPs, hg38=hg38_afr$X.4.number.of.SNPs)
afr_long <- pivot_longer(afr, cols=2:3, names_to="reference", values_to="snps")

## filter af == 0
afr_long = afr_long %>% filter(af > 0)

afrp <- ggplot(afr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=snps, color=reference)) +
  scale_y_log10() +
  ggtitle("AFR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants") + labs(tag="A")

afrp

## Cummulate number of variants
afr_cum = afr_long %>% group_by(reference) %>% mutate(cum_snps = cumsum(snps))

afrp_cum <- ggplot(afr_cum) + 
  geom_line(aes(x=af, y=cum_snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=cum_snps, color=reference)) +
  ggtitle("AFR Allele Cumulative Frequency Distribution") + xlab("max allele frequency") + ylab("Number of Variants") + labs(tag="A")

afrp_cum



## AMR
###############################################################################################

chm_amr <- read.table("singleton_chm13_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_amr <- read.table("singleton_grch38_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

amr<-data.frame(af=chm_amr$X.3.allele.frequency, chm13=chm_amr$X.4.number.of.SNPs, hg38=hg38_amr$X.4.number.of.SNPs)
amr_long <- pivot_longer(amr, cols=2:3, names_to="reference", values_to="snps")
amr_long = amr_long %>% filter(af > 0)

amrp<-ggplot(amr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=snps, color=reference)) +
  scale_y_log10() +
  ggtitle("AMR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants") + labs(tag="B")

amrp

## Cummulate number of variants
amr_cum = amr_long %>% group_by(reference) %>% mutate(cum_snps = cumsum(snps))

amrp_cum <- ggplot(amr_cum) + 
  geom_line(aes(x=af, y=cum_snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=cum_snps, color=reference)) +
  ggtitle("AMR Allele Cumulative Frequency Distribution") + xlab("max allele frequency") + ylab("Number of Variants") + labs(tag="A")

amrp_cum


## EAS
###############################################################################################

chm_eas <- read.table("singleton_chm13_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eas <- read.table("singleton_grch38_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eas<-data.frame(af=chm_eas$X.3.allele.frequency, chm13=chm_eas$X.4.number.of.SNPs, hg38=hg38_eas$X.4.number.of.SNPs)
eas_long <- pivot_longer(eas, cols=2:3, names_to="reference", values_to="snps")
eas_long = eas_long %>% filter(af > 0)

easp<-ggplot(eas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=snps, color=reference)) +
  scale_y_log10() +
  ggtitle("EAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants") + labs(tag="C")

easp

## Cummulate number of variants
eas_cum = eas_long %>% group_by(reference) %>% mutate(cum_snps = cumsum(snps))

easp_cum <- ggplot(eas_cum) + 
  geom_line(aes(x=af, y=cum_snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=cum_snps, color=reference)) +
  ggtitle("EAS Allele Cumulative Frequency Distribution") + xlab("max allele frequency") + ylab("Number of Variants") + labs(tag="A")

easp_cum


## EUR
###############################################################################################

chm_eur <- read.table("singleton_chm13_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eur <- read.table("singleton_grch38_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eur<-data.frame(af=chm_eur$X.3.allele.frequency, chm13=chm_eur$X.4.number.of.SNPs, hg38=hg38_eur$X.4.number.of.SNPs)
eur_long <- pivot_longer(eur, cols=2:3, names_to="reference", values_to="snps")
eur_long = eur_long %>% filter(af > 0)

eurp<-ggplot(eur_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=snps, color=reference)) +
  scale_y_log10() +
  ggtitle("EUR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants") + labs(tag="D")

eurp

## Cummulate number of variants
eur_cum = eur_long %>% group_by(reference) %>% mutate(cum_snps = cumsum(snps))

eurp_cum <- ggplot(eur_cum) + 
  geom_line(aes(x=af, y=cum_snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=cum_snps, color=reference)) +
  ggtitle("EUR Allele Cumulative Frequency Distribution") + xlab("max allele frequency") + ylab("Number of Variants") + labs(tag="A")

eurp_cum


## SAS
###############################################################################################

chm_sas <- read.table("singleton_chm13_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_sas <- read.table("singleton_grch38_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

sas<-data.frame(af=chm_sas$X.3.allele.frequency, chm13=chm_sas$X.4.number.of.SNPs, hg38=hg38_sas$X.4.number.of.SNPs)
sas_long <- pivot_longer(sas, cols=2:3, names_to="reference", values_to="snps")
sas_long = sas_long %>% filter(af > 0)

sasp<-ggplot(sas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=snps, color=reference)) +
  scale_y_log10() +
  ggtitle("SAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants") + labs(tag="E")

sasp

## Cummulate number of variants
sas_cum = sas_long %>% group_by(reference) %>% mutate(cum_snps = cumsum(snps))

sasp_cum <- ggplot(sas_cum) + 
  geom_line(aes(x=af, y=cum_snps, color=reference), alpha=alpha_level, size=size_level) +
  geom_point(aes(x=af, y=cum_snps, color=reference)) +
  ggtitle("SAS Allele Cumulative Frequency Distribution") + xlab("max allele frequency") + ylab("Number of Variants") + labs(tag="A")

sasp_cum


## Plot all in a grid
###############################################################################################

grid.arrange(afrp,amrp,easp,eurp,sasp,ncol=1)


## plot the difference in the number of variants
###############################################################################################

afrd=data.frame(af=chm_afr$X.3.allele.frequency, diff_snps=chm_afr$X.4.number.of.SNPs-hg38_afr$X.4.number.of.SNPs, pop="AFR")
amrd=data.frame(af=chm_amr$X.3.allele.frequency, diff_snps=chm_amr$X.4.number.of.SNPs-hg38_amr$X.4.number.of.SNPs, pop="AMR")
easd=data.frame(af=chm_eas$X.3.allele.frequency, diff_snps=chm_eas$X.4.number.of.SNPs-hg38_eas$X.4.number.of.SNPs, pop="EAS")
eurd=data.frame(af=chm_eur$X.3.allele.frequency, diff_snps=chm_eur$X.4.number.of.SNPs-hg38_eur$X.4.number.of.SNPs, pop="EUR")
sasd=data.frame(af=chm_sas$X.3.allele.frequency, diff_snps=chm_sas$X.4.number.of.SNPs-hg38_sas$X.4.number.of.SNPs, pop="SAS")

all_diff = bind_rows(afrd, amrd, easd, eurd, sasd)

diff_all = ggplot(all_diff) +
  geom_line(aes(x=af, y=diff_snps, color=pop), size=size_level, alpha=alpha_level) +
  ggtitle("Allele Frequency Distribution Difference (CHM13 - GRCh38)") + 
  xlab("allele frequency") + ylab("Difference in Number of Variants") + labs(tag="F")

diff_zoom = ggplot(all_diff) +
  geom_line(aes(x=af, y=diff_snps, color=pop), size=size_level, alpha=alpha_level) + coord_cartesian(ylim=c(-10000,10000)) +
  ggtitle("Zoomed Allele Frequency Distribution Difference (CHM13 - GRCh38)") + 
  xlab("allele frequency") + ylab("Difference in Number of Variants") + labs(tag="G")

grid.arrange(diff_all, diff_zoom, ncol=1)


## extracted "fixed" snps
fixed = all_diff %>% filter(af > 0.4) %>% filter(af < .6)

ggplot(fixed) +
  geom_line(aes(x=af, y=diff_snps, color=pop), size=size_level, alpha=alpha_level) + coord_cartesian(ylim=c(-10000,10000)) +
  ggtitle("Zoomed Allele Frequency Distribution Difference (CHM13 - GRCh38)") + xlab("allele frequency") + ylab("Difference in Number of Variants")

fixed_cnt = all_diff %>% filter(af > 0.495) %>% filter(af < .505)
sum(fixed_cnt$diff_snps)


## make a nice figure that combines it all together
###############################################################################################

lay=rbind(c(1,6),
          c(1,6),
          c(2,6),
          c(2,6),
          c(3,6),
          c(3,7),
          c(4,7),
          c(4,7),
          c(5,7),
          c(5,7))

grid.arrange(afrp,amrp,easp,eurp,sasp, diff_all, diff_zoom, layout_matrix=lay)

png("allele_frequency.png", width=23, height=13, units="in", res=300)
grid.arrange(afrp,amrp,easp,eurp,sasp, diff_all, diff_zoom, layout_matrix=lay)
dev.off()

pdf("allele_frequency.pdf", width=23, height=13)
grid.arrange(afrp,amrp,easp,eurp,sasp, diff_all, diff_zoom, layout_matrix=lay)
dev.off()



## Plot on the same plot (dont do this - doenst show differences very well)
###############################################################################################

afr_df = afr_long %>% add_column(pop="AFR")
amr_df = amr_long %>% add_column(pop="AMR")
eas_df = eas_long %>% add_column(pop="EAS")
eur_df = eur_long %>% add_column(pop="EUR")
sas_df = sas_long %>% add_column(pop="SAS")

all_df=bind_rows(afr_df, amr_df, eas_df, eur_df, sas_df)

all_chm  = all_df %>% filter(reference=="chm13")
all_hg38 = all_df %>% filter(reference=="hg38")

all_hg38p = ggplot(all_hg38) + 
  geom_line(aes(x=af, y=snps, color=pop), size=size_level, alpha=alpha_level) + coord_trans(y="log2") + 
  ggtitle("GRCh38 Allele Frequency Distributions") + xlab("allele frequency") + ylab("Number of Variants")

all_chmp = ggplot(all_chm) + 
  geom_line(aes(x=af, y=snps, color=pop), size=size_level, alpha=alpha_level) + coord_trans(y="log2") + 
  ggtitle("CHM13 Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")

grid.arrange(all_hg38p, all_chmp,ncol=1)



## overall collection - dont do this!
###############################################################################################


chm_all  <- read.table("chm13_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_all <- read.table("grch38_autosomes_summarized.txt.af", header=TRUE, sep='\t')

chm_all_df  = data.frame(af=chm_all$X.3.allele.frequency, num_snps=chm_all$X.4.number.of.SNPs, reference="CHM13")
hg38_all_df = data.frame(af=hg38_all$X.3.allele.frequency, num_snps=hg38_all$X.4.number.of.SNPs, reference="GRCh38")

all_df = bind_rows(hg38_all_df, chm_all_df)

ggplot(all_df) +
  geom_point(aes(x=af, y=num_snps, color=reference)) + coord_trans(y="log2") +
  ggtitle("Allele Frequency Distributions") + ylab("Difference in Number of Variants")




