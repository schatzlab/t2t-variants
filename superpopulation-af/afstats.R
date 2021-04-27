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

chm_afr <- read.table("chm13_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_afr <- read.table("grch38_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

afr <- data.frame(af=chm_afr$X.3.allele.frequency, chm13=chm_afr$X.4.number.of.SNPs, hg38=hg38_afr$X.4.number.of.SNPs)
afr_long <- pivot_longer(afr, cols=2:3, names_to="reference", values_to="snps")

afrp <- ggplot(afr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) + coord_trans(y="log2") + 
  geom_point(aes(x=af, y=snps, color=reference)) +
  ggtitle("AFR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## AMR
###############################################################################################

chm_amr <- read.table("chm13_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_amr <- read.table("grch38_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

amr<-data.frame(af=chm_amr$X.3.allele.frequency, chm13=chm_amr$X.4.number.of.SNPs, hg38=hg38_amr$X.4.number.of.SNPs)
amr_long <- pivot_longer(amr, cols=2:3, names_to="reference", values_to="snps")

amrp<-ggplot(amr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) + coord_trans(y="log2") + 
  geom_point(aes(x=af, y=snps, color=reference)) +
  ggtitle("AMR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## EAS
###############################################################################################

chm_eas <- read.table("chm13_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eas <- read.table("grch38_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eas<-data.frame(af=chm_eas$X.3.allele.frequency, chm13=chm_eas$X.4.number.of.SNPs, hg38=hg38_eas$X.4.number.of.SNPs)
eas_long <- pivot_longer(eas, cols=2:3, names_to="reference", values_to="snps")

easp<-ggplot(eas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) + coord_trans(y="log2") + 
  geom_point(aes(x=af, y=snps, color=reference)) +
  ggtitle("EAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## EUR
###############################################################################################

chm_eur <- read.table("chm13_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eur <- read.table("grch38_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eur<-data.frame(af=chm_eur$X.3.allele.frequency, chm13=chm_eur$X.4.number.of.SNPs, hg38=hg38_eur$X.4.number.of.SNPs)
eur_long <- pivot_longer(eur, cols=2:3, names_to="reference", values_to="snps")

eurp<-ggplot(eur_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) + coord_trans(y="log2") + 
  geom_point(aes(x=af, y=snps, color=reference)) +
  ggtitle("EUR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## SAS
###############################################################################################

chm_sas <- read.table("chm13_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_sas <- read.table("grch38_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

sas<-data.frame(af=chm_sas$X.3.allele.frequency, chm13=chm_sas$X.4.number.of.SNPs, hg38=hg38_sas$X.4.number.of.SNPs)
sas_long <- pivot_longer(sas, cols=2:3, names_to="reference", values_to="snps")

sasp<-ggplot(sas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), alpha=alpha_level, size=size_level) + coord_trans(y="log2") + 
  geom_point(aes(x=af, y=snps, color=reference)) +
  ggtitle("SAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")

## Plot all in a grid
###############################################################################################

grid.arrange(afrp,amrp,easp,eurp,sasp,ncol=1)



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
  ggtitle("Allele Frequency Distribution Difference (CHM13 - GRCh38)") + xlab("allele frequency") + ylab("Difference in Number of Variants")

diff_zoom = ggplot(all_diff) +
  geom_line(aes(x=af, y=diff_snps, color=pop), size=size_level, alpha=alpha_level) + coord_cartesian(ylim=c(-10000,10000)) +
  ggtitle("Zoomed Allele Frequency Distribution Difference (CHM13 - GRCh38)") + xlab("allele frequency") + ylab("Difference in Number of Variants")

grid.arrange(diff_all, diff_zoom, ncol=1)


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




