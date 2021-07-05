library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)
library(tibble)
library(ggallin)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/novel_region_af")

alpha_level = 0.6
size_level  = 1.5


## Non-syntenic regions
###############################################################################

nonsyn_afr <- read.table("non_syntenic_AFR_autosomes.txt.af", header=TRUE, sep='\t')
nonsyn_amr <- read.table("non_syntenic_AMR_autosomes.txt.af", header=TRUE, sep='\t')
nonsyn_eas <- read.table("non_syntenic_EAS_autosomes.txt.af", header=TRUE, sep='\t')
nonsyn_eur <- read.table("non_syntenic_EUR_autosomes.txt.af", header=TRUE, sep='\t')
nonsyn_sas <- read.table("non_syntenic_SAS_autosomes.txt.af", header=TRUE, sep='\t')

nonsyn_all <- rbind(nonsyn_afr %>% add_column(pop="AFR"),
                    nonsyn_amr %>% add_column(pop="AMR"),
                    nonsyn_eas %>% add_column(pop="EAS"),
                    nonsyn_eur %>% add_column(pop="EUR"),
                    nonsyn_sas %>% add_column(pop="SAS"))

head(nonsyn_all)

nonsyn_allp <- ggplot(nonsyn_all) + 
  geom_line(aes(x=X.3.allele.frequency, y=X.4.number.of.SNPs, color=pop), alpha=alpha_level, size=size_level) +
  facet_grid(~pop) +
  scale_y_log10() +
  ggtitle("Non-syntenic Regions - Allele Frequency") + xlab("allele frequency") + ylab("Number of Variants") +
  theme(legend.position = "none")

nonsyn_allp


## Novel regions
###############################################################################

novel_afr <- read.table("novel_AFR_autosomes.txt.af", header=TRUE, sep='\t')
novel_amr <- read.table("novel_AMR_autosomes.txt.af", header=TRUE, sep='\t')
novel_eas <- read.table("novel_EAS_autosomes.txt.af", header=TRUE, sep='\t')
novel_eur <- read.table("novel_EUR_autosomes.txt.af", header=TRUE, sep='\t')
novel_sas <- read.table("novel_SAS_autosomes.txt.af", header=TRUE, sep='\t')

novel_all <- rbind(novel_afr %>% add_column(pop="AFR"),
                    novel_amr %>% add_column(pop="AMR"),
                    novel_eas %>% add_column(pop="EAS"),
                    novel_eur %>% add_column(pop="EUR"),
                    novel_sas %>% add_column(pop="SAS"))

head(novel_all)

novel_allp <- ggplot(novel_all) + 
  geom_line(aes(x=X.3.allele.frequency, y=X.4.number.of.SNPs, color=pop), alpha=alpha_level, size=size_level) +
  facet_grid(~pop) +
  scale_y_log10() +
  ggtitle("Novel Regions - Allele Frequency") + xlab("allele frequency") + ylab("Number of Variants") +
  theme(legend.position = "none")

novel_allp


grid.arrange(nonsyn_allp, novel_allp, ncol=1)

