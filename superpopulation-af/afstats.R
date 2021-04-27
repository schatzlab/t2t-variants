library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 
library(tidyr)
library(scales)

setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/superpopulation-af")


## AFR
###############################################################################################

chm_afr <- read.table("chm13_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_afr <- read.table("grch38_AFR_autosomes_summarized.txt.af", header=TRUE, sep='\t')


plot(chm_afr$X.3.allele.frequency, chm_afr$X.4.number.of.SNPs, log="y")
plot(hg38_afr$X.3.allele.frequency, hg38_afr$X.4.number.of.SNPs, log="y")

plot(chm_afr$X.3.allele.frequency, chm_afr$X.4.number.of.SNPs, log="y")
points(hg38_afr$X.3.allele.frequency, hg38_afr$X.4.number.of.SNPs, log="y", col="red")

plot(hg38_afr$X.3.allele.frequency, chm_afr$X.3.allele.frequency)
plot(hg38_afr$X.4.number.of.SNPs, chm_afr$X.4.number.of.SNPs, log="xy")


afr<-data.frame(af=chm_afr$X.3.allele.frequency, chm13=chm_afr$X.4.number.of.SNPs, hg38=hg38_afr$X.4.number.of.SNPs)
afr_long <- pivot_longer(afr, cols=2:3, names_to="reference", values_to="snps")
head(afr_long)

afrp <- ggplot(afr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), size=1.5) + coord_trans(y="log2") + 
  ggtitle("AFR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")



## AMR
###############################################################################################

chm_amr <- read.table("chm13_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_amr <- read.table("grch38_AMR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

amr<-data.frame(af=chm_amr$X.3.allele.frequency, chm13=chm_amr$X.4.number.of.SNPs, hg38=hg38_amr$X.4.number.of.SNPs)
amr_long <- pivot_longer(amr, cols=2:3, names_to="reference", values_to="snps")
head(amr_long)

amrp<-ggplot(amr_long) + 
  geom_line(aes(x=af, y=snps, color=reference), size=1.5) + coord_trans(y="log2") + 
  ggtitle("AMR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## EAS
###############################################################################################

chm_eas <- read.table("chm13_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eas <- read.table("grch38_EAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eas<-data.frame(af=chm_eas$X.3.allele.frequency, chm13=chm_eas$X.4.number.of.SNPs, hg38=hg38_eas$X.4.number.of.SNPs)
eas_long <- pivot_longer(eas, cols=2:3, names_to="reference", values_to="snps")
head(eas_long)

easp<-ggplot(eas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), size=1.5) + coord_trans(y="log2") + 
  ggtitle("EAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## EUR
###############################################################################################

chm_eur <- read.table("chm13_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_eur <- read.table("grch38_EUR_autosomes_summarized.txt.af", header=TRUE, sep='\t')

eur<-data.frame(af=chm_eur$X.3.allele.frequency, chm13=chm_eur$X.4.number.of.SNPs, hg38=hg38_eur$X.4.number.of.SNPs)
eur_long <- pivot_longer(eur, cols=2:3, names_to="reference", values_to="snps")
head(eur_long)

eurp<-ggplot(eur_long) + 
  geom_line(aes(x=af, y=snps, color=reference), size=1.5) + coord_trans(y="log2") + 
  ggtitle("EUR Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


## SAS
###############################################################################################

chm_sas <- read.table("chm13_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')
hg38_sas <- read.table("grch38_SAS_autosomes_summarized.txt.af", header=TRUE, sep='\t')

sas<-data.frame(af=chm_sas$X.3.allele.frequency, chm13=chm_sas$X.4.number.of.SNPs, hg38=hg38_sas$X.4.number.of.SNPs)
sas_long <- pivot_longer(sas, cols=2:3, names_to="reference", values_to="snps")
head(sas_long)

sasp<-ggplot(sas_long) + 
  geom_line(aes(x=af, y=snps, color=reference), size=1.5) + coord_trans(y="log2") + 
  ggtitle("SAS Allele Frequency Distribution") + xlab("allele frequency") + ylab("Number of Variants")


grid.arrange(afrp,amrp,easp,eurp,sasp,ncol=1)
