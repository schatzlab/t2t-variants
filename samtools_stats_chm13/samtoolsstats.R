## Plot the mapping statistics across samples using the samtools stats matrix
##############################################################################

## Set up environment
# install.packages("gridExtra")
# install.packages("scales")

library(ggplot2)
library(gridExtra)
library(dplyr)
library(scales) 

## this will need to be set to the directory where you have the samtoolstats file
setwd("/Users/mschatz/Dropbox/Documents/Projects/T2T/t2t-variants/samtools_stats_chm13")

## Load tables and fix formatting
stats <- as.data.frame(t(read.table("2021.04.22.samtools.stats.all.txt", header=TRUE, row.names=1)))

stats$Sex                  = as.factor(stats$Sex)
stats$Population_code      = as.factor(stats$Population_code)
stats$Superpopulation_code = as.factor(stats$Superpopulation_code)

stats$error_rate                                = as.double(stats$error_rate)
stats$reads_mapped                              = as.double(stats$reads_mapped)
stats$reads_unmapped                            = as.double(stats$reads_unmapped)
stats$average_length                            = as.double(stats$average_length)
stats$raw_total_sequences                       = as.double(stats$raw_total_sequences)
stats$insert_size_average                       = as.double(stats$insert_size_average)
stats$insert_size_standard_deviation            = as.double(stats$insert_size_standard_deviation)
stats$`percentage_of_properly_paired_reads_(%)` = as.double(stats$`percentage_of_properly_paired_reads_(%)`)


## Population Composition
###############################################################################################################################

## overall population counts

dir.create("plot")

pdf("plot/superpopulation_cnt.pdf", width=11, height=8.5)
ggplot(stats, aes(Superpopulation_code, fill=Superpopulation_code)) + geom_histogram(stat="count")
dev.off()

png("plot/superpopulation_cnt.png", width=11, height=8.5, units="in", res=300)
ggplot(stats, aes(Superpopulation_code, fill=Superpopulation_code)) + geom_histogram(stat="count")
dev.off()



ggplot(stats, aes(Population_code, fill=Superpopulation_code)) + geom_histogram(stat="count") + theme(axis.text.x = element_text(angle = 90))

## overall population counts.... there must be a better way to do this
cnts <- stats %>% count(Population_code, Superpopulation_code)
colors = hue_pal()(5)
pcnt = cnts %>% filter(Superpopulation_code=="AFR"); p1 = ggplot(pcnt, aes(x=Population_code, y=n)) + geom_bar(stat='identity', fill=colors[1]) + xlab("AFR") + ylim(0,max(cnts$n)) + theme(axis.text.x = element_text(angle = 90))
pcnt = cnts %>% filter(Superpopulation_code=="AMR"); p2 = ggplot(pcnt, aes(x=Population_code, y=n)) + geom_bar(stat='identity', fill=colors[2]) + xlab("AMR") + ylim(0,max(cnts$n)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(axis.text.x = element_text(angle = 90))
pcnt = cnts %>% filter(Superpopulation_code=="EAS"); p3 = ggplot(pcnt, aes(x=Population_code, y=n)) + geom_bar(stat='identity', fill=colors[3]) + xlab("EAS") + ylim(0,max(cnts$n)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(axis.text.x = element_text(angle = 90))
pcnt = cnts %>% filter(Superpopulation_code=="EUR"); p4 = ggplot(pcnt, aes(x=Population_code, y=n)) + geom_bar(stat='identity', fill=colors[4]) + xlab("EUR") + ylim(0,max(cnts$n)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(axis.text.x = element_text(angle = 90))
pcnt = cnts %>% filter(Superpopulation_code=="SAS"); p5 = ggplot(pcnt, aes(x=Population_code, y=n)) + geom_bar(stat='identity', fill=colors[5]) + xlab("SAS") + ylim(0,max(cnts$n)) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + theme(axis.text.x = element_text(angle = 90))

pdf("plot/population_cnt.pdf", width=11, height=8.5)
grid.arrange(p1,p2,p3,p4,p5,nrow=1, widths=c(7,4,5,5,5))
dev.off()

png("plot/population_cnt.png", width=11, height=8.5, units="in", res=300)
grid.arrange(p1,p2,p3,p4,p5,nrow=1, widths=c(7,4,5,5,5))
dev.off()


## Coverage Analysis, mapping rate, insert sizes
###############################################################################################################################

## raw reads per superpopulation
pdf("plot/superpopulation_data.pdf", width=11, height=8.5)
ggplot(stats, aes(raw_total_sequences, fill=Superpopulation_code)) + geom_density() + facet_grid(stats$Superpopulation_code)
dev.off()

png("plot/superpopulation_data.png", width=11, height=8.5, units="in", res=300)
ggplot(stats, aes(raw_total_sequences, fill=Superpopulation_code)) + geom_density() + facet_grid(stats$Superpopulation_code)
dev.off()


## raw coverage by superpopulation
pop_cov <- stats %>% group_by(Superpopulation_code) %>% mutate(pop_cov_mean=mean(raw_total_sequences*average_length/3e9)) 

pdf("plot/superpopulation_cov.pdf", width=11, height=8.5)
ggplot(pop_cov, aes(raw_total_sequences*average_length/3e9, fill=Superpopulation_code)) + geom_density() + geom_vline(aes(xintercept=pop_cov_mean), color="black")+ facet_grid(stats$Superpopulation_code)
dev.off()

png("plot/superpopulation_cov.png", width=11, height=8.5, units="in", res=300)
ggplot(pop_cov, aes(raw_total_sequences*average_length/3e9, fill=Superpopulation_code)) + geom_density() + geom_vline(aes(xintercept=pop_cov_mean), color="black")+ facet_grid(stats$Superpopulation_code)
dev.off()



## Mapping Rate
pdf("plot/superpopulation_mapping_rate.pdf", width=11, height=8.5)
ggplot(stats, aes(reads_mapped/raw_total_sequences, fill=Sex)) + geom_density(alpha=0.4) + facet_grid(stats$Superpopulation_code)
dev.off()

png("plot/superpopulation_mapping_rate.png", width=11, height=8.5, units="in", res=300)
ggplot(stats, aes(reads_mapped/raw_total_sequences, fill=Sex)) + geom_density(alpha=0.4) + facet_grid(stats$Superpopulation_code)
dev.off()



## Mapping rates by pop + sex
pop_mapping <- stats %>% group_by(Superpopulation_code, Sex) %>% mutate(pop_mean_mapping=mean(reads_mapped/raw_total_sequences)) 
ggplot(pop_mapping, aes(reads_mapped/raw_total_sequences, fill=Sex)) + geom_density(alpha=0.4) + xlim(.9975, .9995) + geom_vline(aes(xintercept=pop_mean_mapping, color=Sex), size=1.5) + facet_grid(stats$Superpopulation_code) 

## insert size distribution
ggplot(stats, aes(x=insert_size_average, y=insert_size_standard_deviation, color=Superpopulation_code)) + geom_point()
ggplot(stats, aes(x=insert_size_average, fill=Sex)) + geom_density(alpha=0.4) + facet_grid(stats$Superpopulation_code)
ggplot(stats, aes(x=insert_size_standard_deviation, fill=Sex)) + geom_density(alpha=0.4) + facet_grid(stats$Superpopulation_code)


## proper pairs
###############################################################################################################################

## overall proper pairs
p1 <- ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`), fill=D)                    + geom_density(color="black", fill="darkblue") + theme(legend.position='top')
p2 <- ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex))                  + geom_density(alpha=0.4)                      + theme(legend.position='top')
p3 <- ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Superpopulation_code)) + geom_density(alpha=0.4)                      + theme(legend.position='top')
grid.arrange(p1, p2, p3, nrow=3)

## properpairs by population
ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Superpopulation_code)) + geom_density(alpha=0.4) + facet_grid(stats$Superpopulation_code) + theme(legend.position='top') 

## scatterplot of proper pairs by insert size
ggplot(stats, aes(x=insert_size_average, y=`percentage_of_properly_paired_reads_(%)`, color=Superpopulation_code)) + geom_point()

## All in 1 giant messy plot
ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(stats$Population_code)
ggplot(stats, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Superpopulation_code)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(stats$Population_code)

## Per superpopulation
for (pop in levels(stats$Superpopulation_code)) { print(pop) }
pop = "AFR"; superpop = filter(stats, Superpopulation_code==pop); ggplot(superpop, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(superpop$Population_code) + ggtitle(paste("Superpopulation: ", pop)) + xlim(90, 100)
pop = "AMR"; superpop = filter(stats, Superpopulation_code==pop); ggplot(superpop, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(superpop$Population_code) + ggtitle(paste("Superpopulation: ", pop)) + xlim(90, 100)  
pop = "EAS"; superpop = filter(stats, Superpopulation_code==pop); ggplot(superpop, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(superpop$Population_code) + ggtitle(paste("Superpopulation: ", pop)) + xlim(90, 100)  
pop = "EUR"; superpop = filter(stats, Superpopulation_code==pop); ggplot(superpop, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(superpop$Population_code) + ggtitle(paste("Superpopulation: ", pop)) + xlim(90, 100)  
pop = "SAS"; superpop = filter(stats, Superpopulation_code==pop); ggplot(superpop, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex)) + geom_density(alpha=0.4) + theme(legend.position='top') + facet_grid(superpop$Population_code) + ggtitle(paste("Superpopulation: ", pop)) + xlim(90, 100)  

## focus on small insert sizes
small = filter(stats, insert_size_average < 500)
p1 <- ggplot(small, aes(x=`percentage_of_properly_paired_reads_(%)`), fill=D)                    + geom_density(color="black", fill="darkblue") + theme(legend.position='top')
p2 <- ggplot(small, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex))                  + geom_density(alpha=0.4)                      + theme(legend.position='top')
p3 <- ggplot(small, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Superpopulation_code)) + geom_density(alpha=0.4)                      + theme(legend.position='top')
grid.arrange(p1, p2, p3, nrow=3)

## focus on big
big = filter(stats, insert_size_average >= 500)
p1 <- ggplot(big, aes(x=`percentage_of_properly_paired_reads_(%)`), fill=D)                    + geom_density(color="black", fill="darkblue") + theme(legend.position='top')
p2 <- ggplot(big, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Sex))                  + geom_density(alpha=0.4)                      + theme(legend.position='top')
p3 <- ggplot(big, aes(x=`percentage_of_properly_paired_reads_(%)`, fill=Superpopulation_code)) + geom_density(alpha=0.4)                      + theme(legend.position='top')
grid.arrange(p1, p2, p3, nrow=3)


## Error Rate
###############################################################################################################################
 
## error_rate (all)
sex_err <- stats %>% group_by(Sex) %>% mutate(sex_err_mean=mean(error_rate)) 
pop_err <- stats %>% group_by(Superpopulation_code) %>% mutate(pop_err_mean=mean(error_rate)) 
sex_pop_err <- stats %>% group_by(Sex, Superpopulation_code) %>% mutate(sex_pop_err_mean=mean(error_rate)) 

p1 <- ggplot(stats, aes(error_rate, fill=D))                      + geom_density(color="black", fill="darkblue") + geom_vline(aes(xintercept=mean(error_rate)), size=1, color="red") + theme(legend.position='top')
p2 <- ggplot(sex_err, aes(error_rate, fill=Sex))                  + geom_density(alpha=0.4) + geom_vline(aes(xintercept=sex_err_mean, color=Sex), size=1.5) + theme(legend.position='top')
p3 <- ggplot(pop_err, aes(error_rate, fill=Superpopulation_code)) + geom_density(alpha=0.4) + geom_vline(aes(xintercept=pop_err_mean, color=Superpopulation_code), size=1.5) + theme(legend.position='top')

pdf("plot/superpopulation_error_rate.pdf", width=11, height=8.5)
grid.arrange(p1, p2, p3, nrow=3)
dev.off()

png("plot/superpopulation_error_rate.png", width=11, height=8.5, units="in", res=300)
grid.arrange(p1, p2, p3, nrow=3)
dev.off()


pdf("plot/superpopulation_error_rate_sex.pdf", width=11, height=8.5)
ggplot(sex_pop_err, aes(error_rate, fill=Sex)) + geom_density(alpha=0.4) + geom_vline(aes(xintercept=sex_pop_err_mean, color=Sex), size=1.5) + theme(legend.position='top') + facet_grid(stats$Superpopulation_code) 
dev.off()

png("plot/superpopulation_error_rate_sex.png", width=11, height=8.5, units="in", res=300)
ggplot(sex_pop_err, aes(error_rate, fill=Sex)) + geom_density(alpha=0.4) + geom_vline(aes(xintercept=sex_pop_err_mean, color=Sex), size=1.5) + theme(legend.position='top') + facet_grid(stats$Superpopulation_code) 
dev.off()


## error rate small ins
small = filter(stats, insert_size_average < 500)
small_sex_err <- small %>% group_by(Sex) %>% mutate(sex_err_mean=mean(error_rate)) 
small_pop_err <- small %>% group_by(Superpopulation_code) %>% mutate(pop_err_mean=mean(error_rate)) 

p1 <- ggplot(small, aes(error_rate, fill=D))                            + geom_density(color="black", fill="darkblue") + geom_vline(aes(xintercept=mean(error_rate)), size=1, color="red") + theme(legend.position='top')
p2 <- ggplot(small_sex_err, aes(error_rate, fill=Sex))                  + geom_density(alpha=0.4) + geom_vline(aes(xintercept=sex_err_mean, color=Sex), size=1.5) + theme(legend.position='top')
p3 <- ggplot(small_pop_err, aes(error_rate, fill=Superpopulation_code)) + geom_density(alpha=0.4) + geom_vline(aes(xintercept=pop_err_mean, color=Superpopulation_code), size=1.5) + theme(legend.position='top')
grid.arrange(p1, p2, p3, nrow=3)



## reads_unmapped
sex_unmapped <- stats %>% group_by(Sex) %>% mutate(sex_unmapped_mean=mean(reads_unmapped/raw_total_sequences)) 
pop_unmapped <- stats %>% group_by(Superpopulation_code) %>% mutate(pop_unmapped_mean=mean(reads_unmapped/raw_total_sequences)) 

p1 <- ggplot(stats, aes(reads_unmapped/raw_total_sequences, fill=D))                    + geom_density(color="black", fill="darkblue") + geom_vline(aes(xintercept=mean(reads_unmapped/raw_total_sequences)), size=1, color="red") + theme(legend.position='top')
p2 <- ggplot(sex_unmapped, aes(reads_unmapped/raw_total_sequences, fill=Sex))                  + geom_density(alpha=0.4) + geom_vline(aes(xintercept=sex_unmapped_mean, color=Sex), size=1.5) + theme(legend.position='top')
p3 <- ggplot(pop_unmapped, aes(reads_unmapped/raw_total_sequences, fill=Superpopulation_code)) + geom_density(alpha=0.4) + geom_vline(aes(xintercept=pop_unmapped_mean, color=Superpopulation_code), size=1.5) + theme(legend.position='top')

pdf("plot/superpopulation_unmapped_rate.pdf", width=11, height=8.5)
grid.arrange(p1, p2, p3, nrow=3)
dev.off()

png("plot/superpopulation_unmapped_rate.png", width=11, height=8.5)
grid.arrange(p1, p2, p3, nrow=3)
dev.off()

