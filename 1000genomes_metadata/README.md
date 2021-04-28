# 1000 genomes metadata

## sample metadata and download links at ENA

- tab delimited file list for the original 2504 panel at [https://www.ebi.ac.uk/ena/data/view/PRJEB31736](https://www.ebi.ac.uk/ena/data/view/PRJEB31736)


- tab delimited file list for the 698 additional samples at [https://www.ebi.ac.uk/ena/browser/view/PRJEB36890](https://www.ebi.ac.uk/ena/browser/view/PRJEB36890)

## sample metadata from 1000 genomes project

- click on "download the list" from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38

- then rename: "igsr-1000 genomes 30x on grch38.tsv.tsv" -> igsr.metadata.tsv

## index files for download from EBI

```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_698_related_high_coverage.sequence.index
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index


$ less 1000G_2504_high_coverage.sequence.index | grep ftp:// | awk -F '\t' '{print $3,$10,$11}' > 2504.idx
$ cat 1000G_698_related_high_coverage.sequence.index | grep ftp:// | awk -F '\t' '{print $3,$10,$11}' > 698.idx
$ cat 2504.idx 698.idx >> all.idx


$ grep '-' batch2.yaml | awk '{print $2}' | sort > batch2.txt
$ sort all.idx > all.s.idx
$ less batch2.idx | awk '{print $3}' | hashcount
```


## pedigree information

Summary available here: [http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt)

make a simple summary of samples:

```
awk '{if (($3 == "0") && ($4 == "0")){s="singleton"} else if (($3!="0") && ($4!="0")){s="trio"} else {s="one_parent"} print $2,s}' 20130606_g1k_3202_samples_ped_population.txt > trio_status.txt
```

```
$ awk '{print $2}' trio_status.txt | hashcount
one_parent	9	0.28%
singleton	2590	80.86%
trio	604	18.86%
```