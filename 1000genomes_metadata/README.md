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

The pedigree file includes 9 parent sample that were not actually sequenced: 

```
$ skip.pl 1 20130606_g1k_3202_samples_ped_population.txt | awk '{print $2}' | sort > samples.id
$ skip.pl 1 20130606_g1k_3202_samples_ped_population.txt | awk '{print $3;print $4}' | sort | uniq > parents.id
$ comm -2 -3 parents.id samples.id | awk '{if ($1 != 0){print}}' > parents.exclude

$ cat parents.exclude
HG00144
HG00866
HG02567
HG03466
HG04067
NA19105
NA19432
NA19453
NA20909
```

### Fixing the pedigree file

```
$ awk '{printf " -e s/"$1"/0/"}' parents.exclude
$ sed  -e s/HG00144/0/ -e s/HG00866/0/ -e s/HG02567/0/ -e s/HG03466/0/ -e s/HG04067/0/ -e s/NA19105/0/ -e s/NA19432/0/ -e s/NA19453/0/ -e s/NA20909/0/ 20130606_g1k_3202_samples_ped_population.txt > 20130606_g1k_3202_samples_ped_population.txt.corrected
```


make a simple summary of samples:

```
$ awk '{if (($3 == "0") && ($4 == "0")){s="singleton"} else if (($3!="0") && ($4!="0")){s="trio"} else {s="one_parent"} print $2,s}' 20130606_g1k_3202_samples_ped_population.txt.corrected > trio_status.txt

$ skip.pl 1 trio_status.txt | awk '{print $2}' | hashcount | column -t
one_parent  2     0.06%
singleton   2598  81.14%
trio        602   18.80%
```