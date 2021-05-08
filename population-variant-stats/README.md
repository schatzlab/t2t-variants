### CHM13 preprocessing

Download: chromosome table
https://anvil.terra.bio/#workspaces/t2t-nhgri/t2t-vc-processing/data

```
$ tr '\t' '\n' < chromosome.tsv  | grep singleton.bcftools.stats.txt > singleton.stats.urls
$ cat singleton.stats.urls | gsutil -m cp -I .
$ mkdir global population superpopulation
$ mv *AFR* *AMR* *EAS* *EUR* *SAS* superpopulation/
$ mv *.singleton.bcftools.* global/
$ mv *.txt population/

$ cd population
$ for pop in `ls *chr1.* | cut -f7 -d '.' | cut -f1 -d'_'`; do (echo -e "chr\tsnps\tindels"; for file in `ls *$pop*.txt | sort -V `; do chr=`echo $file | cut -f2 -d'.'`; grep '^AF' $file | awk -v chr=$chr '{if ($3 > 0){s += $4; i += $7}} END{print chr"\t"s"\t"i}'; done) > $pop.chr.stats; done


## separate out hi frequency variants
$ for pop in `ls *chr1.* | cut -f7 -d '.' | cut -f1 -d'_'`; do (echo -e "chr\tsnps\tindels\tsnps_hi\tindels_hi"; for file in `ls *$pop*.txt | sort -V `; do chr=`echo $file | cut -f2 -d'.'`; grep '^AF' $file | awk -v chr=$chr '{if ($3 >= 1){shi += $4; ihi += $7} else if ($3 > 0){s += $4; i += $7}} END{print chr"\t"s"\t"i"\t"shi"\t"ihi}'; done) > $pop.chr.stats; done
```




### GRCh38 preprocessing:

Download: grch38_chromosome table
https://anvil.terra.bio/#workspaces/t2t-nhgri/t2t-vc-processing/data

```
$ tr '\t' '\n' < grch38_chromosome.tsv  | grep singleton.bcftools.stats.txt > singleton.stats.urls
$ cat singleton.stats.urls | gsutil -m cp -I .
$ mkdir global population superpopulation
$ mv *AFR* *AMR* *EAS* *EUR* *SAS* superpopulation/
$ mv *.singleton.bcftools.* global/
$ mv *.txt population/
$ cd population
$ for pop in `ls *chr1.* | cut -f5 -d'.' | cut -f1 -d'_'`; do (echo -e "chr\tsnps\tindels"; for file in `ls *$pop*.txt | sort -V `; do chr=`echo $file | cut -f1 -d'.' | rev | cut -f1 -d'_' | rev`; grep '^AF' $file | awk -v chr=$chr '{if ($3 > 0){s += $4; i += $7}} END{print chr"\t"s"\t"i}'; done) > $pop.chr.stats; done
```


