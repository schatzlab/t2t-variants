### CHM13 preprocessing

Download: chromosome table
https://anvil.terra.bio/#workspaces/t2t-nhgri/t2t-vc-processing/data

```
$ tr '\t' '\n' < ../chromosome.tsv | grep per_sample > persample.urls
$ cat persample.urls | gsutil -m cp -I .
$ mkdir population superpopulation
$ mv *AFR* *AMR* *EAS* *EUR* *SAS* superpopulation/
$ mv *.stats.txt population/

$ cd population

$ for i in `/bin/ls *.stats.txt`; do chr=`echo $i | cut -f2 -d'.'`; grep -B1 '^PSC' $i | tr ' ' '_' | awk -v chr=$chr '{print chr,$0}' > $i.psc; done

$ head -1 *chr1.*ACB*.psc | sed 's/chr1/chr/' | sed 's/#_PSC/PSC/' | tr ' ' '\t' > header

$ for pop in `ls *chr1.* | grep psc | cut -f 7 -d '.' | cut -f 1 -d'_'`; do echo $pop; (cat header; grep -h -v '#_PSC' *$pop*.stats.txt.psc) | tr ' ' '\t' > $pop.pop.psc; done

```



### GRCh38 preprocessing:

Download: grch38_chromosome table
https://anvil.terra.bio/#workspaces/t2t-nhgri/t2t-vc-processing/data


```
$ tr '\t' '\n' < ../grch38_chromosome.tsv | grep per_sample > persample.urls
$ cat persample.urls | gsutil -m cp -I .
$ mkdir population superpopulation
$ mv *AFR* *AMR* *EAS* *EUR* *SAS* superpopulation/
$ mv *.stats.txt population/

$ cd population

$ for i in `/bin/ls *.stats.txt`; do chr=`echo $i | tr '._' '  ' | cut -f8 -d' '`; grep -B1 '^PSC' $i | tr ' ' '_' | awk -v chr=$chr '{print chr,$0}' > $i.psc; done

$ head -1 *chr1.*ACB*.psc | sed 's/chr1/chr/' | sed 's/#_PSC/PSC/' | tr ' ' '\t' > header


$ for pop in `ls *chr1.* | grep psc | tr '._' '  ' | awk '{print $13}'`; do echo $pop; (cat header; grep -h -v '#_PSC' *$pop*.stats.txt.psc) | tr ' ' '\t' > $pop.pop.psc; done

```



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

