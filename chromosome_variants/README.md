## download the hg38 summary stats files

```
$ gsutil cp gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/grch38_data/singleton/*.stats.txt grch38/
$ (echo -e "chr\tsnps\tindels"; for i in `ls *.txt | sort -V `; do c=`echo $i | cut -f1 -d'.' | rev | cut -f1 -d'_' | rev`; grep '^AF' $i | awk -v chr=$c '{if ($3 > 0){s += $4; i += $7}} END{print chr"\t"s"\t"i}'; done) > hg38.chr.stats
```


## download the chm13 summary stats files

### first "Download all rows" from the "chromosome" table in https://anvil.terra.bio/#workspaces/t2t-nhgri/t2t-vc-processing/data

```
$ awk '{print $172}' chromosome.tsv  | skip.pl 1 > stats.urls
$ cat stats.urls | gsutil -m cp -I .

$ (echo -e "chr\tsnps\tindels"; for i in `ls *.txt | sort -V `; do c=`echo $i | cut -f 2 -d '.' `; grep '^AF' $i | awk -v chr=$c '{if ($3 > 0){s += $4; i += $7}} END{print chr"\t"s"\t"i}'; done) > chm13.chr.stats
```


