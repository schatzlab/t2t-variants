#!/bin/bash

KMC=kmc
KMCTOOLS=kmc_tools

if [ $# -lt 1 ]
then
  echo "USAGE: analyze_uniqueness.sh path/to/fasta"
  exit
fi

allfasta=$1

echo "Analyzing $allfasta"

if [ ! -r $allfasta.fai ]
then
  echo "indexing $allfasta"
  samtools faidx $allfasta
fi

if [ ! -r explode.done ]
then
  echo "exploding fasta"
  ~/build/sge_mummer/explode_fasta.pl $allfasta . && touch explode.done
fi

if [ ! -r $allfasta.autosomesX.fa ]
then
  echo "concat autosomesX..."

  for i in `seq 1 22`
  do
    cat chr$i.fa >> $allfasta.autosomesX.fa
  done

  cat chrX.fa >> $allfasta.autosomesX.fa
fi

if [ ! -r "$allfasta.autosomesX.fa.fai" ]
then
  echo "indexing autosomesX..."
  samtools faidx $allfasta.autosomesX.fa
fi

## Process each chromosome, skipping random, HLA, chrUn, alt
## Note autosomeX must be considered before individual chromosomes

#for chr in "./chr20.fa" "./chr21.fa"
for chr in `echo "./$allfasta.autosomesX.fa"; grep -v random summary  | grep -v HLA | grep -v chrUn | grep -v alt | grep -v chrY | grep -v chrM | grep -v chrEBV`
do
  echo "Processing $chr..."

  short=`echo $chr | cut -f 2 -d'/' | cut -f1 -d '.'`
  echo "k	$short" > $chr.uniqueness

  for k in 25 50 100 250
  do
    if [ ! -r $chr.$k.kmc_pre ]
    then
      echo "  indexing $chr $k..."
      $KMC -t10 -k$k -ci1 -cs10000000 -fm $chr $chr.$k . >& $chr.$k.kmc.log
    fi

    if [ ! -r $chr.$k.histo ]
    then
      echo "  tally $chr $k..."
      $KMCTOOLS -t10 transform $chr.$k histogram $chr.$k.histo >& $chr.$k.histo.log
    fi

    if [[ $chr == "./$allfasta.autosomesX.fa" ]]
    then
      head -1 $chr.$k.histo | awk -v k=$k '{print k"\t"$2}' >> $chr.uniqueness
    else

      if [ ! -r $chr.$k.autosomeX.kmc_pre ]
      then
        echo "  intersecting $chr.$k with autosomeX..."
        $KMCTOOLS -t10 simple $chr.$k $allfasta.autosomesX.fa.$k intersect $chr.$k.autosomeX -ocright >& $chr.$k.autosomeX.intersect.log
      fi

      if [ ! -r $chr.$k.autosomeX.histo ]
      then
        echo "  tally $chr.$k.autosomeX..."
        $KMCTOOLS -t10 transform $chr.$k.autosomeX histogram $chr.$k.autosomeX.histo >& $chr.$k.autosomeX.histo.log
      fi

      head -1 $chr.$k.autosomeX.histo | awk -v k=$k '{print k"\t"$2}' >> $chr.uniqueness
    fi
  done

  if [ ! -r $chr.badchar ]
  then
    echo "  counting bad characters"
    grep -v '>' $chr | tr -d 'ACGTacgt' | fold -w1 | grep -v '^$' | wc -l > $chr.badchar
  fi

  head $chr.badchar | awk '{print "N\t"$1}' >> $chr.uniqueness

  if [[ $chr == "./$allfasta.autosomesX.fa" ]]
  then
    awk '{s+=$2}END{print "all\t"s}' $allfasta.autosomesX.fa.fai >> $chr.uniqueness
  else
    grep $short'	' $allfasta.fai | awk '{print "all\t"$2}' >> $chr.uniqueness
  fi
done




echo "building final table"

awk '{print $1}' chr1.fa.uniqueness > header

for i in `seq 1 22; echo X`
do
  awk '{print $2}' chr$i.fa.uniqueness > chr$i.fa.uniqueness.clean
done

paste header \
  chr1.fa.uniqueness.clean \
  chr2.fa.uniqueness.clean \
  chr3.fa.uniqueness.clean \
  chr4.fa.uniqueness.clean \
  chr5.fa.uniqueness.clean \
  chr6.fa.uniqueness.clean \
  chr7.fa.uniqueness.clean \
  chr8.fa.uniqueness.clean \
  chr9.fa.uniqueness.clean \
  chr10.fa.uniqueness.clean \
  chr11.fa.uniqueness.clean \
  chr12.fa.uniqueness.clean \
  chr13.fa.uniqueness.clean \
  chr14.fa.uniqueness.clean \
  chr15.fa.uniqueness.clean \
  chr16.fa.uniqueness.clean \
  chr17.fa.uniqueness.clean \
  chr18.fa.uniqueness.clean \
  chr19.fa.uniqueness.clean \
  chr20.fa.uniqueness.clean \
  chr21.fa.uniqueness.clean \
  chr22.fa.uniqueness.clean \
  chrX.fa.uniqueness.clean  > $allfasta.uniqueness

