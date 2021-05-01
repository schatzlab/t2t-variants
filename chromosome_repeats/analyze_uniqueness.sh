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

## Process each chromosome, skipping random, HLA, chrUn, alt
#for chr in `grep -v random summary  | grep -v HLA | grep -v chrUn | grep -v alt`
for chr in `grep chr21.fa summary`
do
  echo "Processing $chr..."

  short=`echo $chr | cut -f 2 -d'/' | cut -f1 -d '.'`
  echo "k	$short" > $chr.uniqueness

  for k in 25 50 100 150 200 250
  do
    if [ ! -r $chr.$k.kmc_pre ]
    then
      echo "indexing $chr $k..."
      $KMC -t10 -k$k -ci1 -cs10000000 -fm $chr $chr.$k . >& $chr.$k.kmc.log
    fi

    if [ ! -r $chr.$k.histo ]
    then
      echo "tally $chr $k..."
      $KMCTOOLS transform $chr.$k histogram $chr.$k.histo >& $chr.$k.histo.log
    fi

    head -1 $chr.$k.histo | awk -v k=$k '{print k"\t"$2}' >> $chr.uniqueness
  done

  if [ ! -r $chr.badchar ]
  then
    echo "counting bad characters"
    grep -v '>' $chr | tr -d 'ACGTacgt' | fold -w1 | grep -v '^$' | wc -l > $chr.badchar
  fi

  head $chr.badchar | awk '{print "N\t"$1}' >> $chr.uniqueness

  grep $short'	' $allfasta.fai | awk '{print "all\t"$2}' >> $chr.uniqueness
done

paste chr*.uniqueness > $allfasta.uniqueness
