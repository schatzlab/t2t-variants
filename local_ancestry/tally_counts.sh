#!/bin/bash

file=$1

if [ ! -r $file ]
then
  echo "USAGE: tally file.ancestry.txt"
  exit 1
fi

if [ ! -r $file.tally ]
then 
  echo "tally $file"
  awk '{print $5}' $file | hashcount > $file.tally
fi

set -xv

echo $file > $file.stats

for i in `awk '{print $1}' $file.tally`
do
  echo $i
  echo -n "$i " >> $file.stats
  grep $i $file | stats -f 4 >> $file.stats
done

