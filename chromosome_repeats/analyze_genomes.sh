#!/bin/bash

if [ ! -r work ]
then
	mkdir work
fi

cd work

if [ ! -r hg38 ]
then
	mkdir hg38
fi

cd hg38

if [ ! -r Homo_sapiens_assembly38.fasta ]
then
    echo "downloading hg38"
	#gsutil cp gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/resources/grch38/Homo_sapiens_assembly38.fasta .
fi

../../analyze_uniqueness.sh Homo_sapiens_assembly38.fasta

cd ..

if [ ! -r chm13 ]
then
	mkdir chm13
fi

cd chm13

if [ ! -r t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta ]
then
    echo "downloading chm13"
#	gsutil cp gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/resources/chm13/t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta .
fi

../../analyze_uniqueness.sh t2t-chm13.20200921.withGRCh38chrY.chrEBV.chrYKI270740v1r.fasta
