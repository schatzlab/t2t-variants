## Download the samtools stats files from the processing bucket and build a unified matrix across all samples
###############################################################################################################

## Samantha's processing bucket
#sambucket=fc-a5ce9ae4-5914-4e1f-9d78-fb9263a7c5f5

## AnVIL_T2T bucket
anvilbucket=fc-47de7dae-e8e6-429c-b760-b4ba49136eee

## from Mike's t2t-nhgri-schatz-rstudio bucket
mikebucket=fc-675b02d2-2b79-47c7-b788-9a053268c650

## download the metadata table of sample ethnicity, gender
gsutil cp gs://$mikebucket/igsr.metadata.tsv .

## Grab the list of samtools.stats.txt files in the bucket
today=`date +%Y.%m.%d`

## download the datatable of samples and extract samtoolsstats urls
#Rscript -e 'library(AnVIL); args = commandArgs(trailingOnly=TRUE); write.table(avtable("ena_sample", "t2t-nhgri", "t2t-nhgri-processing"), args[1], sep="\t")' $today.ena_sample.samantha.txt
#Rscript -e 'library(AnVIL); args = commandArgs(trailingOnly=TRUE); write.table(avtable("ena_sample", "anvil-datastorage", "AnVIL_T2T"), args[1], sep="\t")' $today.ena_sample.anvil.txt
#cat $today.ena_sample.*.txt | tr '\t' '\n' | grep samtools.stats.txt | tr -d '"' > $today.samtoolsstats.urls

## get the list of samtoolsstatsfiles directly from the bucket
gsutil ls gs://$anvilbucket/samtools_stats/*.samtools.stats.txt > $today.samtoolsstats.urls

## confirm number of samples
wc -l $today.samtoolsstats.urls 

## download the set of samtools stats files
mkdir samtoolsstats
cd samtoolsstats
#cat ../$today.samtoolsstats.urls | parallel -t gsutil cp {} .
cat ../$today.samtoolsstats.urls | gsutil -m cp -I ./

## Now build the data matrix

## make a header file with the row names
(head -1 ../igsr.metadata.tsv | tr '\t' '\n' | tr ' ' '_'; grep ^SN NA18572.samtools.stats.txt | cut -f 2- | cut -f1 -d: | tr ' ' '_') > header

## now build a sample-specific set of files with metadata from igsr.metadata file plus samtools mapping rates
for i in `/bin/ls *.stats.txt`; 
do 
  j=`echo $i | cut -f1 -d'.'`
  grep $j ../igsr.metadata.tsv  | tr '\t' '\n' | tr ' ' '_' > $j.clean
  grep ^SN $i | cut -f2 -d':' | awk '{print $1}' >> $j.clean; 
done

## paste all the files together to build the full matrix
paste header *.clean > ../$today.samtools.stats.all.txt

## copy back to mike's bucket for long term storage
gsutil cp ../$today.samtools.stats.all.txt gs://$mikebucket/
#gsutil cp ../$today.ena_sample.txt gs://$mikebucket/
