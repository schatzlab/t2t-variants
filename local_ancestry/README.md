tr '\t' '\n' < population.tsv | grep ancestry.txt > ancestry.urls
mkdir data
cd data/
cat ../ancestry.urls | gsutil -m cp -I .
