tr '\t' '\n' < population.tsv | grep ancestry.txt > ancestry.urls
mkdir data
cd data/
cat ../ancestry.urls | gsutil -m cp -I .
/bin/ls *.ancestry.txt | parallel -t ../tally_counts.sh {}
for i in `/bin/ls *.txt.stats`; do awk '{if (NR==1){split($1, a, "\.")} else {print a[1],a[2],$1,$6,$7,$8}}' $i > $i.clean; done
cat *.clean > ../summary.all.txt
