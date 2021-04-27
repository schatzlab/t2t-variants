# AF Stats

## preprocess the data
```
gsutil -m cp gs://fc-51aefb1c-4e8e-4dcb-a59c-62e318ea351a/*_autosomes_summarized.txt .
for i in `/bin/ls *_summarized.txt`; do echo $i; grep -B1 '^AF' $i | sed 's/^# //' > $i.af; done
```

## Make the plots

See `afstats.R`



