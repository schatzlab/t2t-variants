#!/bin/sh

cut -f 4 igsr.metadata.tsv | sort | uniq > pops
cut -f 6 igsr.metadata.tsv | sort | uniq > super_pops

while read -r pop; do
    grep "${pop}" igsr.metadata.tsv | cut -f 1 > "${pop}".txt
    gsutil cp "${pop}".txt gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/resources/1kgp/"${pop}".txt
    echo "{\"subset_vcf_sample_list.inputVCFgz\":\"\${this.singleton_bgzip}\",\"subset_vcf_sample_list.population\":\"${pop}_singleton\",\"subset_vcf_sample_list.sample_list\":\"gs://fc-47de7dae-e8e6-429c-b760-b4ba49136eee/resources/1kgp/${pop}.txt\"}" > "${pop}"_input.json
    echo "{\"subset_vcf_sample_list.subsetVCF\":\"\${this.${pop}_singleton_vcf}\",\"subset_vcf_sample_list.subset_bgzip\":\"\${this.${pop}_singleton_bgzip}\",\"subset_vcf_sample_list.subset_index\":\"\${this.${pop}_singleton_index}\",\"subset_vcf_sample_list.subset_stats\":\"\${this.${pop}_singleton_stats}\"}" > "${pop}"_output.json
done < super_pops