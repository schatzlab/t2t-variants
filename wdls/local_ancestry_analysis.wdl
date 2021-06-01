version 1.0

workflow local_ancestry_analysis {
    input {
        File inputVCF
        File windowsBed
        File ancestryBed
        String population
        String chromosome
    }

    call vcf_coverage_windows {
        input:
            inputVCF = inputVCF,
            windowsBed = windowsBed,
            population = population,
            chromosome = chromosome
    }

    call separate_by_ancestry {
        input:
            countsFile = countsFile,
            ancestryBed = ancestryBed,
            population = population,
            chromosome = chromosome
    }

    output {
        File countsFile = vcf_coverage_windows.countsFile
        File ancestryCounts = separate_by_ancestry.ancestryCounts
    }
}

task vcf_coverage_windows {
    input {
        File inputVCF
        File windowsBed
        String population
        String chromosome
    }

    command <<<
        awk '/~{chromosome}\t/' "~{windowsBed}" | sort -k2,2n > "~{chromosome}".bed
        bedtools coverage \
            -a "~{chromosome}".bed \
            -b "~{inputVCF}" \
            -sorted \
            -counts > "~{chromosome}.~{population}.counts.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ~{diskGb} SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File countsFile = "~{chromosome}.~{population}.counts.txt"
    }
}

task separate_by_ancestry {
    input {
        File countsFile
        File ancestryBed
        String population
        String chromosome
    }

    command <<<
        # get all windows that overlap
        bedtools \
            intersect \
            -a "~{countsFile}" \
            -b "~{ancestryBed}" \
            -wa -wb | cut -f1,2,3,4,8 | sort | uniq > all.txt

bedtools intersect -a chr1.EUR.counts.txt -b updated_grch38_ancestry.bed -wa -wb | cut -f1,2,3,4,8 | sort | uniq > all.txt

        # get windows that 100% overlap with one population
        bedtools \
            intersect \
            -a "~{countsFile}" \
            -b "~{ancestryBed}" \
            -wa -wb -f 1 | cut -f1,2,3,4,8 | sort | uniq > full.txt

bedtools intersect -a chr1.EUR.counts.txt -b updated_grch38_ancestry.bed -wa -wb -f 1 | cut -f1,2,3,4,8 | sort | uniq > full.txt

        # convert windows that partly overlap to read as MULTI-population
        grep -F -x -v -f full.txt all.txt | cut -f1,2,3,4 | sort -k1,3 | uniq > tmp.txt
        awk -v OFS='\t' '{print $1, $2, $3, $4, "MULTI"}' tmp.txt | sort -k1,3 | uniq > multi.txt

        # get repeated lines (lines assigned to more than one population)
        cat full.txt multi.txt > tmp_merged.txt
        cat all.txt | cut -f1,2,3 | uniq -d > repeats.txt

        # prioritize Neanderthal or MULTI if assigned to more than one pop
        while IFS= read -r line; do
            grep "${line}" all.txt >> all_repeats.txt
        done < repeats.txt
        grep "Neanderthal" all_repeats.txt > repeats_to_include.txt
        grep "MULTI" all_repeats.txt >> repeats_to_include.txt

        # remove non-Neanderthal or non-MULTI lines
        grep -F -x -v -f repeats_to_include.txt all_repeats.txt > repeats_to_exclude.txt

        # produce final output
        grep -F -x -v -f repeats_to_exclude.txt tmp_merged.txt | sort -n -k2 > "~{chromosome}.~{population}.counts.ancestry.txt"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 10 SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File ancestryCounts = "~{chromosome}.~{population}.counts.ancestry.txt"
    }
}