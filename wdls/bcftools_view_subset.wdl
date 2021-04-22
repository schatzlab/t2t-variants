version 1.0

workflow bcftools_view_subset {
    input {
        File inputVCF
        File sample_list
        String population
    }

    call bcftools_view_subset {
        input:
            inputVCF = inputVCF,
            sample_list = sample_list,
            population = population
    }

    output {
        File subsetVCF = bcftools_view_subset.subsetVCF
    }
}

task bcftools_view_subset {
    input {
        File inputVCF
        File sample_list
        String population
    }

    String vcfPrefix = '~{basename(inputVCF,".vcf")}'

    command <<<
        bcftools view \
            -S "~{sample_list}" \
            --force-samples \
            "~{inputVCF}" > "~{vcfPrefix}.~{population}.vcf"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File subsetVCF = "~{vcfPrefix}.~{population}.vcf"
    }
}