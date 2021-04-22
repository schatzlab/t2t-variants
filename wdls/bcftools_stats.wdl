version 1.0

workflow bcftools_stats {
    input {
        File inputVCFgz
        File sample_list
        String population
    }

    call bcftools_stats {
        input:
            inputVCFgz = inputVCFgz,
            sample_list = sample_list,
            population = population
    }

    output {
        File stats = bcftools_stats.stats
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
        File sample_list
        String population
    }

    String vcfPrefix = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools stats \
            "~{inputVCFgz}" \
            --samples-file "~{sample_list}" > "~{vcfPrefix}.~{population}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File stats = "~{vcfPrefix}.~{population}.bcftools.stats.txt"
    }
}