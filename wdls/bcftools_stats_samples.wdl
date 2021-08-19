version 1.0

workflow bcftools_stats {
    input {
        File inputVCFgz
        File samples
    }

    call bcftools_stats {
        input:
            inputVCFgz = inputVCFgz,
            samples = samples
    }

    output {
        File stats = bcftools_stats.stats
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
        File samples
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools stats \
            -v \
            --samples-file "~{samples}" \
            "~{inputVCFgz}" > "~{vcfName}.per_sample.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File stats = "~{vcfName}.per_sample.stats.txt"
    }
}