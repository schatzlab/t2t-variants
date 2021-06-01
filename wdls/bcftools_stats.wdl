version 1.0

workflow bcftools_stats {
    input {
        File inputVCFgz
        File vcfIndex
        File regionsBed
    }

    call bcftools_stats {
        input:
            inputVCFgz = inputVCFgz,
            vcfIndex = vcfIndex,
            regionsBed = regionsBed
    }

    output {
        File stats = bcftools_stats.stats
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
        File vcfIndex
        File regionsBed
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'
    String bedPrefix = '~{basename(regionsBed,".bed")}'

    command <<<
        bcftools stats \
            -v \
            --regions-file "~{regionsBed}" \
            "~{inputVCFgz}" > "~{vcfName}.~{bedPrefix}.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File stats = "~{vcfName}.~{bedPrefix}.stats.txt"
    }
}