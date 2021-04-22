version 1.0

workflow bcftools_stats {
    input {
        File inputVCFgz
    }

    call bcftools_stats {
        input:
            inputVCFgz = inputVCFgz
    }

    output {
        File stats = bcftools_stats.stats
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools stats "~{inputVCFgz}" > "~{vcfName}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File stats = "~{vcfName}.bcftools.stats.txt"
    }
}