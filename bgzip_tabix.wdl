version 1.0

workflow bgzip_tabix {
    input {
        File inputVCF
    }

    call bgzip_tabix {
        input:
            inputVCF = inputVCF
    }

    output {
        File bgzipVCF = bgzip_tabix.bgzipVCF
        File tabix = bgzip_tabix.tabix
    }
}

task bgzip_tabix {
    input {
        File inputVCF
    }

    String vcfName = '~{basename(inputVCF, ".vcf")}'

    command <<<
        bgzip -c "~{inputVCF}" > "~{vcfName}.vcf.gz"
        tabix -p vcf "~{vcfName}.vcf.gz"
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
        File bgzipVCF = "~{vcfName}.vcf.gz"
        File tabix = "~{vcfName}.vcf.gz.tbi"
    }
}