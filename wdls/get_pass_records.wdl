version 1.0

workflow get_pass_records {
    input {
        File inputVCFgz
    }

    call get_pass_records {
        input:
            inputVCFgz = inputVCFgz
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = get_pass_records.passingVCF
    }

    call bcftools_stats {
        input:
            inputVCFgz = bgzip_bcftools_index.bgzipVCF
    }

    output {
        File passVCF = get_pass_records.passingVCF
        File pass_bgzip = bgzip_bcftools_index.bgzipVCF
        File pass_index = bgzip_bcftools_index.index
        File pass_stats = bcftools_stats.stats
    }
}

task get_pass_records {
    input {
        File inputVCFgz
    }

    String vcfPrefix = '~{basename(inputVCFgz, ".vcf.gz")}'

    command <<<
        bcftools view -f PASS "~{inputVCFgz}" > "~{vcfPrefix}.pass.vcf"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 2500 SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File passingVCF = "~{vcfPrefix}.pass.vcf"
    }
}

task bgzip_bcftools_index {
    input {
        File inputVCF
    }

    String vcfName = '~{basename(inputVCF)}'

    command <<<
        bgzip -@ "$(nproc)" -c "~{inputVCF}" > "~{vcfName}.gz"
        bcftools index --threads "$(nproc)" "~{vcfName}.gz" -o "~{vcfName}.gz.csi"
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
        File bgzipVCF = "~{vcfName}.gz"
        File index = "~{vcfName}.gz.csi"
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
    }

    String vcfName = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools stats -v "~{inputVCFgz}" > "~{vcfName}.bcftools.stats.txt"
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
        File stats = "~{vcfName}.bcftools.stats.txt"
    }
}