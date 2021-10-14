version 1.0

workflow filter_ac_variants {
    input {
        File inputVCFgz
    }

    call bcftools_view_ac {
        input:
            inputVCFgz = inputVCFgz
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = bcftools_view_ac.subsetVCF
    }

    call bcftools_stats {
        input:
            inputVCFgz = bgzip_bcftools_index.bgzipVCF
    }

    output {
        File subsetVCF = bcftools_view_ac.subsetVCF
        File subset_bgzip = bgzip_bcftools_index.bgzipVCF
        File subset_index = bgzip_bcftools_index.index
        File subset_stats = bcftools_stats.stats
    }
}

task bcftools_view_ac {
    input {
        File inputVCFgz
    }

    String vcfPrefix = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools view \
            --min-ac 1 \
            "~{inputVCFgz}" > "~{vcfPrefix}.variants.vcf"
    >>>

    Int diskGb = ceil(100.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "16G"
        cpu : 4
    }

    output {
        File subsetVCF = "~{vcfPrefix}.variants.vcf"
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
        maxRetries: 2
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
        bcftools stats "~{inputVCFgz}" > "~{vcfName}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 2
    }

    output {
        File stats = "~{vcfName}.bcftools.stats.txt"
    }
}