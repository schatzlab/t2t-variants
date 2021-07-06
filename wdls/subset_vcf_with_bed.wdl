version 1.0

workflow subset_vcf_regions {
    input {
        File inputVCFgz
        File vcfIndex
        File regionsBed
    }

    call subset_bed {
        input:
            inputVCFgz = inputVCFgz,
            vcfIndex = vcfIndex,
            regionsBed = regionsBed
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = subset_bed.subsetVCF
    }

    call bcftools_stats {
        input:
            inputVCFgz = bgzip_bcftools_index.bgzipVCF
    }

    output {
        File subsetVCF = subset_bed.subsetVCF
        File subset_bgzip = bgzip_bcftools_index.bgzipVCF
        File subset_index = bgzip_bcftools_index.index
        File subset_stats = bcftools_stats.stats
    }
}

task subset_bed {
    input {
        File inputVCFgz
        File vcfIndex
        File regionsBed
    }

    String vcfPrefix = '~{basename(inputVCFgz,".vcf.gz")}'
    String bedPrefix = '~{basename(regionsBed,".bed")}'

    command <<<
        bcftools view \
            -R "~{regionsBed}" \
            "~{inputVCFgz}" > "~{vcfPrefix}.~{bedPrefix}.unsorted.vcf"

        mkdir -p ./sort_tmp/
        bcftools sort -T ./sort_tmp/ "~{vcfPrefix}.~{bedPrefix}.unsorted.vcf" > "~{vcfPrefix}.~{bedPrefix}.vcf"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 2000 SSD"
        memory: "8G"
        cpu : 4
    }

    output {
        File subsetVCF = "~{vcfPrefix}.~{bedPrefix}.vcf"
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
        preemptible: 2
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
        bcftools stats \
            -v \
            "~{inputVCFgz}" > "~{vcfName}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
        preemptible: 2
        maxRetries: 2
    }

    output {
        File stats = "~{vcfName}.bcftools.stats.txt"
    }
}