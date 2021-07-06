version 1.0

workflow filter_af_variants {
    input {
        File inputVCFgz
        File samples
    }

    call bcftools_filter_af {
        input:
            inputVCFgz = inputVCFgz
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = bcftools_filter_af.subsetVCF
    }

    call bcftools_stats {
        input:
            inputVCFgz = bgzip_bcftools_index.bgzipVCF,
            vcfIndex = bgzip_bcftools_index.index,
            samples = samples
    }

    output {
        File subsetVCF = bcftools_filter_af.subsetVCF
        File subset_bgzip = bgzip_bcftools_index.bgzipVCF
        File subset_index = bgzip_bcftools_index.index
        File subset_stats = bcftools_stats.stats
    }
}

task bcftools_filter_af {
    input {
        File inputVCFgz
    }

    String vcfPrefix = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        bcftools view \
            -q 0.495 -Q 0.505 \
            "~{inputVCFgz}" -o tmp.1.vcf

        bcftools annotate \
            -x INFO/NEGATIVE_TRAIN_SITE,INFO/POSITIVE_TRAIN_SITE \
            tmp.1.vcf -o tmp.2.vcf

        bcftools norm \
            -a --atom-overlaps . \
            -m -both \
            tmp.2.vcf -o tmp.3.vcf

        bcftools view \
            -q 0.495 -Q 0.505  \
            tmp.3.vcf \
            -o "~{vcfPrefix}.af_0.5.vcf"
    >>>

    Int diskGb = ceil(5.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/bcftools:v1.12"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File subsetVCF = "~{vcfPrefix}.af_0.5.vcf"
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
    }

    output {
        File bgzipVCF = "~{vcfName}.gz"
        File index = "~{vcfName}.gz.csi"
    }
}

task bcftools_stats {
    input {
        File inputVCFgz
        File vcfIndex
        File samples
    }

    String vcfName = '~{basename(inputVCFgz)}'
    String indexName = '~{basename(vcfIndex)}'
    String vcfPrefix = '~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        mv "~{inputVCFgz}" .
        mv "~{vcfIndex}" .

        bcftools stats \
            -v \
            ./"~{vcfName}" \
            --samples-file "~{samples}" > "~{vcfPrefix}.bcftools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputVCFgz, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File stats = "~{vcfPrefix}.bcftools.stats.txt"
    }
}