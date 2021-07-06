version 1.0

workflow mendelian_analysis {
    input {
        File refFasta
        File refIndex
        File refDict
        File inputVCF
        File sample_list
        File pedigree
        String population
    }

    call bcftools_view_subset {
        input:
            inputVCF = inputVCF,
            sample_list = sample_list,
            population = population
    }

    call bgzip_bcftools_index {
        input:
            inputVCF = bcftools_view_subset.subsetVCF
    }

    call mendelian_analysis {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            inputVCFgz = bgzip_bcftools_index.bgzipVCF,
            pedigree = pedigree
    }

    output {
        File subsetVCF = bcftools_view_subset.subsetVCF
        File subset_bgzip = bgzip_bcftools_index.bgzipVCF
        File subset_index = bgzip_bcftools_index.index
        File mendelian_violations = mendelian_analysis.MVtable
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
    }

    output {
        File subsetVCF = "~{vcfPrefix}.~{population}.vcf"
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

task mendelian_analysis {
    input {
        File refFasta
        File refIndex
        File refDict
        File inputVCFgz
        File pedigree
    }

    String vcfBase='~{basename(inputVCFgz,".vcf.gz")}'

    command <<<
        gatk IndexFeatureFile -I "~{inputVCFgz}" -O "~{inputVCFgz}.tbi"

        gatk \
            --java-options -Xmx16G \
            VariantEval \
            -R "~{refFasta}" \
            --eval "~{inputVCFgz}" \
            --pedigree "~{pedigree}" \
            -no-ev -no-st --lenient \
            -ST Family \
            -EV MendelianViolationEvaluator \
            -O "~{vcfBase}.MVs.byFamily.table"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 1000 SSD"
        memory: "128G"
        cpu : 8
    }

    output {
        File MVtable = "~{vcfBase}.MVs.byFamily.table"
    }
}