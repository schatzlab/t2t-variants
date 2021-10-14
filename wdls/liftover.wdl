version 1.0

workflow liftoverVCF {
    input {
        File inputVCF
        File chain
        File targetFasta
        File targetDict
    }

    call liftoverVCF {
        input:
            inputVCF = inputVCF,
            chain = chain,
            targetFasta = targetFasta,
            targetDict = targetDict
    }

    call bgzip_bcftools_index as lifted {
        input:
            inputVCF = liftoverVCF.liftedVCF
    }

    call bgzip_bcftools_index as rejected {
        input:
            inputVCF = liftoverVCF.rejectedVCF
    }

    call bcftools_stats as stats_lifted {
        input:
            inputVCFgz = lifted.bgzipVCF
    }

    call bcftools_stats as stats_rejected {
        input:
            inputVCFgz = rejected.bgzipVCF
    }

    output {
        File liftedVCF = liftoverVCF.liftedVCF
        File lifted_bgzip = lifted.bgzipVCF
        File lifted_index = lifted.index
        File lifted_stats = stats_lifted.stats

        File rejectedVCF = liftoverVCF.rejectedVCF
        File rejected_bgzip = rejected.bgzipVCF
        File rejected_index = rejected.index
        File rejected_stats = stats_rejected.stats
    }
}

task liftoverVCF {
    input {
        File inputVCF
        File chain
        File targetFasta
        File targetDict
    }

    String vcfBase='~{basename(inputVCF,".vcf")}'
    String target='~{basename(targetFasta,".fasta")}'

    command <<<
        gatk \
            --java-options -Xmx64G \
            LiftoverVcf \
            -I "~{inputVCF}" \
            --CHAIN "~{chain}" \
            -R "~{targetFasta}" \
            -O "~{vcfBase}.~{target}.liftover.vcf" \
            --REJECT "~{vcfBase}.~{target}.rejected.vcf" \
            --WARN_ON_MISSING_CONTIG \
            --MAX_RECORDS_IN_RAM 100000
    >>>

    Int diskGb = ceil(10.0 * size(inputVCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "128G"
        cpu : 12
    }

    output {
        File liftedVCF = "~{vcfBase}.~{target}.liftover.vcf"
        File rejectedVCF = "~{vcfBase}.~{target}.rejected.vcf"
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