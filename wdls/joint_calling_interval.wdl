version 1.0

workflow joint_calling_interval {
    input {
        File refFasta
        File refIndex
        File refDict
        File genomicsDBtar
        String interval
        String chromosome
        String start
        String marginedStart
        String end
        String marginedEnd
    }

    call genotypeInterval {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            genomicsDBtar = genomicsDBtar,
            interval = interval,
            chromosome = chromosome,
            start = marginedStart,
            end = marginedEnd
    }

    call trimMarginVCF {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            intervalVCF = genotypeInterval.intervalVCF,
            interval = interval,
            chromosome = chromosome,
            start = start,
            end = end
    }

    output {
        File genotypedIntervalVCF = trimMarginVCF.genotypedIntervalVCF
        File genotypedIntervalVCFgz = trimMarginVCF.genotypedIntervalVCFgz
        File genotypedIntervalVCFtabix = trimMarginVCF.genotypedIntervalVCFtabix
    }
}

task genotypeInterval {
    input {
        File refFasta
        File refIndex
        File refDict
        File genomicsDBtar
        String interval
        String chromosome
        String start
        String end
    }

    command <<<
        tar -xf "~{genomicsDBtar}"

        gatk \
            --java-options -Xmx8G \
            GenotypeGVCFs \
            -R "~{refFasta}" \
            -O "~{interval}.margined.genotyped.vcf" \
            -L "~{chromosome}:~{start}-~{end}" \
            -V "gendb://~{interval}" \
            --only-output-calls-starting-in-intervals \
            --merge-input-intervals
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 50 SSD"
        memory: "12G"
        cpu : 4
        preemptible: 3
        maxRetries: 3
    }

    output {
        File intervalVCF = "~{interval}.margined.genotyped.vcf"
    }
}

task trimMarginVCF {
    input {
        File refFasta
        File refIndex
        File refDict
        File intervalVCF
        String interval
        String chromosome
        String start
        String end
    }

    String vcfName = '~{basename(intervalVCF, ".vcf")}'

    command <<<
        bgzip -c "~{intervalVCF}" > "~{vcfName}.vcf.gz"
        tabix -p vcf "~{vcfName}.vcf.gz"

        gatk \
            --java-options -Xmx8G \
            SelectVariants \
            -R "~{refFasta}" \
            -V "~{vcfName}.vcf.gz" \
            -L "~{chromosome}:~{start}-~{end}" \
            -O "~{interval}.genotyped.vcf"

        bgzip -c "~{interval}.genotyped.vcf" > "~{interval}.genotyped.vcf.gz"
        tabix -p vcf "~{interval}.genotyped.vcf.gz"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 50 SSD"
        memory: "12G"
        cpu : 4
        preemptible: 3
        maxRetries: 3
    }

    output {
        File genotypedIntervalVCF = "~{interval}.genotyped.vcf"
        File genotypedIntervalVCFgz = "~{interval}.genotyped.vcf.gz"
        File genotypedIntervalVCFtabix = "~{interval}.genotyped.vcf.gz.tbi"
    }
}