version 1.0

workflow t2t_interval_calling {
    input {
        File refFasta
        File refIndex
        File refDict
        String interval
        String chromosome
        String start
        String end
        String dbBucket
    }

    call genotypeInterval {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            interval = interval,
            chromosome = chromosome,
            start = start,
            end = end,
            dbBucket = dbBucket
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
    }
}

task genotypeInterval {
    input {
        File refFasta
        File refIndex
        File refDict
        String interval
        String chromosome
        String start
        String end
        String dbBucket
    }

    command <<<
        gatk \
            --java-options -Xmx8G \
            GenotypeGVCFs \
            -R "~{refFasta}" \
            -O "~{interval}.margined.genotyped.vcf" \
            -L "~{chromosome}:~{start}-~{end}" \
            -V "gendb.gs://~{dbBucket}/genomics_db/~{interval}"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 10 SSD"
        memory: "12G"
        cpu : 4
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
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 10 SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File genotypedIntervalVCF = "~{interval}.genotyped.vcf"
    }
}