version 1.0

workflow variant_recalibration {
    input {
        File refFasta
        File refIndex
        File refDict
        File VCF
        File dbsnp
        File dbsnp_index
        File kg_mills
        File kg_mills_index
        File hapmap
        File hapmap_index
        File kg_omni
        File kg_omni_index
        File kg_snps
        File kg_snps_index
        String chromosome
    }

    call createSNPrecalibration {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            VCF = VCF,
            hapmap = hapmap,
            hapmap_index = hapmap_index,
            kg_omni = kg_omni,
            kg_omni_index = kg_omni_index,
            kg_snps = kg_snps,
            kg_snps_index = kg_snps_index,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            chromosome = chromosome
    }

    call createIndelRecalibration {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            VCF = VCF,
            dbsnp = dbsnp,
            dbsnp_index = dbsnp_index,
            kg_mills = kg_mills,
            kg_mills_index = kg_mills_index,
            chromosome = chromosome
    }

    call applySNPrecalibration {
        input:
            VCF = VCF,
            recal = createSNPrecalibration.recal,
            tranches = createSNPrecalibration.tranches,
            chromosome = chromosome
    }

    call applyIndelRecalibration {
        input:
            VCF = applySNPrecalibration.recalibrated_vcf,
            recal = createIndelRecalibration.recal,
            tranches = createIndelRecalibration.tranches,
            chromosome = chromosome
    }

    call finalizeVCF {
        input:
            VCF = applyIndelRecalibration.recalibrated_vcf
    }

    output {
        File recalibratedVCF = applyIndelRecalibration.recalibrated_vcf
        File recalibratedVCFgz = finalizeVCF.bgzipVCF
        File recalibratedVCFtabix = finalizeVCF.tabix
    }
}

task createSNPrecalibration {
    input {
        File refFasta
        File refIndex
        File refDict
        File VCF
        File hapmap
        File hapmap_index
        File kg_omni
        File kg_omni_index
        File kg_snps
        File kg_snps_index
        File dbsnp
        File dbsnp_index
        String chromosome
    }

    String fastaName = '~{basename(refFasta)}'
    String hapmap_name = '~{basename(hapmap)}'
    String kg_omni_name = '~{basename(kg_omni)}'
    String kg_snps_name = '~{basename(kg_snps)}'
    String dbsnp_name = '~{basename(dbsnp)}'

    command <<<
        mv "~{refFasta}" .
        mv "~{refIndex}" .
        mv "~{refDict}" .

        mv "~{hapmap}" .
        mv "~{hapmap_index}" .
        mv "~{kg_omni}" .
        mv "~{kg_omni_index}" .
        mv "~{kg_snps}" .
        mv "~{kg_snps_index}" .
        mv "~{dbsnp}" .
        mv "~{dbsnp_index}" .

        gatk \
            VariantRecalibrator \
            --java-options -Xmx8G \
            -R "~{fastaName}" \
            -V "~{VCF}" \
            --mode SNP \
            -O "~{chromosome}.recalibrate_snp.recal" \
            --tranches-file "~{chromosome}.recalibrate_snp.tranches" \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "~{dbsnp_name}" \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "~{hapmap_name}" \
            --resource:omni,known=false,training=true,truth=true,prior=12.0 "~{kg_omni_name}" \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 "~{kg_snps_name}" \
            -an QD -an FS -an MQ -an ReadPosRankSum -an MQRankSum -an SOR -an DP \
            -tranche 100.0 -tranche 99.8 -tranche 99.6 -tranche 99.4 -tranche 99.2 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
            --max-gaussians 2
    >>>

    Int diskGb = ceil(20.0 * size(VCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File recal = "~{chromosome}.recalibrate_snp.recal"
        File tranches = "~{chromosome}.recalibrate_snp.tranches"
    }
}

task createIndelRecalibration {
    input {
        File refFasta
        File refIndex
        File refDict
        File VCF
        File kg_mills
        File kg_mills_index
        File dbsnp
        File dbsnp_index
        String chromosome
    }

    String fastaName = '~{basename(refFasta)}'
    String kg_mills_name = '~{basename(kg_mills)}'
    String dbsnp_name = '~{basename(dbsnp)}'

    command <<<
        mv "~{refFasta}" .
        mv "~{refIndex}" .
        mv "~{refDict}" .

        mv "~{kg_mills}" .
        mv "~{kg_mills_index}" .
        mv "~{dbsnp}" .
        mv "~{dbsnp_index}" .

        gatk \
            VariantRecalibrator \
            --java-options -Xmx8G \
            -R "~{fastaName}" \
            -V "~{VCF}" \
            --mode INDEL \
            -O "~{chromosome}.recalibrate_indel.recal" \
            --tranches-file "~{chromosome}.recalibrate_indel.tranches" \
            --resource:mills,known=false,training=true,truth=true,prior=15.0 "~{kg_mills_name}" \
            -an QD -an FS -an ReadPosRankSum -an SOR -an DP \
            -tranche 100.0 -tranche 99.0 -tranche 95.0 -tranche 92.0 -tranche 90.0 \
            --max-gaussians 1
    >>>

    Int diskGb = ceil(10.0 * size(VCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File recal = "~{chromosome}.recalibrate_indel.recal"
        File tranches = "~{chromosome}.recalibrate_indel.tranches"
    }
}

task applySNPrecalibration {
    input {
        File VCF
        File recal
        File tranches
        String chromosome
    }

    String recalName = '~{basename(recal)}'

    command <<<
        mv "~{recal}" .

        gatk \
            IndexFeatureFile \
            --input "~{recalName}"

        gatk \
            ApplyVQSR \
            -V "~{VCF}" \
            -O "~{chromosome}.recalibrated.snp.vcf" \
            -mode SNP \
            --recal-file "~{recalName}" \
            --tranches-file "~{tranches}" \
            -truth-sensitivity-filter-level 99.8
    >>>

    Int diskGb = ceil(10.0 * size(VCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File recalibrated_vcf = "~{chromosome}.recalibrated.snp.vcf"
    }
}

task applyIndelRecalibration {
    input {
        File VCF
        File recal
        File tranches
        String chromosome
    }
    
    String recalName = '~{basename(recal)}'

    command <<<
        mv "~{recal}" .

        gatk \
            IndexFeatureFile \
            --input "~{recalName}"

        gatk \
            ApplyVQSR \
            -V "~{VCF}" \
            -O "~{chromosome}.recalibrated.snp_indel.vcf" \
            -mode INDEL \
            --recal-file "~{recalName}" \
            --tranches-file "~{tranches}" \
            -truth-sensitivity-filter-level 99.0
    >>>

    Int diskGb = ceil(10.0 * size(VCF, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File recalibrated_vcf = "~{chromosome}.recalibrated.snp_indel.vcf"
    }
}

task finalizeVCF {
    input {
        File VCF
    }

    String vcfName = '~{basename(VCF, ".vcf")}'

    command <<<
        bgzip -c "~{VCF}" > "~{vcfName}.vcf.gz"
        tabix -p vcf "~{vcfName}.vcf.gz"
    >>>

    Int diskGb = ceil(2.0 * size(VCF, "G"))

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