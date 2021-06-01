version 1.0

workflow mendelian_analysis {
    input {
        File refFasta
        File refIndex
        File refDict
        File inputVCFgz
        File pedigree
    }

    call mendelian_analysis {
        input:
            refFasta = refFasta,
            refIndex = refIndex,
            refDict = refDict,
            inputVCFgz = inputVCFgz,
            pedigree = pedigree
    }

    output {
        File mendelian_violations = mendelian_analysis.MVtable
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