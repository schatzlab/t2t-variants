version 1.0

workflow bcftools_index {
    input {
        File inputVCF
    }

    call bcftools_index {
        input:
            inputVCF = inputVCF
    }

    output {
        File index = bcftools_index.index
    }
}

task bcftools_index {
    input {
        File inputVCF
    }

    String vcfName = '~{basename(inputVCF)}'

    command <<<
        bcftools index --threads "$(nproc)" "~{inputVCF}" -o "~{vcfName}.csi"
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
        File index = "~{vcfName}.csi"
    }
}