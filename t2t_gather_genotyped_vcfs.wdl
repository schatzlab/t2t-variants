version 1.0

workflow t2t_gather_all_VCFs {
    input {
        Array[File] VCFs
    }

    call gatherVCFs {
        input:
            VCFs = VCFs
    }

    output {
        File genotypedVCF = gatherVCFs.genotypedVCF
    }
}

task gatherVCFs {
    input {
        Array[File] VCFs
    }

    command <<<
        gatk \
            GatherVcfs
            -I ~{sep=' ' VCFs} \
            -O "1kgp_3202_samples.genotyped.vcf"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 100 SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File genotypedVCF = "1kgp_3202_samples.genotyped.vcf"
    }
}