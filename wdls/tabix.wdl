version 1.0

workflow tabix {
    input {
        File inputVCFgz
    }

    call tabix {
        input:
            inputVCFgz = inputVCFgz
    }

    output {
        File tbi = tabix.index
    }
}

task tabix {
    input {
        File inputVCFgz
    }

    String vcfName = '~{basename(inputVCFgz)}'

    command <<<
        mv "~{inputVCFgz}" .
        tabix -p vcf ./"~{vcfName}"
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
        File index = "~{vcfName}.tbi"
    }
}