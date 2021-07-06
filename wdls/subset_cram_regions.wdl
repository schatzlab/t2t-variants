version 1.0

workflow subset_cram_regions {
    input {
        File cram
        File regionsBed
    }

    call subset_cram {
        input:
            cram = cram,
            regionsBed = regionsBed
    }

    call index_bam {
        input:
            inputBAM = subset_cram.subsetBAM
    }

    output {
        File outputBAM = subset_cram.subsetBAM
        File outputBAI = index_bam.bamIndex
    }
}

task subset_cram {
    input {
        File cram
        File regionsBed
    }

    String cramPrefix = '~{basename(cram,".cram")}'
    String bedPrefix = '~{basename(regionsBed,".bed")}'

    command <<<
        samtools view \
            -L "~{regionsBed}" \
            "~{cram}" -o "~{cramPrefix}.~{bedPrefix}.bam"

        samtools sort \
            "~{cramPrefix}.~{bedPrefix}.bam" \
            -@ "$(nproc)" \
            -o "~{cramPrefix}.~{bedPrefix}.sorted.bam"
    >>>

    Int diskGb = ceil(5.0 * size(cram, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 12
        preemptible: 2
        maxRetries: 2
    }

    output {
        File subsetBAM = "~{cramPrefix}.~{bedPrefix}.sorted.bam"
    }
}

task index_bam {
    input {
        File inputBAM
    }

    String bamName = '~{basename(inputBAM)}'

    command <<<
        samtools index -@ "$(nproc)" "~{inputBAM}" "~{bamName}.bai"
    >>>

    Int diskGb = ceil(2.0 * size(inputBAM, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 2
        maxRetries: 2
    }

    output {
        File bamIndex = "~{bamName}.bai"
    }
}
