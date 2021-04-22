version 1.0

workflow downloadAspera {
    input {
        String accessionId
    }

    call downloadFileAspera {
        input:
            accessionId = accessionId
    }

    output {
        File fastq1 = downloadFileAspera.fastq1
        File fastq2 = downloadFileAspera.fastq2
    }
}

task downloadFileAspera {
    input {
        String accessionId
    }

    command <<<
        enaDataGet -as /root/aspera_settings.ini --format fastq "~{accessionId}"

        mv "~{accessionId}/~{accessionId}_1.fastq.gz" .
        mv "~{accessionId}/~{accessionId}_2.fastq.gz" .
    >>>

    runtime {
        docker : "szarate/ascp"
        disks : "local-disk 60 SSD"
        memory: "4G"
        cpu : 4
        preemptible: 3
        maxRetries: 3
    }

    output {
        File fastq1 = "~{accessionId}_1.fastq.gz"
        File fastq2 = "~{accessionId}_2.fastq.gz"
    }
}