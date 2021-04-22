version 1.0

workflow gunzip {
    input {
        File inputGZ
    }

    call gunzip {
        input:
            inputGZ = inputGZ
    }

    output {
        File output_file = gunzip.unzipped
    }
}

task gunzip {
    input {
        File inputGZ
    }

    String inputPrefix = '~{basename(inputGZ,".gz")}'

    command <<<
        gunzip -c "~{inputGZ}" > "~{inputPrefix}"
    >>>

    runtime {
        docker : "ubuntu:20.04"
        disks : "local-disk 2000 SSD"
        memory: "4G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File unzipped = "~{inputPrefix}"
    }
}