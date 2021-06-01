version 1.0

workflow concatenate_files {
    input {
        Array[File] inputFiles
        String filename
    }

    call concatenate_files {
        input:
            inputFiles = inputFiles,
            filename = filename
    }

    output {
        File concatenatedFile = concatenate_files.ofn
    }
}

task concatenate_files {
    input {
        Array[File] inputFiles
        String filename
    }

    command <<<
        cat ~{sep=' ' inputFiles} > "~{filename}"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 20 SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File ofn = "~{filename}"
    }
}