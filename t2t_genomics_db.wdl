version 1.0

workflow t2t_genomics_db {
    input {
        String dbBucket
        String interval
        String chromosome
        String marginedStart
        String marginedEnd
    }

    call generateGenomicsDB {
        input:
            dbBucket = dbBucket,
            interval = interval,
            chromosome = chromosome,
            start = marginedStart,
            end = marginedEnd
    }

    output {
        File genomicsDBtar = generateGenomicsDB.genomicsDB
    }
}

task generateGenomicsDB {
    input {
        String dbBucket
        String interval
        String chromosome
        String start
        String end
    }

    command <<<
        gatk \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            GenomicsDBImport \
            --sample-name-map "gs://~{dbBucket}/sample_maps/~{chromosome}_sample_map.tsv" \
            --genomicsdb-workspace-path "./~{interval}" \
            --reader-threads $(nproc) \
            -L "~{chromosome}:~{start}-~{end}" \
            --batch-size 50
        
        tar -cf "~{interval}.tar" "./~{interval}"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 100 SSD"
        memory: "32G"
        cpu : 4
    }

    output {
        File genomicsDB = "~{interval}.tar"
    }
}