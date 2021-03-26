version 1.0

workflow t2t_genomics_db {
    input {
        String dbBucket
        String interval
        String chromosome
        String start
        String end
    }

    call generateGenomicsDB {
        input:
            dbBucket = dbBucket,
            interval = interval,
            chromosome = chromosome,
            start = start,
            end = end
    }

    output {
        String done = "done"
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
            --overwrite-existing-genomicsdb-workspace true \
            --genomicsdb-workspace-path "gs://~{dbBucket}/genomics_db/~{interval}" \
            --reader-threads $(nproc) \
            -L "~{chromosome}:~{start}-~{end}" \
            --batch-size 50
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 100 SSD"
        memory: "32G"
        cpu : 4
    }
}