version 1.0

workflow t2t_genomics_db {
    input {
        File starterVCF
        File sampleMappings
        Boolean updateDB
        String chromosome
        String dbBucket
    }

    if (!updateDB) {
        call startGenomicsDB {
            input:
                starterVCF = starterVCF,
                chromosome = chromosome,
                dbBucket = dbBucket
        }
    }

    if (updateDB) {
        call updateGenomicsDB {
            input:
                sampleMappings = sampleMappings,
                chromosome = chromosome,
                dbBucket = dbBucket
        }
    }
}

task startGenomicsDB {
    input {
        File starterVCF
        String chromosome
        String dbBucket
    }

    String vcfName = '~{basename(starterVCF,".vcf")}'

    command <<<
        bgzip -c "~{starterVCF}" > "~{vcfName}".vcf.gz
        tabix -p vcf "~{vcfName}".vcf.gz

        gatk \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            GenomicsDBImport \
            -V "~{vcfName}".vcf.gz \
            --genomicsdb-workspace-path "gs://~{dbBucket}/genomics_db/~{chromosome}" \
            --reader-threads $(nproc) \
            -L "~{chromosome}"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 10 SSD"
        memory: "12G"
        cpu : 16
    }
}

task updateGenomicsDB {
    input {
        File sampleMappings
        String chromosome
        String dbBucket
    }

    command <<<
        gatk \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            GenomicsDBImport \
            --sample-name-map "~{sampleMappings}" \
            --genomicsdb-update-workspace-path "gs://~{dbBucket}/genomics_db/~{chromosome}" \
            --reader-threads $(nproc) \
            -L "~{chromosome}"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 10 SSD"
        memory: "12G"
        cpu : 16
    }
}