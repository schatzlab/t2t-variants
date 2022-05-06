version 1.0

workflow t2t_concat_vcfs_chromosome {
    input {
        Array[File] VCFs
        Array[File] indexes
        String chromosome
    }

    call concatChromosome {
        input:
            VCFs = VCFs,
            indexes = indexes,
            chromosome = chromosome
    }

    output {
        File chromosomeVCF = concatChromosome.chromosomeVCF
        File chromosomeVCF_gz = concatChromosome.chromosomeVCF_gz
        File chromosomeVCF_tbi = concatChromosome.chromosomeVCF_tbi
    }
}

task concatChromosome {
    input {
        Array[File] VCFs
        Array[File] indexes
        String chromosome
    }

    command <<<
        bcftools concat -a -D ~{sep=' ' VCFs} --threads "$(nproc)" -o "~{chromosome}.genotyped.vcf"
        bgzip -c "~{chromosome}.genotyped.vcf" > "~{chromosome}.genotyped.vcf.gz"
        tabix "~{chromosome}.genotyped.vcf.gz"
    >>>

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk 8000 SSD"
        memory: "12G"
        cpu : 4
    }

    output {
        File chromosomeVCF = "~{chromosome}.genotyped.vcf"
        File chromosomeVCF_gz = "~{chromosome}.genotyped.vcf.gz"
        File chromosomeVCF_tbi = "~{chromosome}.genotyped.vcf.gz.tbi"
    }
}