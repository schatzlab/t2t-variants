version 1.0

workflow t2t_haplotype_calling {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String sampleName
        String sex
        Boolean enableSpark
    }

    call haplotypeCaller as hc1 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr1",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc2 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr2",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc3 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr3",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc4 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr4",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc5 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr5",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc6 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr6",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc7 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr7",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc8 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr8",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc9 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr9",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc10 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr10",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc11 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr11",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc12 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr12",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc13 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr13",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc14 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr14",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc15 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr15",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc16 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr16",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc17 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr17",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc18 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr18",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc19 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr19",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc20 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr20",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    call haplotypeCaller as hc21 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr21",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }
    
    call haplotypeCaller as hc22 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr22",
            ploidy = 2,
            sampleName = sampleName,
            enableSpark = enableSpark
    }

    if (sex == "female") {
        call haplotypeCaller as hcFemale {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                chrom = "chrX",
                ploidy = 2,
                sampleName = sampleName,
                enableSpark = enableSpark
        }
    }

    if (sex == "male") {
        call haplotypeCaller as hcMaleX {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                chrom = "chrX",
                ploidy = 1,
                sampleName = sampleName,
                enableSpark = enableSpark
        }
        call haplotypeCaller as hcMaleY {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                chrom = "chrY",
                ploidy = 1,
                sampleName = sampleName,
                enableSpark = enableSpark
        }
    }

    output {
        File chr1_hcVCF = hc1.hcVCF
        File chr2_hcVCF = hc2.hcVCF
        File chr3_hcVCF = hc3.hcVCF
        File chr4_hcVCF = hc4.hcVCF
        File chr5_hcVCF = hc5.hcVCF
        File chr6_hcVCF = hc6.hcVCF
        File chr7_hcVCF = hc7.hcVCF
        File chr8_hcVCF = hc8.hcVCF
        File chr9_hcVCF = hc9.hcVCF
        File chr10_hcVCF = hc10.hcVCF
        File chr11_hcVCF = hc11.hcVCF
        File chr12_hcVCF = hc12.hcVCF
        File chr13_hcVCF = hc13.hcVCF
        File chr14_hcVCF = hc14.hcVCF
        File chr15_hcVCF = hc15.hcVCF
        File chr16_hcVCF = hc16.hcVCF
        File chr17_hcVCF = hc17.hcVCF
        File chr18_hcVCF = hc18.hcVCF
        File chr19_hcVCF = hc19.hcVCF
        File chr20_hcVCF = hc20.hcVCF
        File chr21_hcVCF = hc21.hcVCF
        File chr22_hcVCF = hc22.hcVCF
        File? female_hcVCF = hcFemale.hcVCF
        File? maleX_hcVCF = hcMaleX.hcVCF
        File? maleY_hcVCF = hcMaleY.hcVCF
    }
}

task haplotypeCaller {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String chrom
        String sampleName
        Int ploidy
        Boolean enableSpark
    }

    String tool = if (enableSpark) then "HaplotypeCallerSpark" else "HaplotypeCaller"
    String parallelism = if (enableSpark) then "--spark-runner LOCAL --spark-master $(nproc)" else ""

    String cramName = '~{basename(cram)}'
    String fastaName = '~{basename(refFasta)}'
    String indexName = '~{basename(fastaIndex)}'
    String dictName = '~{basename(fastaDict)}'

    command <<<
        mv "~{cram}" .
        mv "~{cramIndex}" .

        mv "~{refFasta}" .
        mv "~{fastaIndex}" .
        mv "~{fastaDict}" .

        gatk \
            "~{tool}" \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            -R "~{fastaName}" \
            -I ./"~{cramName}" \
            -L "~{chrom}" \
            -pairHMM AVX_LOGLESS_CACHING \
            -O "~{sampleName}.~{chrom}.hc.vcf" \
            -ERC GVCF \
            -ploidy "~{ploidy}" \
            -A Coverage \
            -A DepthPerAlleleBySample \
            -A DepthPerSampleHC \
            -A InbreedingCoeff \
            -A MappingQualityRankSumTest \
            -A MappingQualityZero \
            -A QualByDepth \
            -A ReadPosRankSumTest \
            -A RMSMappingQuality \
            -A StrandBiasBySample
    >>>

    Int diskGb = ceil(2.0 * size(cram, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "24G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File hcVCF = "~{sampleName}.~{chrom}.hc.vcf"
    }
}