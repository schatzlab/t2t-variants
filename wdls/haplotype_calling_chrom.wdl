version 1.0

workflow haplotype_calling_chrom {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        File? parXbed
        File? nonparXbed
        File? parYbed
        String sampleName
        String chrXname = "chrX"
        String chrYname = "chrY"
        String sex
    }

    call diploidHC as hc1 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr1",
            sampleName = sampleName
    }

    call diploidHC as hc2 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr2",
            sampleName = sampleName
    }

    call diploidHC as hc3 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr3",
            sampleName = sampleName
    }

    call diploidHC as hc4 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr4",
            sampleName = sampleName
    }

    call diploidHC as hc5 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr5",
            sampleName = sampleName
    }

    call diploidHC as hc6 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr6",
            sampleName = sampleName
    }

    call diploidHC as hc7 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr7",
            sampleName = sampleName
    }

    call diploidHC as hc8 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr8",
            sampleName = sampleName
    }

    call diploidHC as hc9 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr9",
            sampleName = sampleName
    }

    call diploidHC as hc10 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr10",
            sampleName = sampleName
    }

    call diploidHC as hc11 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr11",
            sampleName = sampleName
    }

    call diploidHC as hc12 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr12",
            sampleName = sampleName
    }

    call diploidHC as hc13 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr13",
            sampleName = sampleName
    }

    call diploidHC as hc14 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr14",
            sampleName = sampleName
    }

    call diploidHC as hc15 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr15",
            sampleName = sampleName
    }

    call diploidHC as hc16 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr16",
            sampleName = sampleName
    }

    call diploidHC as hc17 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr17",
            sampleName = sampleName
    }

    call diploidHC as hc18 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr18",
            sampleName = sampleName
    }

    call diploidHC as hc19 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr19",
            sampleName = sampleName
    }

    call diploidHC as hc20 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr20",
            sampleName = sampleName
    }

    call diploidHC as hc21 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr21",
            sampleName = sampleName
    }
    
    call diploidHC as hc22 {
        input:
            refFasta = refFasta,
            fastaIndex = fastaIndex,
            fastaDict = fastaDict,
            cram = cram,
            cramIndex = cramIndex,
            chrom = "chr22",
            sampleName = sampleName
    }

    if (sex == "female") {
        call diploidHC as hc_XX_X {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                chrom = chrXname,
                sampleName = sampleName
        }
    }

    if (sex == "male") {
        call sexchrHC as hc_XY_X_nonPAR {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                region = chrXname,
                regionName = "chrX_non_PAR",
                excludeInterval = parXbed,
                ploidy = 1,
                sampleName = sampleName
        }
        call sexchrHC as hc_XY_X_PAR {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                region = chrXname,
                regionName = "chrX_PAR",
                excludeInterval = nonparXbed,
                ploidy = 2,
                sampleName = sampleName
        }
        call sexchrHC as hc_XY_Y_nonPAR {
            input:
                refFasta = refFasta,
                fastaIndex = fastaIndex,
                fastaDict = fastaDict,
                cram = cram,
                cramIndex = cramIndex,
                region = chrYname,
                regionName = "chrY",
                excludeInterval = parYbed,
                ploidy = 1,
                sampleName = sampleName
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
        File? XX_X_hcVCF = hc_XX_X.hcVCF
        File? XY_X_non_PAR_hcVCF = hc_XY_X_nonPAR.hcVCF
        File? XY_X_PAR_hcVCF = hc_XY_X_PAR.hcVCF
        File? XY_Y_nonPAR_hcVCF = hc_XY_Y_nonPAR.hcVCF
    }
}

task diploidHC {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String chrom
        String sampleName
    }

    String cramName = '~{basename(cram)}'
    String fastaName = '~{basename(refFasta)}'
    String refBaseFasta = '~{basename(refFasta, ".fasta")}'
    String refBase = '~{basename(refBaseFasta, ".fa")}'

    command <<<
        mkdir -p inputs/
        mv "~{cram}" inputs/
        mv "~{cramIndex}" inputs/

        mv "~{refFasta}" inputs/
        mv "~{fastaIndex}" inputs/
        mv "~{fastaDict}" "inputs/~{refBase}.dict"

        gatk HaplotypeCaller \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            -R inputs/"~{fastaName}" \
            -I inputs/"~{cramName}" \
            -L "~{chrom}" \
            -pairHMM AVX_LOGLESS_CACHING \
            -O "~{sampleName}.~{chrom}.hc.vcf" \
            -ERC GVCF \
            -ploidy 2 \
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
        docker : "szarate/t2t_variants:v0.0.2"
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

task sexchrHC {
    input {
        File refFasta
        File fastaIndex
        File fastaDict
        File cram
        File cramIndex
        String region
        File? regionBed
        String regionName
        File? excludeInterval
        String sampleName
        Int ploidy
    }

    String excludeName = '~{(if defined(excludeInterval) then basename(select_first([excludeInterval])) else "")}'
    String regions = '~{(if defined(regionBed) then basename(select_first([regionBed])) else "")}'

    String cramName = '~{basename(cram)}'
    String fastaName = '~{basename(refFasta)}'
    String refBaseFasta = '~{basename(refFasta, ".fasta")}'
    String refBase = '~{basename(refBaseFasta, ".fa")}'

    command <<<
        mkdir -p inputs/

        mv "~{cram}" inputs/
        mv "~{cramIndex}" inputs/

        mv "~{refFasta}" inputs/
        mv "~{fastaIndex}" inputs/
        mv "~{fastaDict}" "inputs/~{refBase}.dict"

        echo "~{excludeName}"
        echo "~{regions}"

        if [[ -n "~{regions}" ]]; then
            mv "~{regionBed}" inputs/
            specifyRegion="--intervals inputs/~{regions}"
        else
            specifyRegion="--intervals ~{region}"
        fi

        if [[ -n "~{excludeName}" ]]; then
            mv "~{excludeInterval}" inputs/
            specifyExclude="--exclude-intervals inputs/~{excludeName}"
        else
            specifyExclude=""
        fi

        gatk HaplotypeCaller \
            --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=$(nproc) -Djava.io.tmpdir=/dev/shm" \
            -R inputs/"~{fastaName}" \
            -I inputs/"~{cramName}" \
            ${specifyRegion} ${specifyExclude} \
            -pairHMM AVX_LOGLESS_CACHING \
            -O "~{sampleName}.~{regionName}.hc.vcf" \
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
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "24G"
        cpu : 2
        preemptible: 3
        maxRetries: 3
    }

    output {
        File hcVCF = "~{sampleName}.~{regionName}.hc.vcf"
    }
}