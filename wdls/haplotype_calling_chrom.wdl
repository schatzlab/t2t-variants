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
        File chr1_hcVCF_gz = hc1.hcVCF_gz
        File chr1_hcVCF_gz_tbi = hc1.hcVCF_gz_tbi
        File chr2_hcVCF = hc2.hcVCF
        File chr2_hcVCF_gz = hc2.hcVCF_gz
        File chr2_hcVCF_gz_tbi = hc2.hcVCF_gz_tbi
        File chr3_hcVCF = hc3.hcVCF
        File chr3_hcVCF_gz = hc3.hcVCF_gz
        File chr3_hcVCF_gz_tbi = hc3.hcVCF_gz_tbi
        File chr4_hcVCF = hc4.hcVCF
        File chr4_hcVCF_gz = hc4.hcVCF_gz
        File chr4_hcVCF_gz_tbi = hc4.hcVCF_gz_tbi
        File chr5_hcVCF = hc5.hcVCF
        File chr5_hcVCF_gz = hc5.hcVCF_gz
        File chr5_hcVCF_gz_tbi = hc5.hcVCF_gz_tbi
        File chr6_hcVCF = hc6.hcVCF
        File chr6_hcVCF_gz = hc6.hcVCF_gz
        File chr6_hcVCF_gz_tbi = hc6.hcVCF_gz_tbi
        File chr7_hcVCF = hc7.hcVCF
        File chr7_hcVCF_gz = hc7.hcVCF_gz
        File chr7_hcVCF_gz_tbi = hc7.hcVCF_gz_tbi
        File chr8_hcVCF = hc8.hcVCF
        File chr8_hcVCF_gz = hc8.hcVCF_gz
        File chr8_hcVCF_gz_tbi = hc8.hcVCF_gz_tbi
        File chr9_hcVCF = hc9.hcVCF
        File chr9_hcVCF_gz = hc9.hcVCF_gz
        File chr9_hcVCF_gz_tbi = hc9.hcVCF_gz_tbi
        File chr10_hcVCF = hc10.hcVCF
        File chr10_hcVCF_gz = hc10.hcVCF_gz
        File chr10_hcVCF_gz_tbi = hc10.hcVCF_gz_tbi
        File chr11_hcVCF = hc11.hcVCF
        File chr11_hcVCF_gz = hc11.hcVCF_gz
        File chr11_hcVCF_gz_tbi = hc11.hcVCF_gz_tbi
        File chr12_hcVCF = hc12.hcVCF
        File chr12_hcVCF_gz = hc12.hcVCF_gz
        File chr12_hcVCF_gz_tbi = hc12.hcVCF_gz_tbi
        File chr13_hcVCF = hc13.hcVCF
        File chr13_hcVCF_gz = hc13.hcVCF_gz
        File chr13_hcVCF_gz_tbi = hc13.hcVCF_gz_tbi
        File chr14_hcVCF = hc14.hcVCF
        File chr14_hcVCF_gz = hc14.hcVCF_gz
        File chr14_hcVCF_gz_tbi = hc14.hcVCF_gz_tbi
        File chr15_hcVCF = hc15.hcVCF
        File chr15_hcVCF_gz = hc15.hcVCF_gz
        File chr15_hcVCF_gz_tbi = hc15.hcVCF_gz_tbi
        File chr16_hcVCF = hc16.hcVCF
        File chr16_hcVCF_gz = hc16.hcVCF_gz
        File chr16_hcVCF_gz_tbi = hc16.hcVCF_gz_tbi
        File chr17_hcVCF = hc17.hcVCF
        File chr17_hcVCF_gz = hc17.hcVCF_gz
        File chr17_hcVCF_gz_tbi = hc17.hcVCF_gz_tbi
        File chr18_hcVCF = hc18.hcVCF
        File chr18_hcVCF_gz = hc18.hcVCF_gz
        File chr18_hcVCF_gz_tbi = hc18.hcVCF_gz_tbi
        File chr19_hcVCF = hc19.hcVCF
        File chr19_hcVCF_gz = hc19.hcVCF_gz
        File chr19_hcVCF_gz_tbi = hc19.hcVCF_gz_tbi
        File chr20_hcVCF = hc20.hcVCF
        File chr20_hcVCF_gz = hc20.hcVCF_gz
        File chr20_hcVCF_gz_tbi = hc20.hcVCF_gz_tbi
        File chr21_hcVCF = hc21.hcVCF
        File chr21_hcVCF_gz = hc21.hcVCF_gz
        File chr21_hcVCF_gz_tbi = hc21.hcVCF_gz_tbi
        File chr22_hcVCF = hc22.hcVCF
        File chr22_hcVCF_gz = hc22.hcVCF_gz
        File chr22_hcVCF_gz_tbi = hc22.hcVCF_gz_tbi
        File? XX_X_hcVCF = hc_XX_X.hcVCF
        File? XX_X_hcVCF_gz = hc_XX_X.hcVCF_gz
        File? XX_X_hcVCF_gz_tbi = hc_XX_X.hcVCF_gz_tbi
        File? XY_X_non_PAR_hcVCF = hc_XY_X_nonPAR.hcVCF
        File? XY_X_non_PAR_hcVCF_gz = hc_XY_X_nonPAR.hcVCF_gz
        File? XY_X_non_PAR_hcVCF_gz_tbi = hc_XY_X_nonPAR.hcVCF_gz_tbi
        File? XY_X_PAR_hcVCF = hc_XY_X_PAR.hcVCF
        File? XY_X_PAR_hcVCF_gz = hc_XY_X_PAR.hcVCF_gz
        File? XY_X_PAR_hcVCF_gz_tbi = hc_XY_X_PAR.hcVCF_gz_tbi
        File? XY_Y_nonPAR_hcVCF = hc_XY_Y_nonPAR.hcVCF
        File? XY_Y_nonPAR_hcVCF_gz = hc_XY_Y_nonPAR.hcVCF_gz
        File? XY_Y_nonPAR_hcVCF_gz_tbi = hc_XY_Y_nonPAR.hcVCF_gz_tbi
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

        bgzip -c "~{sampleName}.~{chrom}.hc.vcf" > "~{sampleName}.~{chrom}.hc.vcf.gz"
        tabix "~{sampleName}.~{chrom}.hc.vcf.gz"
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
        File hcVCF_gz = "~{sampleName}.~{chrom}.hc.vcf.gz"
        File hcVCF_gz_tbi = "~{sampleName}.~{chrom}.hc.vcf.gz.tbi"
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
        
        bgzip -c "~{sampleName}.~{regionName}.hc.vcf" > "~{sampleName}.~{regionName}.hc.vcf.gz"
        tabix "~{sampleName}.~{regionName}.hc.vcf.gz"
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
        File hcVCF_gz = "~{sampleName}.~{regionName}.hc.vcf.gz"
        File hcVCF_gz_tbi = "~{sampleName}.~{regionName}.hc.vcf.gz.tbi"
    }
}