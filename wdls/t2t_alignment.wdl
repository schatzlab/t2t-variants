version 1.0

workflow t2t_alignment {
    input {
        File targetRef
        String sampleName
        File inputFastq1
        File inputFastq2
        File bwaIndexTar
        Int dedupDistance = 100
    }

    call splitIntoLanes as splitFQ1 {
        input:
            fastq = inputFastq1,
            sampleName = sampleName
    }

    call splitIntoLanes as splitFQ2 {
        input:
            fastq = inputFastq2,
            sampleName = sampleName
    }

    Array[Pair[File, File]] matchedLanes = zip(splitFQ1.fastqs, splitFQ2.fastqs)
    
    scatter (pairedLanes in matchedLanes) {
        call alignLane {
            input:
                read1 = pairedLanes.left,
                read2 = pairedLanes.right,
                targetRef = targetRef,
                sampleName = sampleName,
                bwaIndexTar = bwaIndexTar
        }

        call fixmateLane {
            input:
                bam = alignLane.bam
        }

        call sortBamLane {
            input:
                bam = fixmateLane.fmBam
        }
    }

    call gatherMergeBam {
        input:
            laneBams = sortBamLane.sortBam,
            sampleName = sampleName
    }

    call markDuplicates {
        input:
            bam = gatherMergeBam.mergedBam,
            sampleName = sampleName,
            dedupDistance = dedupDistance
    }

    call bamtoCram {
        input:
            targetRef = targetRef,
            inputBam = markDuplicates.dedupBam,
            sampleName = sampleName
    }

    call samtoolsIndex as indexCRAM {
        input:
            alignmentFile = bamtoCram.cram,
            isBAM = false
    }

    call mosdepthStats {
        input:
            inputCram = bamtoCram.cram,
            cramIndex =  indexCRAM.alignmentIndex,
            targetRef = targetRef,
            sampleName = sampleName
    }

    call samtoolsStats {
        input:
            inputCram = bamtoCram.cram,
            cramIndex =  indexCRAM.alignmentIndex,
            targetRef = targetRef,
            sampleName = sampleName
    }

    output {
        File cram = bamtoCram.cram
        File cramIndex = indexCRAM.alignmentIndex
        File mosdepth_globalDist = mosdepthStats.globalDist
        File mosdepth_regionsDist = mosdepthStats.regionsDist
        File mosdepth_summary = mosdepthStats.summary
        File mosdepth_regionsBed = mosdepthStats.regionsBed
        File mosdepth_regionsBedIndex = mosdepthStats.regionsBedIndex
        File samtools_stats = samtoolsStats.stats
    }
}

task splitIntoLanes {
    input {
        File fastq
        String sampleName
    }

    command <<<
        echo "Number of lines in original FASTQ: $(zcat ~{fastq} | wc -l)"
        sample="~{sampleName}"

        fastq_name="~{fastq}"
        tmp="${fastq_name%.fastq.gz}"
        read="${tmp##*_}"

        bgzip -cd -@ 8 "~{fastq}" | sed 's/^@[^ ]* /@/g' | awk 'BEGIN {FS = ":"} {lane=$3"."$4 ; print $0 | "bgzip -@ 2 > sample_read_lane_"lane".fastq.gz"; for (i = 1; i <= 3; i++) {getline ; print $0 | "bgzip -@ 2 > sample_read_lane_"lane".fastq.gz"}}'

        for file in *.fastq.gz; do
            echo "Number of lines in $file: $(zcat $file | wc -l)"
            mv $file $(echo $file | sed "s/sample/$sample/; s/read/$read/")
        done

        rm "~{fastq}"
    >>>

    Int diskGb = ceil(10.0 * size(fastq, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        Array[File] fastqs = glob("*.fastq.gz")
    }
}

task alignLane {
    input {
        File read1
        File read2
        File targetRef
        String sampleName
        File bwaIndexTar
    }

    String bamBase='~{basename(read1,".fastq.gz")}'

    String fastaName='~{basename(targetRef)}'
    String read1Name='~{basename(read1)}'
    String read2Name='~{basename(read2)}'
    String tarName='~{basename(bwaIndexTar)}'

    command <<<
        fastqName="~{read1}"
        lane="${fastqName%.fastq.gz}"
        lane="${lane##*_lane}"

        bam_name="~{sampleName}.${lane}.bam"

        mv "~{targetRef}" .
        mv "~{read1}" .
        mv "~{read2}" .
        mv "~{bwaIndexTar}" .

        tar -xzf "~{tarName}"

        bwa mem -Y \
            -K 100000000 \
            -t "$(nproc)" \
            -R "@RG\tID:~{sampleName}_${lane}\tPL:illumina\tPM:Unknown\tLB:~{sampleName}\tDS:GRCh38\tSM:~{sampleName}\tCN:NYGenome\tPU:${lane}" \
            "./~{fastaName}" \
            "./~{read1Name}" \
            "./~{read2Name}" | samtools view -Shb -o "~{bamBase}.bam" -
    >>>

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk 1000 SSD"
        memory: "64G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File bam = "~{bamBase}.bam"
    }
}

task fixmateLane {
    input {
        File bam
    }

    String bamBase='~{basename(bam,".bam")}'

    command <<<
        samtools fixmate -m -@ "12" "~{bam}" "~{bamBase}.fm.bam"
    >>>

    Int diskGb = ceil(3.0 * size(bam, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File fmBam = "~{bamBase}.fm.bam"
    }
}

task sortBamLane {
    input {
        File bam
    }

    String bamBase='~{basename(bam,".bam")}'

    command <<<
        samtools sort -o "~{bamBase}.sort.bam" -@ 12 "~{bam}"
    >>>

    Int diskGb = ceil(5.0 * size(bam, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File sortBam = "~{bamBase}.sort.bam"
    }
}

task gatherMergeBam {
    input {
        Array[File] laneBams
        String sampleName
    }

    command <<<
        samtools merge -@ "$(nproc)" "~{sampleName}.merged.bam" ~{sep=' ' laneBams}
    >>>

    Int diskGb = ceil(3.0 * length(laneBams) * size(laneBams[0], "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk 1500 SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File mergedBam = "~{sampleName}.merged.bam"
    }
}

task markDuplicates {
    input {
        File bam
        String sampleName
        Int dedupDistance
    }

    command <<<
        samtools markdup -@ "$(nproc)" \
            -S -s \
            -f "~{sampleName}_metrics.txt" \
            -d ~{dedupDistance} \
            "~{bam}" \
            "~{sampleName}.dedup.bam"

        find / -name '*~{sampleName}*'
    >>>

    Int diskGb = ceil(5.0 * size(bam, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File dedupBam = "~{sampleName}.dedup.bam"
    }
}

task samtoolsIndex {
    input {
        File alignmentFile
        Boolean isBAM
    }

    String alignmentName = basename(alignmentFile)

    String outputSuffix = if (isBAM) then "bai" else "crai"

    command <<<
        samtools index -@ "$(nproc)" "~{alignmentFile}" "~{alignmentName}.~{outputSuffix}"
    >>>

    Int diskGb = ceil(2.0 * size(alignmentFile, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File alignmentIndex = "~{alignmentName}.~{outputSuffix}"
    }
}

task bamtoCram {
    input {
        File targetRef
        File inputBam
        String sampleName
    }

    command <<<
        samtools view -@ "$(nproc)" -C -T "~{targetRef}" -o "~{sampleName}.cram" "~{inputBam}"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File cram = "~{sampleName}.cram"
    }
}

task mosdepthStats {
    input {
        File inputCram
        File cramIndex
        File targetRef
        String sampleName
    }

    String cramName = basename(inputCram)

    command <<<
        mv "~{inputCram}" .
        mv "~{cramIndex}" .

        mosdepth -n --fast-mode --by 500 \
            -t "$(nproc)" \
            --fasta "~{targetRef}" \
            "~{sampleName}" \
            "~{cramName}"
    >>>

    Int diskGb = ceil(2.0 * size(inputCram, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File globalDist = "~{sampleName}.mosdepth.global.dist.txt"
        File regionsDist = "~{sampleName}.mosdepth.region.dist.txt"
        File summary = "~{sampleName}.mosdepth.summary.txt"
        File regionsBed = "~{sampleName}.regions.bed.gz"
        File regionsBedIndex = "~{sampleName}.regions.bed.gz.csi"
    }
}

task samtoolsStats {
    input {
        File inputCram
        File cramIndex
        File targetRef
        String sampleName
    }

    command <<<
        samtools stats -r "~{targetRef}" \
            --reference "~{targetRef}" \
            -@ "$(nproc)" \
            "~{inputCram}" > "~{sampleName}.samtools.stats.txt"
    >>>

    Int diskGb = ceil(2.0 * size(inputCram, "G"))

    runtime {
        docker : "szarate/t2t_variants:v0.0.2"
        disks : "local-disk ${diskGb} SSD"
        memory: "12G"
        cpu : 16
        preemptible: 3
        maxRetries: 3
    }

    output {
        File stats = "~{sampleName}.samtools.stats.txt"
    }
}