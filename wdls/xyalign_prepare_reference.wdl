version 1.0

workflow xyalign_prepare_reference {
    input {
        File refFasta
        File refMask
        String chrXname
        String chrYname
    }

    call prepare_reference {
        input:
            refFasta = refFasta,
            refMask = refMask,
            chrXname = chrXname,
            chrYname = chrYname
    }

    output {
        File XXref = prepare_reference.XXref
        File XXrefIndex = prepare_reference.XXrefIndex
        File XXrefDict = prepare_reference.XXrefDict
        File XYref = prepare_reference.XYref
        File XYrefIndex = prepare_reference.XYrefIndex
        File XYrefDict = prepare_reference.XYrefDict
        File out_tar = prepare_reference.out_tar
    }
}

task prepare_reference {
    input {
        File refFasta
        File refMask
        String chrXname
        String chrYname
    }

    String refBase = '~{basename(refFasta,".fasta")}'
    String ref = '~{basename(refBase,".fa")}'

    command <<<
        mkdir out/

        xyalign --PREPARE_REFERENCE --ref "~{refFasta}" --xx_ref_out "~{ref}.XX.fasta" --xy_ref_out "~{ref}.XY.fasta" --x_chromosome "~{chrXname}" --y_chromosome "~{chrYname}" --reference_mask "~{refMask}" --output_dir "out"
        
        tar -czvf prepare_ref.out.tar.gz out/
    >>>

    Int diskGb = ceil(4.0 * size(refFasta, "G"))

    runtime {
        docker : "szarate/xyalign"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File XXref = "~{refBase}.XX.fasta"
        File XXrefIndex = "~{refBase}.XX.fasta.fai"
        File XXrefDict = "~{refBase}.XX.fasta.dict"
        File XYref = "~{refBase}.XY.fasta"
        File XYrefIndex = "~{refBase}.XY.fasta.fai"
        File XYrefDict = "~{refBase}.XY.fasta.dict"
        File out_tar = "prepare_ref.out.tar.gz"
    }
}