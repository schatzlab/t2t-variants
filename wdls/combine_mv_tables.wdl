version 1.0

workflow combine_mv_tables {
    input {
        File concatMVtable
        File mvTableHeader
    }

    call combine_mv_tables {
        input:
            concatMVtable = concatMVtable,
            mvTableHeader = mvTableHeader
    }

    output {
        File mvTable = combine_mv_tables.table
    }
}

task combine_mv_tables {
    input {
        File concatMVtable
        File mvTableHeader
    }

    String tablePrefix = '~{basename(concatMVtable,".txt")}'

    command <<<
        # 5: nVar, 6: nSkipped, 9: nLowQual, 12: nViolations
        # can add more fields later with header if desired; just keeping these for now
        # echo -e "Family\tnVariants\tnSkipped\tnLowQual\tnViolations" > "~{tablePrefix}.table"
        awk -v OFS='\t' '{a[$4]+=$5;b[$4]+=$6;c[$4]+=$9;d[$4]+=$12}END{for(i in a) print i,a[i],b[i],c[i],d[i]}' "~{concatMVtable}" | grep -v 'all' | grep -v 'Family' >> "~{tablePrefix}.txt"
        grep -v "	0	0	0	0" "~{tablePrefix}.txt" > "~{tablePrefix}.table.txt"
    >>>

    Int diskGb = ceil(2.0 * size(concatMVtable, "G"))

    runtime {
        docker : "szarate/t2t_variants"
        disks : "local-disk ${diskGb} SSD"
        memory: "4G"
        cpu : 2
    }

    output {
        File table = "~{tablePrefix}.table.txt"
    }
}
