process EditSummaryFileHappy {
    tag "$meta.id"
    label 'EditSummaryFileHappy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        // meta should have the keys 'id', 'query' and 'truth'
        tuple val(meta), path(summary_csv)

    output:
        path("*_INDEL_PASS.summary.csv"), emit: indel_pass_csv
        path("*_INDEL_ALL.summary.csv"), emit: indel_all_csv
        path("*_SNP_PASS.summary.csv"), emit: snp_pass_csv
        path("*_SNP_ALL.summary.csv"), emit: snp_all_csv

    script:
        """
        # Add samplenames as columns (header and row values) at start of line
        sed '1s/^/samples,sample_truth,sample_query,/; 2,\$s/^/${meta.truth}_${meta.query},${meta.truth},${meta.query},/' ${summary_csv} > ${summary_csv}.tmp
        
        # Split file including header (first line)
        awk -F',' 'FNR==1{hdr=\$0;next} {
            print hdr>"${meta.truth}_${meta.query}_"\$4"_"\$5".summary.csv";
            print \$0>>"${meta.truth}_${meta.query}_"\$4"_"\$5".summary.csv"
        }' ${summary_csv}.tmp
        
        # Remove tmp files
        rm ${summary_csv}.tmp     
        """
}
