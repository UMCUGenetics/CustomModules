process EditSummaryFileHappy {
    tag {"EditSummaryFileHappy"}
    label 'EditSummaryFileHappy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        // meta should have the keys 'id', 'query' and 'truth'
        tuple(val(meta), path(summary_csv))

    output:
        path("INDEL_PASS_summary.csv"), emit: indel_pass_summary_csv
        path("INDEL_ALL_summary.csv"), emit: indel_all_summary_csv
        path("SNP_PASS_summary.csv"), emit: snp_pass_summary_csv
        path("SNP_ALL_summary.csv"), emit: snp_all_summary_csv

    script:
        """
        # Add samplenames as columns (header and row values) at start of line
        sed '1s/^/samples,sample_truth,sample_query,/; 2,\$s/^/${meta.truth}_${meta.query},${meta.truth},${meta.query},/' ${summary_csv} > ${summary_csv}.tmp
        
        # Split file including header (first line)
        awk -F',' 'FNR==1{hdr=\$0;next} {print hdr>\$3"_"\$4"_summary.csv"; print \$0>>\$3"_"\$4"_summary.csv"}' ${summary_csv}.tmp
        
        # Remove tmp files
        rm ${summary_csv}.tmp     
        """
}