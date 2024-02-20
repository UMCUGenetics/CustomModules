process EditSummaryFileHappy {
    tag {"EditSummaryFileHappy"}
    label 'EditSummaryFileHappy'
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        tuple(val(meta), path(summary_csv))

    output:
        path("INDEL_PASS_summary.txt", includeInputs: false), emit: indel_pass_summary_txt
        path("INDEL_ALL_summary.txt", includeInputs: false), emit: indel_all_summary_txt
        path("SNP_PASS_summary.txt", includeInputs: false), emit: snp_pass_summary_txt
        path("SNP_ALL_summary.txt", includeInputs: false), emit: snpp_all_summary_txt

    script:
        """
        # Add samplenames as columns
        sed '1s/$/,sample_query,sample_truth/; 2,$s/$/,${meta.query},${meta.truth}' ${summary_csv} > ${summary_csv}.tmp
        
        # Split file and add column
        awk '{print $0 > $3_$4"_summary.txt"}' ${summary_csv}.tmp
        
        # Remove tmp files
        rm *tmp     
        """
}