process Fraction {
    // Custom process to parse subsample fraction from Mosdepth summary file
    tag {"Fraction ${sample_id}"}
    label 'Fraction'
    shell = ['/bin/bash', '-eo', 'pipefail']

    input:
        tuple(val(sample_id), path(mosdepth_summary))

    output:
        tuple(val(sample_id), stdout)

    script:
        """
        cat ${mosdepth_summary} |tail -n1 |cut -f4 | awk '{if ($params.downsample_coverage/\$0 >1) {print 1} else {print $params.downsample_coverage/\$0}}' | tr -d '\n'
        """
}
