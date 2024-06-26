process SampleUDFDx {
    // Custom process to run clarity_epp export sample_udf_dx
    tag {"ClarityEpp SampleUDFDx ${sample_id}"}
    label 'ClarityEpp'
    label 'ClarityEpp_SampleUDFDx'
    shell = ['/bin/bash', '-eo', 'pipefail']
    cache = false  //Disable cache to force a clarity export restarting the workflow.

    input:
        val(sample_id)

    output:
        tuple(val(sample_id), stdout)

    script:
        """
        source ${params.clarity_epp_path}/venv/bin/activate
        python ${params.clarity_epp_path}/clarity_epp.py export sample_udf_dx \
        -a ${sample_id} -u '$params.udf' -c '$params.column_name' | cut -f 2 | grep -v $params.column_name | tr -d '\n'
        """
}
